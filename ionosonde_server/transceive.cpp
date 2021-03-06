// tx_worker.cpp
// 
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/stream.hpp>

#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <fstream>
#include <csignal>
#include <thread>
#include <mutex>

#include "global_variables.h"
#include "sounder.hpp"

void tx_worker(
        unsigned int bufflen,
        const uhd::tx_streamer::sptr & tx_stream,
        uhd::time_spec_t start_time,
        std::complex<int16_t> *vec_ptr,
        int end
);

void rx_worker(
        const uhd::rx_streamer::sptr & rx_stream,
        unsigned int samps_per_pulse,
        std::vector<std::complex<int16_t> *> &recv_ptr
);

void transceive(
        uhd::usrp::multi_usrp::sptr usrp,
        uhd::tx_streamer::sptr tx_stream,
        uhd::rx_streamer::sptr rx_stream,
        unsigned int npulses,
        float pulse_time,
        std::vector<std::complex<int16_t> > *txbuff0,
        std::vector<std::complex<int16_t> > *txbuff1,
        float tx_ontime,
        std::complex<int16_t> **outdata,
        size_t samps_per_pulse) {

    BOOST_LOG_TRIVIAL(trace) << "Entered function transceive";
    BOOST_LOG_TRIVIAL(debug) << "Samples per pulse: " << samps_per_pulse;

    uhd::time_spec_t start_time = usrp->get_time_now() + 0.05;
    BOOST_LOG_TRIVIAL(info) << "Time start: " << start_time.get_full_secs() << start_time.get_frac_secs();

    usrp->set_gpio_attr("RXA", "CTRL", 0x0, TRTRIG_BIT); //GPIO mode
    usrp->set_gpio_attr("RXA", "DDR", TRTRIG_BIT, TRTRIG_BIT); //Direction out


    std::vector<std::complex<int16_t> > buff(samps_per_pulse, 0);
    BOOST_LOG_TRIVIAL(info) << "Receive buffer size:  " << buff.size();
    BOOST_LOG_TRIVIAL(info) << "Transmit buffer size: " << txbuff0->size();
    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = npulses * samps_per_pulse;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec =
            start_time + 22 / usrp->get_rx_rate(); //Digital hardware delay is 22 samples long.  Found by experiment.

    //loop for every pulse in the sequence
    std::vector<std::complex<int16_t> *> rx_dptr;
    rx_dptr.resize(usrp->get_rx_num_channels());
    BOOST_LOG_TRIVIAL(info) << "npulses: " << npulses;
    BOOST_LOG_TRIVIAL(info) << "pulse_time: " << pulse_time;

    boost::thread_group rx_threads;
    boost::thread_group tx_threads;
    for (int ipulse = 0; ipulse < npulses; ipulse++) {
        BOOST_LOG_TRIVIAL(trace) << "pulse number: " << ipulse;

        for (size_t ichan = 0; ichan < usrp->get_rx_num_channels(); ichan++) {
            rx_dptr[ichan] = ipulse * samps_per_pulse + outdata[ichan];
        }

        usrp->set_command_time(start_time - 50e-6, 0);
        usrp->set_gpio_attr("RXA", "OUT", T_BIT, 0x8100);

        if (ipulse == 0) {
            BOOST_LOG_TRIVIAL(info) << "time spec: " << stream_cmd.time_spec.get_real_secs();
            BOOST_LOG_TRIVIAL(info) << "Issuing stream command to start collecting samples";
            usrp->issue_stream_cmd(stream_cmd);
        }

        usrp->set_command_time(start_time + tx_ontime, 0);
        usrp->set_gpio_attr("RXA", "OUT", R_BIT, 0x8100);


        // Segfault sometime after this
        std::vector<std::complex<int16_t> *> vec_ptr;
        vec_ptr.resize(1);

        if (ipulse % 2 == 0) {
            vec_ptr[0] = &txbuff0->front();
        } else  {
            vec_ptr[0] = &txbuff1->front();
        }

        if (ipulse != npulses - 1) {
            tx_threads.create_thread(boost::bind(tx_worker,
                                                 txbuff0->size(), tx_stream, start_time, vec_ptr[0], 0));
        } else {
            tx_threads.create_thread(boost::bind(tx_worker,
                                                 txbuff0->size(), tx_stream, start_time, vec_ptr[0], 1));
        }

        rx_threads.join_all();
        rx_threads.create_thread(boost::bind(rx_worker,
                                             rx_stream, samps_per_pulse, rx_dptr));


        start_time += pulse_time;
    }
    tx_threads.join_all();
    rx_threads.join_all();
}

/***********************************************************************
 * tx_worker function
 * A function to be used in its own thread for transmitting.  Push all
 * tx values into the USRP buffer as USRP buffer space is available,
 * but allow other actions to occur concurrently.
 **********************************************************************/
 // TODO: make send call blocking
 std::mutex tx_send_mutex;
void tx_worker(
        const unsigned int bufflen,
        const uhd::tx_streamer::sptr & tx_stream,
        const uhd::time_spec_t start_time,
        std::complex<int16_t> *vec_ptr,
        const int end
) {
    BOOST_LOG_TRIVIAL(trace) << "Entered function tx_worker";
    unsigned int acc_samps = 0;

    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.has_time_spec = true;
    md.time_spec = start_time;

    size_t spb = tx_stream->get_max_num_samps();
    if (spb > bufflen) spb = bufflen;


    std::lock_guard<std::mutex> lock(tx_send_mutex);
    while (acc_samps < bufflen - spb) {
        size_t nsamples = tx_stream->send(vec_ptr, spb, md);
        vec_ptr += spb;
        acc_samps += nsamples;
        md.start_of_burst = false;
        md.has_time_spec = false;
        //BOOST_LOG_TRIVIAL(trace) << "Sent " << acc_samps << " tx packets";
    }
    // Now on the last packet
    if (end) md.end_of_burst = true;
    spb = bufflen - acc_samps;
    size_t nsamples = tx_stream->send(vec_ptr, spb, md);
    BOOST_LOG_TRIVIAL(trace) << "Sent last tx packet";
}

void rx_worker(
        const uhd::rx_streamer::sptr & rx_stream,
        const unsigned int samps_per_pulse,
        std::vector<std::complex<int16_t> *> &recv_ptr
) {
    BOOST_LOG_TRIVIAL(trace) << "Entered function rx_worker";
    uhd::rx_metadata_t rxmd;
    float timeout = 4.1;
    rxmd.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;
    size_t nrx_samples = rx_stream->recv(recv_ptr, samps_per_pulse, rxmd, timeout);
    BOOST_LOG_TRIVIAL(trace) << "Received " << nrx_samples << " packets";
    if (rxmd.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
        BOOST_LOG_TRIVIAL(error) << rxmd.error_code;
    }
}
