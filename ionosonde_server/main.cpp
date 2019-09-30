// main.cpp
// 
// Program for USRP radar operations
// Adapted from example txrx_loopback_to_file.cpp
// Alex Morris
// 06 Dec 2013
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>

#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>


#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <csignal>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <unistd.h>
#include <fftw3.h>

#include "sounder.hpp"
#include "lp_filter.hpp"
#include "matched_filter.hpp"

//#include "c_utils.h"
#include "utils.hpp"
#include "global_variables.h"

namespace po = boost::program_options;
typedef std::complex <int16_t> sc16;

namespace logging = boost::log;

void init_logging()
{
    logging::core::get()->set_filter
            (
                    logging::trivial::severity >= logging::trivial::info
            );
}

/***********************************************************************
 * Signal interrupt handlers
 **********************************************************************/
static bool stop_signal_called = false;

void sig_int_handler(int) { stop_signal_called = true; }

int verbose = 1;
int debug = 0;

void choose_pcode(std::vector<float> *pcode0, std::vector<float> *pcode1, size_t max_code_length) {


    if (max_code_length < 4) {
        *pcode0 = RECT;
        *pcode1 = RECT;
    } else if (max_code_length < 8) {
        *pcode0 = GOLAY_4_0;
        *pcode1 = GOLAY_4_1;
    } else if (max_code_length < 10) {
        *pcode0 = GOLAY_8_0;
        *pcode1 = GOLAY_8_1;
    } else if (max_code_length < 16) {
        *pcode0 = GOLAY_10_0;
        *pcode1 = GOLAY_10_1;
    } else {
        *pcode0 = GOLAY_16_0;
        *pcode1 = GOLAY_16_1;
    }


}

void set_frontend_parms(int freq, uhd::usrp::multi_usrp::sptr usrp) {
    uint32_t ctrl_bits = 0x00;

    /*Configure control bits for low pass filter*/
    if (freq > 12e6) {
        ctrl_bits |= LPF_32; //Use low pass cutoff of 32 MHz
    } else if (freq < 12e6 && freq > 6e6) {
        ctrl_bits |= LPF_16; //Use low pass cutoff of 16 MHz
    } else if (freq < 6e6 && freq > 3e6) {
        ctrl_bits |= LPF_8; //Use low pass cutoff of 8 MHz
    } else if (freq < 3e6) {
        ctrl_bits |= LPF_4; //Use low pass cutoff of 4 MHz
    }

    /*Configure control bits for hi pass filter*/
    if (freq < 2e6) {
        ctrl_bits |= HPF_1; //Use hipass cutoff of 1 MHz
    } else if (freq > 2e6 && freq < 6e6) {
        ctrl_bits |= HPF_2; //Use hi pass cutoff of 2 MHz
    } else if (freq > 6e6 && freq > 12e6) {
        ctrl_bits |= HPF_4; //Use hi pass cutoff of 4 MHz
    } else if (freq > 12e6) {
        ctrl_bits |= HPF_8; //Use hi pass cutoff of 8 MHz
    }

    /*Configure attenuation bits*/
    // ..?

    /* Send bits to USRP after configuring GPIO*/
    usrp->set_gpio_attr("RXA", "CTRL", 0x0, 0xfc); //GPIO mode
    usrp->set_gpio_attr("RXA", "DDR", 0xfc, 0xfc); //Direction out
    usrp->set_gpio_attr("RXA", "OUT", 0x0, 0xfc);
    usrp->set_gpio_attr("RXA", "OUT", ctrl_bits, 0xfc);

    /*Configure TX front-end to use the proper output port*/
    if (freq < XOVER_FREQ) {
        if (verbose) std::cout << "Using antenna A\n";
        usrp->set_tx_subdev_spec(std::string("A:A"));
    } else {
        if (verbose) std::cout << "Using antenna B\n";
        usrp->set_tx_subdev_spec(std::string("A:B"));
    }

}

/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char *argv[]) {

    init_logging();

    BOOST_LOG_TRIVIAL(info) << "Starting Server";

    std::ofstream myfile;
    myfile.open("example.txt");

    std::ofstream perfile;
    perfile.open("periodogram.txt");

    uhd::set_thread_priority_safe();

    int return_status = 0;

    //universal variables to be set by po
    std::string ref, otw, type;
    double freq;

    std::vector<float> pcode0;
    std::vector<float> pcode1;
    std::vector<float> fake_pcode0;
    std::vector<float> fake_pcode1;
    float *pcode_ptrs[2];
    float *fake_pcode_ptrs[2];

    //status flags
    int new_seq_flag = 1;

    //periodogram variables
    size_t center_freq_khz;
    size_t span_khz;
    std::vector<float> periodogram;

    //transmit variables
    std::string tx_subdev = "A:A";
    unsigned int bufflen;
    unsigned int samps_per_sym;
    unsigned int symboltime_usec, ipp_usec, max_code_length;
    boost::thread_group transmit_thread;
    std::vector <std::complex<float>> tx_raw_buff0;
    std::vector <std::complex<float>> tx_raw_buff1;
    std::vector <std::complex<int16_t> > tx_filt_buff0;
    std::vector <std::complex<int16_t> > tx_filt_buff1;
    std::vector <std::complex<float>> filter_taps;
    int ntaps;
    float txsamprate, txbw, tx_ontime;
    size_t tx_ontime_usec;

    //receive variables
    std::string rx_subdev = "A:A A:B";
    size_t spb;
    double rx_rate;
    unsigned int nave = 1;
    std::vector < std::vector < std::complex < int16_t > > > rawvecs;
    std::vector < std::complex < int16_t > * > rawvec_ptrs;
    float ptime_eff;
    unsigned int nsamps_per_pulse;

    //data processing variables
    int dmrate, slowdim, fastdim;
    size_t nranges, first_range_km, last_range_km;
    size_t first_inx;
    size_t filter_delay;
    float bandwidth;
    size_t osr;

    std::vector < std::vector < std::complex < float > > > outvecs[2][2];
    std::vector < std::complex < float > * > outvec_ptrs[2][2];

    std::vector < std::vector < std::complex < float > > > filtvecs[2][2];
    std::vector < std::complex < float > * > filtvec_ptrs[2][2];
    std::complex<float> **filtvec_dptr[2][2];

    std::vector < std::vector < std::complex < float > > > ffvecs[2];
    std::vector < std::complex < float > * > ffvec_ptrs[2];
    std::vector<float> fpow[2];
    std::vector<float> fvel[2];

    //socket-related variables
    int sock, msgsock, rval, rfds, efds;
    int msg;
    struct soundingParms2 parms, actual_parms;
    struct periodogramParms lparms;
    char usrpmsg;

    //create a usrp device
    std::string args = "addr=192.168.10.2";
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
    BOOST_LOG_TRIVIAL(info) << "Created USRP object at IP address: " << args;

    //Lock device to the motherboard clock and 
    //then initialize to an arbitrary time
    usrp->set_clock_source("internal");
    usrp->set_time_now(uhd::time_spec_t(0.0));

    //create a transmit streamer
    uhd::stream_args_t tx_stream_args("sc16", "sc16");
    usrp->set_tx_subdev_spec(tx_subdev);
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(tx_stream_args);

    usrp->set_rx_subdev_spec(rx_subdev);

    //create a receive streamer
    BOOST_LOG_TRIVIAL(info) << "Number of channels: " << usrp->get_rx_num_channels();

    uhd::stream_args_t rx_stream_args("sc16", "sc16");
    for (size_t rx_chan = 0; rx_chan < usrp->get_rx_num_channels(); rx_chan++) {
        rx_stream_args.channels.push_back(rx_chan); //linear mapping
    }

    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(rx_stream_args);

    BOOST_LOG_TRIVIAL(info) << boost::format("Using Device: %s") % usrp->get_pp_string();

    //USRP is initialized;
    //now execute the tx/rx operations per arguments passed by the tcp socket
    sock = tcpsocket(HOST_PORT);

    BOOST_LOG_TRIVIAL(info) << boost::format("socket: %i\n") % sock << std::endl;

    while (true) {
        listen(sock, 1);

        BOOST_LOG_TRIVIAL(info) <<"Waiting for client connection...";

        msgsock = accept(sock, nullptr, nullptr);

        BOOST_LOG_TRIVIAL(info) << boost::format("Listening on socket %i\n") % msgsock;

        while (true) {
            rval = recv_data(msgsock, &usrpmsg, sizeof(usrpmsg));

            if (rval <= 0) BOOST_LOG_TRIVIAL(info) << "breaking...";

            if (rval <= 0) myfile.close();
            if (rval <= 0) perfile.close();
            if (rval <= 0) break;

            switch (usrpmsg) {
                case EXIT:
                    BOOST_LOG_TRIVIAL(info) << boost::format("Client done");
                    break;

                case LISTEN:
                    BOOST_LOG_TRIVIAL(info) << boost::format("Starting Listening.");
                    if (verbose > -1) perfile << boost::format("Starting Listening.") << std::endl;

                    rval = recv_data(msgsock, &lparms, sizeof(lparms));
                    center_freq_khz = (lparms.end_freq_khz + lparms.start_freq_khz) / 2;
                    span_khz = lparms.end_freq_khz - lparms.start_freq_khz;
                    if (verbose) {
                        perfile << boost::format("start freq: ") << lparms.start_freq_khz << std::endl;
                        perfile << boost::format("end freq: ") << lparms.end_freq_khz << std::endl;
                        perfile << boost::format("bandwidth: ") << lparms.bandwidth_khz << std::endl;
                        perfile << boost::format("center_freq_khz: ") << center_freq_khz << std::endl;
                        perfile << boost::format("span_khz: ") << span_khz << std::endl;
                    }
                    BOOST_LOG_TRIVIAL(info) << boost::format("Start frequency (khz): ") << lparms.start_freq_khz;
                    BOOST_LOG_TRIVIAL(info) << boost::format("End frequency (khz): ") << lparms.end_freq_khz;
                    BOOST_LOG_TRIVIAL(info) << boost::format("Bandwidth (khz): ") << lparms.bandwidth_khz;
                    BOOST_LOG_TRIVIAL(info) << boost::format("Center frequency (khz): ") << center_freq_khz;
                    BOOST_LOG_TRIVIAL(info) << boost::format("Span (khz): ") << span_khz;

                    periodogram.resize(span_khz);

                    capture_spectrum(
                            usrp,
                            rx_stream,
                            lparms.start_freq_khz,
                            lparms.end_freq_khz,
                            lparms.bandwidth_khz,
                            &periodogram.front());

                    if (verbose > -1) {
                        for (int k = 0; k < span_khz; k++) {
                            perfile << k << "\t";
                            perfile << periodogram[k] << std::endl;
                        }

                    }

                    send(msgsock, &periodogram.front(), periodogram.size() * sizeof(float), 0);
                    return_status = 0;
                    send(msgsock, &return_status, sizeof(return_status), 0);
                    break;

                case SEND:

                    // 170525 GLB trying to determine why needs server reset
                    //usrp->set_time_now(uhd::time_spec_t(0.0));

                    BOOST_LOG_TRIVIAL(info) << "Starting sounding.";

                    rval = recv_data(msgsock, &parms, sizeof(parms));

                    osr = parms.over_sample_rate;

                    symboltime_usec = parms.range_res_km / 1.5e-1;
                    dmrate = (size_t) ceil(symboltime_usec * RX_RATE / (osr * 1e6));
                    // Decimation rate?
                    BOOST_LOG_TRIVIAL(info) << boost::format("dmrate: ") << dmrate << std::endl;

                    symboltime_usec = osr * 1e6 * dmrate / RX_RATE;

                    ipp_usec = (size_t) ceil(2 * parms.last_range_km / 3.0e-1 / symboltime_usec / 100);
                    ipp_usec *= 100;
                    ipp_usec *= symboltime_usec;

                    if (verbose)
                        myfile << boost::format("Symboltime in microseconds: ") << symboltime_usec << std::endl;
                    BOOST_LOG_TRIVIAL(info) << boost::format("Symboltime in microseconds: ") << symboltime_usec;
                    if (verbose) myfile << boost::format("Interpulse period: ") << ipp_usec << std::endl;
                    BOOST_LOG_TRIVIAL(info) << boost::format("Interpulse period: ") << ipp_usec;

                    nsamps_per_pulse = (unsigned int) (ipp_usec * RX_RATE / 1e6);

                    if (verbose) myfile << boost::format("Samples per pulse: ") << nsamps_per_pulse << std::endl;
                    BOOST_LOG_TRIVIAL(info) << boost::format("Samples per pulse: ") << nsamps_per_pulse;

                    max_code_length = (size_t) floor(2 * parms.first_range_km / 3.0e-1 / symboltime_usec);

                    if (verbose) myfile << boost::format("Max code length: ") << max_code_length << std::endl;
                    BOOST_LOG_TRIVIAL(info) << boost::format("Max code length: ") << max_code_length;

                    choose_pcode(&pcode0, &pcode1, max_code_length);


                    for (int i = 0; i < pcode0.size(); i++) {
                        if (verbose) myfile << pcode0[i] << " " << pcode1[i] << std::endl;
                        BOOST_LOG_TRIVIAL(info) << pcode0[i] << " " << pcode1[i];
                    }

                    freq = 1e3 * parms.freq_khz;
                    set_frontend_parms(freq, usrp);

                    //configure the USRP according to the arguments from the client
                    for (size_t i = 0; i < usrp->get_rx_num_channels(); i++) {
                        usrp->set_rx_freq(freq, i);
                    }
                    usrp->set_tx_freq(freq);
                    usrp->set_rx_rate(RX_RATE);
                    usrp->set_tx_rate(TX_RATE);


                    BOOST_LOG_TRIVIAL(info) << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq() / 1e6);

                    BOOST_LOG_TRIVIAL(info) << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq() / 1e6);

                    if (new_seq_flag == 1) {
                        usrp->set_tx_rate(TX_RATE);
                        BOOST_LOG_TRIVIAL(info) << boost::format("Actual TX Rate: %f Ksps...") % (usrp->get_tx_rate() / 1e3);
                        BOOST_LOG_TRIVIAL(info) << boost::format("Actual TX Rate: %f Ksps...") % (usrp->get_tx_rate() / 1e3);
                        usrp->set_rx_rate(RX_RATE);
                        BOOST_LOG_TRIVIAL(info) << boost::format("Actual RX Rate: %f Ksps...") % (usrp->get_rx_rate() / 1e3);
                        BOOST_LOG_TRIVIAL(info) << boost::format("Actual RX Rate: %f Ksps...") % (usrp->get_rx_rate() / 1e3);
                    }
                    new_seq_flag = 0;

                    stop_signal_called = false;

                    //Prepare lp filter taps for filtering the tx samples
                    samps_per_sym = symboltime_usec * TX_RATE / 1e6;
                    ntaps = 4 * samps_per_sym;
                    filter_taps.resize(ntaps, 0.0);
                    txbw = 1 / (2.e-6 * symboltime_usec);
                    txsamprate = usrp->get_tx_rate();
                    for (int i = 0; i < ntaps; i++) {
                        double x = 2 * (2 * M_PI * ((float) i / ntaps) - M_PI);
                        filter_taps[i] = std::complex<float>(1 * (0.54 - 0.46
                                                                         * cos((2 * M_PI * ((float) (i) + 0.5)) /
                                                                               ntaps)) * sin(x) / (x), 0);
                    }
                    filter_taps[ntaps / 2] = std::complex<float>(1, 0);

                    if (debug) {
                        for (int i = 0; i < ntaps; i++) {
                            printf("filter_taps %i: %f, %f\n", i, filter_taps[i].real(), filter_taps[i].imag());
                        }
                    }

                    //prepare raw tx information
                    pcode_ptrs[0] = &pcode0.front();
                    pcode_ptrs[1] = &pcode1.front();
                    fake_pcode0 = pcode0;
                    fake_pcode1 = pcode1;
                    fake_pcode_ptrs[0] = &fake_pcode0.front();
                    fake_pcode_ptrs[1] = &fake_pcode1.front();

                    bufflen = pcode0.size() * samps_per_sym;
                    tx_raw_buff0.clear();
                    tx_raw_buff1.clear();
                    tx_raw_buff0.resize(bufflen + ntaps, 0L);
                    tx_raw_buff1.resize(bufflen + ntaps, 0L);
                    for (size_t isym = 0; isym < pcode0.size(); isym++) {
                        tx_raw_buff0[isym * samps_per_sym + ntaps / 2 + samps_per_sym / 2] = std::complex<float>(
                                pcode0[isym] * 15000, 0x0000);
                        tx_raw_buff1[isym * samps_per_sym + ntaps / 2 + samps_per_sym / 2] = std::complex<float>(
                                pcode1[isym] * 15000, 0x0000);
                    }

                    //filter the raw tx vector
                    tx_filt_buff0.resize(bufflen, 0);
                    tx_filt_buff1.resize(bufflen, 0);
                    for (int i = 0; i < bufflen; i++) {
                        std::complex<float> temp0(0, 0);
                        std::complex<float> temp1(0, 0);
                        for (int j = 0; j < ntaps; j++) {
                            temp0 += filter_taps[j] * tx_raw_buff0[i + j];
                            temp1 += filter_taps[j] * tx_raw_buff1[i + j];
                        }
                        tx_filt_buff0[i] = std::complex<int16_t>((int16_t) temp0.real(), (int16_t) temp0.imag());
                        tx_filt_buff1[i] = std::complex<int16_t>((int16_t) temp1.real(), (int16_t) temp1.imag());
                    }

                    tx_ontime = (float) tx_filt_buff0.size() / (TX_RATE) + 100e-6;

                    BOOST_LOG_TRIVIAL(debug) << "tx_ontime: " << tx_ontime << std::endl;
                    if (verbose) myfile << "tx_ontime: " << tx_ontime << std::endl;

                    tx_filt_buff0.resize((ipp_usec * TX_RATE / 1e6), 0);
                    tx_filt_buff1.resize((ipp_usec * TX_RATE / 1e6), 0);

                    BOOST_LOG_TRIVIAL(info) << "ipp_usec: " << ipp_usec << std::endl;
                    BOOST_LOG_TRIVIAL(info) << "TX_RATE: " << TX_RATE << std::endl;
                    BOOST_LOG_TRIVIAL(info) << "tx_filt_buffx size: " << ipp_usec * TX_RATE << std::endl;
                    BOOST_LOG_TRIVIAL(info) << "tx_filt_buffx size: " << (size_t)(ipp_usec * TX_RATE / 1e6) << std::endl;

                    //prepare rx information
                    ptime_eff = 3e8 / (2 * MAX_VELOCITY * freq);

                    BOOST_LOG_TRIVIAL(info) << "ptime_eff: " << ptime_eff;
                    if (verbose) myfile << "ptime_eff: " << ptime_eff << std::endl;

                    nave = 2;
                    ptime_eff /= 2;
                    while (ptime_eff > 1.e-6 * ipp_usec && parms.num_pulses / nave > 1) {
                        nave *= 2;
                        ptime_eff /= 2;
                    }

                    if (verbose) myfile << "nave: " << nave << std::endl;
                    BOOST_LOG_TRIVIAL(info) << "nave: " << nave;

                    ptime_eff = nave * 1e-6 * ipp_usec;

                    if (verbose) myfile << "ptime_eff: " << ptime_eff << std::endl;
                    BOOST_LOG_TRIVIAL(info) << "ptime_eff: " << ptime_eff;

                    slowdim = parms.num_pulses / nave;

                    if (verbose) myfile << "slowdim: " << slowdim << std::endl;
                    BOOST_LOG_TRIVIAL(info) << "slowdim: " << slowdim;


                    rawvecs.resize(usrp->get_rx_num_channels());
                    rawvec_ptrs.resize(usrp->get_rx_num_channels());

                    // 170525 GLB trying to figure out why have to restart server for consistent results between soundings
                    //usrp->set_time_now(uhd::time_spec_t(100.0));

                    for (size_t tloop = 1; tloop < 31; tloop++) {

                        BOOST_LOG_TRIVIAL(info) << "Transmit Loop: " << tloop;

                        for (size_t i = 0; i < rawvecs.size(); i++) {
                            rawvecs[i].resize(nsamps_per_pulse * parms.num_pulses * tloop);
                            rawvec_ptrs[i] = (tloop - 1) * (nsamps_per_pulse * parms.num_pulses) + &rawvecs[i].front();
                        }

                        transceive(
                                usrp,
                                tx_stream,
                                rx_stream,
                                parms.num_pulses,
                                1.e-6 * ipp_usec,
                                &tx_filt_buff0,
                                &tx_filt_buff1,
                                tx_ontime,
                                &rawvec_ptrs.front(),
                                nsamps_per_pulse
                        );
                    }
                    slowdim *= 30;


                    BOOST_LOG_TRIVIAL(info) << "Done receiving, waiting for transmit thread...";
                    if (verbose) myfile << "Done receiving, waiting for transmit thread.." << std::endl;

                    transmit_thread.join_all();

                    if (return_status) {
                        BOOST_LOG_TRIVIAL(error) << "This is a bad record..\n";
                        //Do something..?
                    }

                    parms.range_res_km = 1.5e-1 * symboltime_usec;
                    fastdim = nsamps_per_pulse / dmrate;

                    BOOST_LOG_TRIVIAL(info) << "Done Receiving";
                    if (verbose) myfile << "Done rxing 0\n";

                    memset(&actual_parms, 0, sizeof(actual_parms));
                    actual_parms.freq_khz = freq / 1e3;
                    actual_parms.num_pulses = parms.num_pulses;
                    actual_parms.first_range_km = -1;
                    actual_parms.last_range_km = -1;
                    actual_parms.range_res_km = parms.range_res_km;

                    send(msgsock, &actual_parms, sizeof(actual_parms), 0);
                    send(msgsock, &return_status, sizeof(return_status), 0);
                    BOOST_LOG_TRIVIAL(info) << "Finished SEND Command";
                    break;

                case PROCESS:
                    BOOST_LOG_TRIVIAL(info) << "Starting Processing";
                    //dmrate = parms.symboltime_usec * parms.rxrate_khz / (osr*1000);
                    //dmrate = (size_t) (1.e-6*symboltime_usec * RX_RATE / osr);
                    bandwidth = 1 / (2.e-6 * symboltime_usec);

                    BOOST_LOG_TRIVIAL(info) << "Symbol Time: " << symboltime_usec << " usec";
                    BOOST_LOG_TRIVIAL(info) << "fastdim: " << fastdim;
                    BOOST_LOG_TRIVIAL(info) << "dmrate: " << dmrate;
                    BOOST_LOG_TRIVIAL(info) << "bandwidth: " << bandwidth;

                    if (verbose) myfile << "symbol time in usec: " << symboltime_usec << std::endl;
                    if (verbose) myfile << "fastdim: " << fastdim << std::endl;
                    if (verbose) myfile << "dmrate: " << dmrate << std::endl;
                    if (verbose) myfile << "bandwidth: " << bandwidth << std::endl;

                    outvecs[0][0].resize(slowdim);
                    outvecs[0][1].resize(slowdim);
                    outvecs[1][0].resize(slowdim);
                    outvecs[1][1].resize(slowdim);
                    outvec_ptrs[0][0].resize(slowdim);
                    outvec_ptrs[0][1].resize(slowdim);
                    outvec_ptrs[1][0].resize(slowdim);
                    outvec_ptrs[1][1].resize(slowdim);
                    for (int i = 0; i < slowdim; i++) {
                        outvecs[0][0][i].clear();
                        outvecs[0][1][i].clear();
                        outvecs[1][0][i].clear();
                        outvecs[1][1][i].clear();
                        outvecs[0][0][i].resize(nsamps_per_pulse, 0L);
                        outvecs[0][1][i].resize(nsamps_per_pulse, 0L);
                        outvecs[1][0][i].resize(nsamps_per_pulse, 0L);
                        outvecs[1][1][i].resize(nsamps_per_pulse, 0L);
                        outvec_ptrs[0][0][i] = &outvecs[0][0][i].front();
                        outvec_ptrs[0][1][i] = &outvecs[0][1][i].front();
                        outvec_ptrs[1][0][i] = &outvecs[1][0][i].front();
                        outvec_ptrs[1][1][i] = &outvecs[1][1][i].front();
                    }

                    for (int i = 0; i < slowdim; i++) {
                        for (int j = 0; j < nave; j++) {
                            for (int k = 0; k < nsamps_per_pulse; k++) {
                                if (j % 2 == 0) {
                                    outvecs[0][0][i][k] +=
                                            std::complex<int16_t>(1, 0) *
                                            (rawvecs[0][i * nave * nsamps_per_pulse + j * nsamps_per_pulse + k]);
                                    outvecs[1][0][i][k] +=
                                            std::complex<int16_t>(1, 0) *
                                            (rawvecs[1][i * nave * nsamps_per_pulse + j * nsamps_per_pulse + k]);

                                }
                                if (j % 2 == 1) {
                                    outvecs[0][1][i][k] +=
                                            std::complex<int16_t>(1, 0) *
                                            (rawvecs[0][i * nave * nsamps_per_pulse + j * nsamps_per_pulse + k]);
                                    outvecs[1][1][i][k] +=
                                            std::complex<int16_t>(1, 0) *
                                            (rawvecs[1][i * nave * nsamps_per_pulse + j * nsamps_per_pulse + k]);

                                }
                            }

                        }
                    }

                    filtvecs[0][0].resize(slowdim);
                    filtvecs[0][1].resize(slowdim);
                    filtvecs[1][0].resize(slowdim);
                    filtvecs[1][1].resize(slowdim);
                    filtvec_ptrs[0][0].resize(slowdim);
                    filtvec_ptrs[0][1].resize(slowdim);
                    filtvec_ptrs[1][0].resize(slowdim);
                    filtvec_ptrs[1][1].resize(slowdim);
                    for (int i = 0; i < slowdim; i++) {
                        filtvecs[0][0][i].resize(nsamps_per_pulse / dmrate);

                        filtvecs[0][1][i].resize(nsamps_per_pulse / dmrate);
                        filtvecs[1][0][i].resize(nsamps_per_pulse / dmrate);
                        filtvecs[1][1][i].resize(nsamps_per_pulse / dmrate);
                        filtvec_ptrs[0][0][i] = &filtvecs[0][0][i].front();
                        filtvec_ptrs[0][1][i] = &filtvecs[0][1][i].front();
                        filtvec_ptrs[1][0][i] = &filtvecs[1][0][i].front();
                        filtvec_ptrs[1][1][i] = &filtvecs[1][1][i].front();
                    }

                    for (int mode = 0; mode < 2; mode++) {
                        for (int pcode = 0; pcode < 2; pcode++) {
                            rval = lp_filter(
                                    outvec_ptrs[mode][pcode],
                                    filtvec_ptrs[mode][pcode],
                                    slowdim,
                                    nsamps_per_pulse,
                                    RX_RATE,
                                    bandwidth,
                                    dmrate);
                        }
                    }
                    filtvec_dptr[0][0] = &filtvec_ptrs[0][0].front();
                    filtvec_dptr[0][1] = &filtvec_ptrs[0][1].front();
                    filtvec_dptr[1][0] = &filtvec_ptrs[1][0].front();
                    filtvec_dptr[1][1] = &filtvec_ptrs[1][1].front();


                    ffvec_ptrs[0].resize(slowdim);
                    ffvec_ptrs[1].resize(slowdim);
                    ffvecs[0].resize(slowdim);
                    ffvecs[1].resize(slowdim);
                    for (int i = 0; i < slowdim; i++) {
                        ffvecs[0][i].resize(nsamps_per_pulse / dmrate);
                        ffvecs[1][i].resize(nsamps_per_pulse / dmrate);
                        ffvec_ptrs[0][i] = &ffvecs[0][i].front();
                        ffvec_ptrs[1][i] = &ffvecs[1][i].front();
                    }

                    /* Performed matched filtering for both o- and x-mode samples*/
                    for (int mode = 0; mode < 2; mode++) {
                        rval = matched_filter(
                                filtvec_dptr[mode],
                                ffvec_ptrs[mode],
                                fake_pcode_ptrs,
                                pcode0.size(),
                                slowdim,
                                nsamps_per_pulse / dmrate,
                                osr);
                    }

                    first_inx = osr * pcode0.size();
                    filter_delay = osr * pcode0.size() / 2;
                    BOOST_LOG_TRIVIAL(info) << "first inx plus filter delay: " << first_inx + filter_delay;
                    if (verbose) {
                        myfile << "first inx plus filter delay: " << first_inx + filter_delay << std::endl;
                    }

                    fpow[0].resize(nsamps_per_pulse / dmrate, 0);
                    fpow[1].resize(nsamps_per_pulse / dmrate, 0);
                    fvel[0].resize(nsamps_per_pulse / dmrate, 0);
                    fvel[1].resize(nsamps_per_pulse / dmrate, 0);
                    for (int mode = 0; mode < 2; mode++) {
                        rval = doppler_process(
                                &ffvec_ptrs[mode].front(),
                                &fpow[mode].front(),
                                &fvel[mode].front(),
                                slowdim,
                                nsamps_per_pulse / dmrate,
                                1,
                                first_inx + filter_delay);
                    }
                    send(msgsock, &return_status, sizeof(return_status), 0);
                    for (int i = 0; i < nsamps_per_pulse / dmrate; i++) {

                        fpow[0][i] = 10 * log10(fpow[0][i]);
                        fpow[1][i] = 10 * log10(fpow[1][i]);
                        fvel[0][i] = 3e5 * fvel[0][i] / (8. * ptime_eff * slowdim * parms.freq_khz);
                        fvel[1][i] = 3e5 * fvel[1][i] / (8. * ptime_eff * slowdim * parms.freq_khz);
                    }
                    BOOST_LOG_TRIVIAL(info) << "Finished PROCESS Command";
                    break;

                case GET_DATA:
                    BOOST_LOG_TRIVIAL(info) << "Starting Get Data";
                    nranges = fastdim - filter_delay;

                    send(msgsock, &nranges, sizeof(nranges), 0);

                    first_range_km = (pcode0.size() * 1.e-6 * symboltime_usec) * 1.5e5;
                    first_range_km = 0;

                    send(msgsock, &first_range_km, sizeof(first_range_km), 0);

                    last_range_km = (1.e-6 * ipp_usec - pcode0.size() * 1.e-6 * symboltime_usec) * 1.5e5;

                    send(msgsock, &last_range_km, sizeof(last_range_km), 0);
                    send(msgsock,
                         &fpow[0].front() + filter_delay,
                         nranges * sizeof(float), 0);
                    send(msgsock,
                         &fpow[1].front() + filter_delay,
                         nranges * sizeof(float), 0);
                    send(msgsock,
                         &fvel[0].front() + filter_delay,
                         nranges * sizeof(float), 0);
                    send(msgsock,
                         &fvel[1].front() + filter_delay,
                         nranges * sizeof(float), 0);
                    BOOST_LOG_TRIVIAL(info) << "Finished GET_DATA Command";
                    break;

                default:
                    BOOST_LOG_TRIVIAL(fatal) << "Not a valid usrpmsg.  Exiting.\n";
                    return 1;
            }

        }
    }
    //return EXIT_SUCCESS;
}
