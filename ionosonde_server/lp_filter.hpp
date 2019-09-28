//
// Created by bgklug on 9/27/19.
//

#ifndef IONOSONDE_SERVER_LP_FILTER_HPP
#define IONOSONDE_SERVER_LP_FILTER_HPP

#include <vector>
#include <complex>

int lp_filter(
        std::vector<std::complex <float> *> indata,
        std::vector<std::complex < float> *> outdata,
        int slowdim,
        int fastdim,
        float samprate,
        float bw,
        int decimrate
);

#endif //IONOSONDE_SERVER_LP_FILTER_HPP
