//
// Created by bgklug on 9/27/19.
//

#ifndef IONOSONDE_SERVER_MATCHED_FILTER_HPP
#define IONOSONDE_SERVER_MATCHED_FILTER_HPP

#include <vector>
#include <complex>

int matched_filter(
        std::complex<float> ***indata,
        std::vector<std::complex <float> *> outdata,
        float **pcode,
        int pcode_length,
        int slowdim,
        int fastdim,
        int osr
);

#endif //IONOSONDE_SERVER_MATCHED_FILTER_HPP
