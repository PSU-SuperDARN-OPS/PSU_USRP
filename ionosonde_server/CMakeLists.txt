cmake_minimum_required(VERSION 3.5.1...3.14)
project(ionosonde_server VERSION 1.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(FFTW)

find_package(Boost "1.46" REQUIRED
        program_options
        thread)

find_package(UHD "3.9.7")


set(CMAKE_CXX_STANDARD 14)


add_executable(ionosonde_server
        main.cpp
        capture_spectrum.cpp
        lp_filter.cpp
        doppler_process.cpp
        matched_filter.cpp
        transceive.cpp
        utils.cpp utils.hpp
        sounder.cpp sounder.hpp
        global_variables.h
        )

target_link_libraries(ionosonde_server fftw3)