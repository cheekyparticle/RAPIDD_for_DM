cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_FLAGS_INIT "-Werror -Wall -Wextra -Wpedantic")

list(APPEND CMAKE_REQUIRED_LIBRARIES m)

# set the project name
project(RAPIDD VERSION 1.0)
configure_file(RAPIDDConfig.h.in RAPIDDConfig.h)

add_library(RAPIDD SHARED pcubature.c hcubature.c source/halo.c source/effFormFact.c source/difCrosSec.c source/WriteData.c source/coeffs_eft.c source/binning_general.c source/difRateGen.c source/phys_consts.c source/efficiency_curve.c source/specific_exps.c)

find_package(GSL REQUIRED)
target_link_libraries(RAPIDD PUBLIC GSL::gsl GSL::gslcblas)
# add the executable
add_executable(rpd main.c)
target_link_libraries(rpd PUBLIC RAPIDD GSL::gsl GSL::gslcblas m)