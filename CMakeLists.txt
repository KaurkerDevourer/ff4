cmake_minimum_required(VERSION 3.18)

project(fast_f4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)


add_library(util 
    lib/util/monomial.cpp
    lib/util/monomial.h
    lib/util/polynomial.cpp
    lib/util/polynomial.h
    lib/util/rational.cpp
    lib/util/rational.h
    lib/util/term.cpp
    lib/util/term.h)

add_library(algo 
    lib/algo/buchberger.cpp
    lib/algo/buchberger.h)

# add_library(external
#     external/GroebnerBasis/GroebnerLib/includes/PolynomialSet.hpp
#     external/GroebnerBasis/GroebnerLib/includes/Rational.hpp
#     external/GroebnerBasis/GroebnerLib/includes/F4GB.hpp)

add_executable(ut-tests tests/main.cpp)
set_target_properties(ut-tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(ut-tests PUBLIC util algo)

add_executable(benchmark benchmarks/main.cpp)
set_target_properties(benchmark PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(benchmark PUBLIC util algo)
