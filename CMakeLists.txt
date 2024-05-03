cmake_minimum_required(VERSION 3.18)

project(fast_f4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED 20)
set(CMAKE_CXX_EXTENSIONS false)


add_library(util
    lib/util/critical_pair.hpp
    lib/util/comp.hpp
    lib/util/monomial.hpp
    lib/util/polynomial.hpp
    lib/util/prime_field.hpp
    lib/util/rational.cpp
    lib/util/rational.h
    lib/util/term.cpp
    lib/util/term.h)

add_library(algo
    lib/algo/util/groebner_basis_util.cpp
    lib/algo/util/groebner_basis_util.h
    lib/algo/util/matrix_reduction.cpp
    lib/algo/util/matrix_reduction.h
    lib/algo/f4.cpp
    lib/algo/f4.h
    lib/algo/buchberger.cpp
    lib/algo/buchberger.h
    lib/algo/buchberger_with_criterias.cpp
    lib/algo/buchberger_with_criterias.h
    lib/algo/buchberger_with_criterias_version2.cpp
    lib/algo/buchberger_with_criterias_version2.h)

add_library(external
    external/GroebnerBasisFork/GroebnerLib/includes/PolynomialSet.hpp
    external/GroebnerBasisFork/GroebnerLib/includes/Rational.hpp
    external/GroebnerBasisFork/GroebnerLib/includes/Modular.hpp
    external/GroebnerBasisFork/GroebnerLib/includes/F4GB.hpp
    external/GroebnerBasisFork/GroebnerLib/srcs/Lib.cpp
    external/GroebnerBasisFork/GroebnerLib/srcs/Rational.cpp
    external/GroebnerBasisFork/GroebnerLib/srcs/Term.cpp)

add_executable(ut-tests tests/main.cpp)
set_target_properties(ut-tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(ut-tests PUBLIC util algo)

add_executable(benchmark benchmarks/main.cpp)
set_target_properties(benchmark PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(benchmark PUBLIC util algo external)
