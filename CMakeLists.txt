cmake_minimum_required(VERSION 3.18)

project(ff4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED 20)
set(CMAKE_CXX_EXTENSIONS false)


add_library(util
    lib/util
    lib/util/rational.cpp
    lib/util/term.cpp)

add_library(algo
    INTERFACE
    lib/algo/util/
    lib/algo/)

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

add_executable(internal_benchmark benchmarks/internal/main.cpp)
set_target_properties(internal_benchmark PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(internal_benchmark PUBLIC util algo external)

# External benchmarks. To run, install openf4.
# add_executable(external_benchmark benchmarks/external/main.cpp)
# set_target_properties(external_benchmark PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
# add_library(openf4 SHARED IMPORTED)
# if (APPLE)
#    set_target_properties(openf4 PROPERTIES IMPORTED_LOCATION /usr/local/lib/libopenf4.a)
#    target_compile_options(external_benchmark PUBLIC  -I /usr/local/include/)
#    target_link_options(external_benchmark PUBLIC -L/usr/local/lib/ -lopenf4 -Wl)
#elseif (UNIX)
#    set_target_properties(openf4 PROPERTIES IMPORTED_LOCATION /usr/local/lib/libopenf4.so)
#endif()
#target_link_libraries(external_benchmark PUBLIC util algo external openf4)
