#include "benchmarking.h"
#include "groebner_basis.cpp"
#include <iostream>

int main() {
    #ifdef NDEBUG
        //nothing
    #else
        std::cout << "DEBUG MODE" << std::endl;
    #endif
    benchmark_cyclic4();
    benchmark_katsura4();
    benchmark_sym3_3();
}
