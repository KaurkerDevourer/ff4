#include "../benchmarking.h"
#include "groebner_basis.cpp"
#include <iostream>

int main() {
    #ifdef NDEBUG
        //nothing
    #else
        std::cout << "DEBUG MODE" << std::endl;
    #endif
    // benchmark_buchberger_cyclic4_rational();
    // benchmark_buchberger_cyclic4();
    // benchmark_buchberger_katsura4();
    // benchmark_buchberger_sym3_3();
    // benchmark_buchberger_katsura5();

    // benchmark_cyclic4();
    // benchmark_katsura4();
    // benchmark_sym3_3();

    benchmark_cyclic5();
    benchmark_cyclic6();
    benchmark_cyclic7();

    benchmark_katsura5();
    benchmark_katsura9();
    benchmark_katsura10();
    benchmark_katsura11();
    benchmark_katsura12();
}
