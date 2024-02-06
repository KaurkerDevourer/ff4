#include <iostream>

#include "buchberger.h"

int main() {
    TPolynomial a = {
        {{3}, 1}, 
        {{1, 1} - 2}
    };
    TPolynomial b = {
        {{2, 1}, 1}, 
        {{1}, 1},
        {{0, 2} -2},
    };
    TPolynomials test = {a, b};
    DoProcess(test);
    std::cout << test.size() << std::endl;
}