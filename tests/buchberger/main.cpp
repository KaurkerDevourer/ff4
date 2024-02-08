#include <iostream>

#include "../../lib/algo/buchberger.h"

int main() {
    NUtils::TMonomials amon;
    amon.push_back(NUtils::Monomial({3}, NUtils::Rational(1)));
    amon.push_back(NUtils::Monomial({1, 1}, NUtils::Rational(-2)));

    NUtils::Polynomial a(std::move(amon));

    NUtils::TMonomials bmon;
    bmon.push_back(NUtils::Monomial({2, 1}, NUtils::Rational(1)));
    bmon.push_back(NUtils::Monomial({1}, NUtils::Rational(1)));
    bmon.push_back(NUtils::Monomial({0, 2}, NUtils::Rational(-2)));

    NUtils::Polynomial b(std::move(bmon));

    NUtils::TPolynomials test = {a, b};
    NAlgo::Buchberger::DoProcess(test);
    std::cout << "Size of Graebner basis by Buchberger: " << test.size() << std::endl;
}