#include <iostream>

#include "../lib/algo/buchberger.h"
#include "../lib/util/rational.h"

void test_buchberger() {
    using namespace NUtils;
    TMonomials<Rational> amon;
    amon.push_back(Monomial(TTerm({3}), Rational(1)));
    amon.push_back(Monomial(TTerm({1, 1}), Rational(-2)));

    Polynomial<Rational> a(std::move(amon));

    TMonomials<Rational> bmon;
    bmon.push_back(Monomial(TTerm({2, 1}), Rational(1)));
    bmon.push_back(Monomial(TTerm({1}), Rational(1)));
    bmon.push_back(Monomial(TTerm({0, 2}), Rational(-2)));

    Polynomial<Rational> b(std::move(bmon));

    TPolynomials<Rational> test = {a, b};
    std::cout << test << std::endl;
    NAlgo::Buchberger::FindGroebnerBasis(test);
    std::cout << "Size of Groebner basis by Buchberger: " << test.size() << std::endl;
    std::cout << test << std::endl;
}
