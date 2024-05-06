#include "../lib/util/polynomial.h"
#include "testing.h"
#include <iostream>
#include <cassert>

void test_polynomial() {
    using namespace FF4::NUtils;
    std::vector<Monomial<Rational>> amon;
    amon.push_back(Monomial(TTerm({3}), Rational(1)));
    amon.push_back(Monomial(TTerm({1, 1}), Rational(-2)));

    Polynomial<Rational, LexComp> a(std::move(amon));
    Polynomial<Rational, LexComp> b = a * Rational(2);
    ASSERT_EQUAL(b.GetMonomials()[0], Monomial(TTerm({3}), Rational(2)));
    ASSERT_EQUAL(b.GetMonomials()[1], Monomial(TTerm({1, 1}), Rational(-4)));
    b /= Rational(2);
    ASSERT_EQUAL(a, b);

    Monomial c(TTerm({0, 2}), Rational(3));
    Polynomial<Rational, LexComp> d = a * c;
    ASSERT_EQUAL(d.GetMonomials()[0], Monomial(TTerm({3, 2}), Rational(3)));
    ASSERT_EQUAL(d.GetMonomials()[1], Monomial(TTerm({1, 3}), Rational(-6)));


    std::vector<Monomial<Rational>> emon;
    emon.push_back(Monomial(TTerm({1}), Rational(1)));
    emon.push_back(Monomial(TTerm({0, 1}), Rational(1)));
    Polynomial<Rational, LexComp> e(std::move(emon));

    std::vector<Monomial<Rational>> fmon;
    fmon.push_back(Monomial(TTerm({1}), Rational(1)));
    fmon.push_back(Monomial(TTerm({0, 1}), Rational(-1)));
    Polynomial<Rational, LexComp> f(std::move(fmon));

    std::vector<Monomial<Rational>> kmon;
    kmon.push_back(Monomial(TTerm({2}), Rational(1)));
    kmon.push_back(Monomial(TTerm({0, 2}), Rational(-1)));
    Polynomial<Rational, LexComp> k(std::move(kmon));
    ASSERT_EQUAL(e * f, k);

    std::cout << "Successfully tested Polynomial" << std::endl;
}
