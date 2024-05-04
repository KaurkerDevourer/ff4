#include <iostream>

#include "../lib/algo/f4.h"
#include "../lib/util/rational.h"
#include "../lib/algo/util/groebner_basis_util.h"

void test_f4() {
    using namespace NUtils;
    // cyclic-4
    {
        TMonomials<Rational> amon;
        amon.push_back(Monomial(TTerm({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(TTerm({0}), Rational(-1)));

        Polynomial<Rational> a(std::move(amon));

        TMonomials<Rational> bmon;
        bmon.push_back(Monomial(TTerm({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational> b(std::move(bmon));

        TMonomials<Rational> cmon;
        cmon.push_back(Monomial(TTerm({1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational> c(std::move(cmon));

        TMonomials<Rational> dmon;
        dmon.push_back(Monomial(TTerm({1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational> d(std::move(dmon));

        TPolynomials<Rational> test = {a, b, c, d};
        std::cout << "F4: " << test << std::endl;
        NAlgo::F4::FindGroebnerBasis(test);
        assert(NAlgo::NUtil::CheckBasisIsGroebner(test));
        std::cout << "Size of Groebner basis by F4: " << test.size() << std::endl;
        std::cout << test << std::endl;
    }
}
