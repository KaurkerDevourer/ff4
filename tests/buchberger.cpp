#include <iostream>

#include "../lib/algo/buchberger.h"
#include "../lib/algo/buchberger_with_criterias.h"
#include "../lib/algo/buchberger_with_criterias_version2.h"
#include "../lib/util/rational.h"
#include "../lib/algo/util/groebner_basis_util.h"
#include "../lib/util/comp.hpp"

void test_buchberger() {
    using namespace NUtils;
    // cyclic-4-
    {
        TMonomials<Rational> amon;
        amon.push_back(Monomial(TTerm({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(TTerm({0}), Rational(-1)));

        Polynomial<Rational, LexComp> a(std::move(amon));

        TMonomials<Rational> bmon;
        bmon.push_back(Monomial(TTerm({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> b(std::move(bmon));

        TMonomials<Rational> cmon;
        cmon.push_back(Monomial(TTerm({1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> c(std::move(cmon));

        TMonomials<Rational> dmon;
        dmon.push_back(Monomial(TTerm({1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational, LexComp> d(std::move(dmon));

        TPolynomials<Rational, LexComp> test = {a, b, c, d};
        std::cout << "Buchberger: " << test << std::endl;
        NAlgo::Buchberger::FindGroebnerBasis(test);
        assert(NAlgo::NUtil::CheckBasisIsGroebner(test));
        std::cout << "Size of Groebner basis by Buchberger: " << test.size() << std::endl;
        std::cout << test << std::endl;
    }

    {
        TMonomials<Rational> amon;
        amon.push_back(Monomial(TTerm({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(TTerm({0}), Rational(-1)));

        Polynomial<Rational, LexComp> a(std::move(amon));

        TMonomials<Rational> bmon;
        bmon.push_back(Monomial(TTerm({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> b(std::move(bmon));

        TMonomials<Rational> cmon;
        cmon.push_back(Monomial(TTerm({1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> c(std::move(cmon));

        TMonomials<Rational> dmon;
        dmon.push_back(Monomial(TTerm({1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational, LexComp> d(std::move(dmon));

        TPolynomials<Rational, LexComp> test = {a, b, c, d};
        std::cout << "Buchberger with criterias: " << test << std::endl;
        NAlgo::BuchbergerWithCreterias::FindGroebnerBasis(test);
        assert(NAlgo::NUtil::CheckBasisIsGroebner(test));
        std::cout << "Size of Groebner basis by Buchberger with criterias: " << test.size() << std::endl;
        std::cout << test << std::endl;
    }

    {
        TMonomials<Rational> amon;
        amon.push_back(Monomial(TTerm({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(TTerm({0}), Rational(-1)));

        Polynomial<Rational, LexComp> a(std::move(amon));

        TMonomials<Rational> bmon;
        bmon.push_back(Monomial(TTerm({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> b(std::move(bmon));

        TMonomials<Rational> cmon;
        cmon.push_back(Monomial(TTerm({1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> c(std::move(cmon));

        TMonomials<Rational> dmon;
        dmon.push_back(Monomial(TTerm({1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational, LexComp> d(std::move(dmon));

        TPolynomials<Rational, LexComp> test = {a, b, c, d};
        std::cout << "Buchberger with criterias version 2: " << test << std::endl;
        NAlgo::BuchbergerWithCreteriasVersion2::FindGroebnerBasis(test);
        assert(NAlgo::NUtil::CheckBasisIsGroebner(test));
        std::cout << "Size of Groebner basis by Buchberger with criterias version 2: " << test.size() << std::endl;
        std::cout << test << std::endl;
    }
}
