#include "../external/GroebnerBasisFork/GroebnerLib/includes/PolynomialSet.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/Rational.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/F4GB.hpp"
#include "benchmarking.h"
#include "../lib/algo/buchberger.h"
#include "../lib/algo/buchberger_with_criterias.h"
#include "../lib/util/rational.h"

using namespace NUtils;

namespace  {
    void FakeTestForRazogrev(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        for (int i = 0; i < 400; i++) {
            ideal.MakeGroebnerBasis();
        }
    }

    void FindGroebnerBasisF4(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        gb::inplace_calculate_f4_gb(ideal);
    }

    void FindGroebnerBasisLib(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        ideal.MakeGroebnerBasis();
    }

    void FindGroebnerBasisFlex(TPolynomials<Rational>& F) {
        NAlgo::Buchberger::FindGroebnerBasis(F);
    }

    void FindGroebnerBasisFlex2(TPolynomials<Rational>& F) {
        NAlgo::BuchbergerWithCreterias::FindGroebnerBasis(F);
    }
}

void benchmark_cyclic4() {
    {
        gb::Polynomial<gb::fields::Rational> i1({
            {{{1, 1, 1, 1}}, 1},
            {{{0}}, -1},
        });
        gb::Polynomial<gb::fields::Rational> i2({
            {{{1, 1, 1}}, 1},
            {{{1, 1, 0, 1}}, 1},
            {{{1, 0, 1, 1}}, 1},
            {{{0, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Rational> i3({
            {{{1, 1}}, 1},
            {{{1, 0, 0, 1}}, 1},
            {{{0, 1, 1, 0}}, 1},
            {{{0, 0, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Rational> i4({
            {{{1}}, 1},
            {{{0, 1}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0, 0, 0, 1}}, 1},
        });
        gb::PolynomialSet<gb::fields::Rational> ideal({i1, i2, i3, i4});
        test_time(FakeTestForRazogrev, "FakeTestForRazogrev ").call(ideal);
    }
    {
        gb::Polynomial<gb::fields::Rational> i1({
            {{{1, 1, 1, 1}}, 1},
            {{{0}}, -1},
        });
        gb::Polynomial<gb::fields::Rational> i2({
            {{{1, 1, 1}}, 1},
            {{{1, 1, 0, 1}}, 1},
            {{{1, 0, 1, 1}}, 1},
            {{{0, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Rational> i3({
            {{{1, 1}}, 1},
            {{{1, 0, 0, 1}}, 1},
            {{{0, 1, 1, 0}}, 1},
            {{{0, 0, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Rational> i4({
            {{{1}}, 1},
            {{{0, 1}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0, 0, 0, 1}}, 1},
        });
        gb::PolynomialSet<gb::fields::Rational> ideal({i1, i2, i3, i4});
        test_time(FindGroebnerBasisLib, "GroebnerBasisLibBuchberger_cyclic4 ").call(ideal);
    }

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
        test_time(FindGroebnerBasisFlex, "buchberger_cyclic4 ").call(test);
    }
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
        test_time(FindGroebnerBasisFlex2, "buchberger_with_criterion_cyclic4 ").call(test);
    }

    {
        gb::Polynomial<gb::fields::Rational> i1({
            {{{1, 1, 1, 1}}, 1},
            {{{0}}, -1},
        });
        gb::Polynomial<gb::fields::Rational> i2({
            {{{1, 1, 1}}, 1},
            {{{1, 1, 0, 1}}, 1},
            {{{1, 0, 1, 1}}, 1},
            {{{0, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Rational> i3({
            {{{1, 1}}, 1},
            {{{1, 0, 0, 1}}, 1},
            {{{0, 1, 1, 0}}, 1},
            {{{0, 0, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Rational> i4({
            {{{1}}, 1},
            {{{0, 1}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0, 0, 0, 1}}, 1},
        });
        gb::PolynomialSet<gb::fields::Rational> ideal({i1, i2, i3, i4});
        test_time(FindGroebnerBasisF4, "GroebnerBasisLibF4_cyclic4 ").call(ideal);
    }
}
