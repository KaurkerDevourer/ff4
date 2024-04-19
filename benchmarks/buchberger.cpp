#include "../external/GroebnerBasisFork/GroebnerLib/includes/PolynomialSet.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/Rational.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/Modular.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/F4GB.hpp"
#include "benchmarking.h"
#include "../lib/algo/buchberger.h"
#include "../lib/algo/buchberger_with_criterias.h"
#include "../lib/util/rational.h"
#include "../lib/util/prime_field.h"

using namespace NUtils;

namespace  {
    void FindGroebnerBasisF4(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        for (int i = 0; i < 1000; i++) {
            gb::PolynomialSet<gb::fields::Rational> ideal2 = ideal;
            gb::inplace_calculate_f4_gb(ideal2);
        }
    }

    void FindGroebnerBasisF4Modular(gb::PolynomialSet<gb::fields::Modular<31>>& ideal) {
        for (int i = 0; i < 1000; i++) {
            gb::PolynomialSet<gb::fields::Modular<31>> ideal2 = ideal;
            gb::inplace_calculate_f4_gb(ideal2);
        }
    }

    void FindGroebnerBasisLib(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        for (int i = 0; i < 1000; i++) {
            ideal.MakeGroebnerBasis();
        }
    }

    void FindGroebnerBasisFlex(TPolynomials<Rational>& F) {
        for (int i = 0; i < 1000; i++) {
            TPolynomials<Rational> F2 = F;
            NAlgo::Buchberger::FindGroebnerBasis(F2);
        }
    }

    void FindGroebnerBasisFlex2(TPolynomials<Rational>& F) {
        for (int i = 0; i < 1000; i++) {
            TPolynomials<Rational> F2 = F;
            NAlgo::BuchbergerWithCreterias::FindGroebnerBasis(F2);
        }
    }

    void FindGroebnerBasisFlex2PrimeField(TPolynomials<PrimeField<31>>& F) {
        for (int i = 0; i < 1000; i++) {
            TPolynomials<PrimeField<31>> F2 = F;
            NAlgo::BuchbergerWithCreterias::FindGroebnerBasis(F2);
        }
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


void benchmark_katsura4() {
    {
        TMonomials<PrimeField<31>> amon;
        amon.push_back(Monomial(TTerm({2}), PrimeField<31>(1)));
        amon.push_back(Monomial(TTerm({1}), PrimeField<31>(-1)));
        amon.push_back(Monomial(TTerm({0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(TTerm({0, 0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(TTerm({0, 0, 0, 2}), PrimeField<31>(2)));

        Polynomial<PrimeField<31>> a(std::move(amon));

        TMonomials<PrimeField<31>> bmon;
        bmon.push_back(Monomial(TTerm({1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(TTerm({0, 1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(TTerm({0, 1}), PrimeField<31>(-1)));
        bmon.push_back(Monomial(TTerm({0, 0, 1, 1}), PrimeField<31>(2)));

        Polynomial<PrimeField<31>> b(std::move(bmon));

        TMonomials<PrimeField<31>> cmon;
        cmon.push_back(Monomial(TTerm({1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(TTerm({0, 2}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>> c(std::move(cmon));

        TMonomials<PrimeField<31>> dmon;
        dmon.push_back(Monomial(TTerm({1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(TTerm({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>> d(std::move(dmon));

        TPolynomials<PrimeField<31>> test = {a, b, c, d};
        test_time(FindGroebnerBasisFlex2PrimeField, "buchberger_with_criterion_katsura4 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<31>> i1({
            {{{2}}, 1},
            {{{1}}, -1},
            {{{0, 2}}, 2},
            {{{0, 0, 2}}, 2},
            {{{0, 0, 0, 2}}, 2},
        });

        gb::Polynomial<gb::fields::Modular<31>> i2({
            {{{1, 1}}, 2},
            {{{0, 1, 1}}, 2},
            {{{0, 1}}, -1},
            {{{0, 0, 1, 1}}, 2},
        });

        gb::Polynomial<gb::fields::Modular<31>> i3({
            {{{1, 0, 1}}, 2},
            {{{0, 2}}, 1},
            {{{0, 1, 0, 1}}, 2},
            {{{0, 0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i4({
            {{{1}}, 1},
            {{{0, 1}}, 2},
            {{{0, 0, 1}}, 2},
            {{{0, 0, 0, 1}}, 2},
            {{{0}}, -1},
        });

        gb::PolynomialSet<gb::fields::Modular<31>> ideal({i1, i2, i3, i4});
        test_time(FindGroebnerBasisF4Modular, "GroebnerBasisLibF4_katsura4 ").call(ideal);
    }
}

void benchmark_sym3_3() {
    {
        TMonomials<PrimeField<31>> amon;
        amon.push_back(Monomial(TTerm({1}), PrimeField<31>(1)));
        amon.push_back(Monomial(TTerm({0, 1, 3}), PrimeField<31>(1)));
        amon.push_back(Monomial(TTerm({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>> a(std::move(amon));

        TMonomials<PrimeField<31>> bmon;
        bmon.push_back(Monomial(TTerm({3, 0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(TTerm({0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(TTerm({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>> b(std::move(bmon));

        TMonomials<PrimeField<31>> cmon;
        cmon.push_back(Monomial(TTerm({1, 0, 3}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>> c(std::move(cmon));

        TPolynomials<PrimeField<31>> test = {a, b, c};
        test_time(FindGroebnerBasisFlex2PrimeField, "buchberger_with_criterion_sym3-3 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<31>> i1({
            {{{1}}, 1},
            {{{0, 1, 3}}, 1},
            {{{0}}, -2},
        });

        gb::Polynomial<gb::fields::Modular<31>> i2({
            {{{3, 0, 1}}, 1},
            {{{0, 1}}, 1},
            {{{0}}, -2},
        });

        gb::Polynomial<gb::fields::Modular<31>> i3({
            {{{1, 0, 3}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0}}, -2},
        });

        gb::PolynomialSet<gb::fields::Modular<31>> ideal({i1, i2, i3});
        test_time(FindGroebnerBasisF4Modular, "GroebnerBasisLibF4_sym3-3 ").call(ideal);
    }
    
}