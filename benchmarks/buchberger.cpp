#include "../external/GroebnerBasisFork/GroebnerLib/includes/PolynomialSet.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/Rational.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/F4GB.hpp"
#include "benchmarking.h"
#include "../lib/algo/buchberger.h"
#include "../lib/util/rational.h"

using namespace NUtils;

namespace  {
    void FindGroebnerBasisF4(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        for (int i = 0; i < 1000; i++) {
            gb::inplace_calculate_f4_gb(ideal);
        }
    }

    void FindGroebnerBasisLib(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        for (int i = 0; i < 1000; i++) {
            ideal.MakeGroebnerBasis();
        }
    }

    void FindGroebnerBasisFlex(TPolynomials<Rational>& F) {
        for (int i = 0; i < 1000; i++) {
            NAlgo::Buchberger::FindGroebnerBasis(F);
        }
    }
}

void benchmark_buchberger() {
    TMonomials<Rational> amon;
    amon.push_back(Monomial({3}, Rational(1)));
    amon.push_back(Monomial({1, 1}, Rational(-2)));

    Polynomial a(std::move(amon));

    TMonomials<Rational> bmon;
    bmon.push_back(Monomial({2, 1}, Rational(1)));
    bmon.push_back(Monomial({1}, Rational(1)));
    bmon.push_back(Monomial({0, 2}, Rational(-2)));

    Polynomial b(std::move(bmon));

    TPolynomials<Rational> test = {a, b};
    test_time(FindGroebnerBasisFlex, "buchberger_small ").call(test);

    gb::Polynomial<gb::fields::Rational> i1({  // HW 07, ex 01
        {{{3}}, 1},
        {{{1, 1}}, -2},
    });
    gb::Polynomial<gb::fields::Rational> i2({
        {{{2, 1}}, 1},
        {{{1}}, 1},
        {{{0, 2}}, -2},
    });

    gb::PolynomialSet<gb::fields::Rational> ideal2({i1, i2});
    test_time(FindGroebnerBasisLib, "GroebnerBasisLibBuchberger_small ").call(ideal2);

    gb::PolynomialSet<gb::fields::Rational> ideal({i1, i2});
    test_time(FindGroebnerBasisF4, "GroebnerBasisLibF4_small ").call(ideal);
}
