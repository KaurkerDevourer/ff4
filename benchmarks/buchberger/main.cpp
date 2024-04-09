#include "../../GroebnerBasis/GroebnerLib/includes/PolynomialSet.hpp"
#include "../../GroebnerBasis/GroebnerLib/includes/Rational.hpp"
#include "../../GroebnerBasis/GroebnerLib/includes/F4GB.hpp"
#include "../benchmarking.h"
#include "../../lib/algo/buchberger.h"

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

    void FindGroebnerBasisFlex(NUtils::TPolynomials& F) {
        for (int i = 0; i < 1000; i++) {
            NAlgo::Buchberger::FindGroebnerBasis(F);
        }
    }
}

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
    test_time(FindGroebnerBasisFlex, "buchberger_small ").call(test);
    //FakeBenchmark(NAlgo::Buchberger::FindGroebnerBasis, 10, "buchberger_small", test);

    gb::Polynomial<gb::fields::Rational> i1({  // HW 07, ex 01
        {{{3}}, 1},
        {{{1, 1}}, -2},
    });
    gb::Polynomial<gb::fields::Rational> i2({
        {{{2, 1}}, 1},
        {{{1}}, 1},
        {{{0, 2}}, -2},
    });
    //gb::PolynomialSet<gb::fields::Rational> ideal({i1, i2});
    //FakeBenchmark(GroebnerBasisLibF4::FindGroebnerBasisF4, 10, "GroebnerBasisLibF4_small", ideal);


    gb::PolynomialSet<gb::fields::Rational> ideal2({i1, i2});
    test_time(GroebnerBasisLibBuchberger::FindGroebnerBasisLib, "GroebnerBasisLibBuchberger_small ").call(ideal2);
    //FakeBenchmark(GroebnerBasisLibBuchberger::FindGroebnerBasis, 10, "GroebnerBasisLibBuchberger_small", ideal2);
}
