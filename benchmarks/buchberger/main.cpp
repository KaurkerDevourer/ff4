#include "../../GroebnerBasis/GroebnerLib/includes/PolynomialSet.hpp"
#include "../../GroebnerBasis/GroebnerLib/includes/Rational.hpp"
#include "../../GroebnerBasis/GroebnerLib/includes/F4GB.hpp"
#include "../benchmarking.h"
#include "../../lib/algo/buchberger.h"

namespace GroebnerBasisLibF4 {
    void FindGroebnerBasis(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        gb::inplace_calculate_f4_gb(ideal);
    }
}

namespace GroebnerBasisLibBuchberger {
    void FindGroebnerBasis(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        ideal.MakeGroebnerBasis();
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
    FakeBenchmark(NAlgo::Buchberger::FindGroebnerBasis, 10, "buchberger_small", test);

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
    //FakeBenchmark(GroebnerBasisLibF4::FindGroebnerBasis, 10, "GroebnerBasisLibF4_small", ideal);


    gb::PolynomialSet<gb::fields::Rational> ideal2({i1, i2});
    FakeBenchmark(GroebnerBasisLibBuchberger::FindGroebnerBasis, 10, "GroebnerBasisLibBuchberger_small", ideal2);
}
