#include "../external/GroebnerBasisFork/GroebnerLib/includes/PolynomialSet.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/Rational.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/Modular.hpp"
#include "../external/GroebnerBasisFork/GroebnerLib/includes/F4GB.hpp"
#include "benchmarking.h"
#include "../lib/algo/buchberger.h"
#include "../lib/algo/buchberger_with_criteria.h"
#include "../lib/algo/f4.h"
#include "../lib/util/rational.h"
#include "../lib/util/prime_field.h"
#include <libopenf4.h>

using namespace FF4;
using namespace FF4::NUtils;

namespace  {
    #define TimesToRun 1000
    void FindGroebnerBasisF4Lib(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                gb::PolynomialSet<gb::fields::Rational> ideal2 = ideal;
                gb::inplace_calculate_f4_gb(ideal2);
            }
        #else
            gb::inplace_calculate_f4_gb(ideal);
            std::cout << ideal << std::endl;
        #endif
    }

    void FindGroebnerBasisF4LibModular(gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp>& ideal) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp> ideal2 = ideal;
                gb::inplace_calculate_f4_gb(ideal2);
            }
        #else
            gb::inplace_calculate_f4_gb(ideal);
            std::cout << ideal << std::endl;
        #endif
    }

    void FindGroebnerBasisLib(gb::PolynomialSet<gb::fields::Rational>& ideal) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                ideal.MakeGroebnerBasis();
            }
        #else
            ideal.MakeGroebnerBasis();
        #endif
    }

    void FindGroebnerBasis(TPolynomials<Rational, LexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<Rational, LexComp> F2 = F;
                NAlgo::Buchberger::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::Buchberger::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisWithCriterias(TPolynomials<Rational, LexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<Rational, LexComp> F2 = F;
                NAlgo::BuchbergerWithCreteria::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::BuchbergerWithCreteria::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisWithCriteriasPrimeField(TPolynomials<PrimeField<31>, GrevLexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<PrimeField<31>, GrevLexComp> F2 = F;
                NAlgo::BuchbergerWithCreteria::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::BuchbergerWithCreteria::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisF4(TPolynomials<Rational, LexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<Rational, LexComp> F2 = F;
                NAlgo::F4::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::F4::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisF4PrimeField(TPolynomials<PrimeField<31>, GrevLexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<PrimeField<31>, GrevLexComp> F2 = F;
                NAlgo::F4::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::F4::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisOpenF4PrimeField(std::vector<std::string> variableName, std::vector<std::string> polynomialList) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                std::vector<std::string> basis = groebnerBasisF4(31, variableName.size(), variableName, polynomialList, 1, 0);
            }
        #else
            std::vector<std::string> basis = groebnerBasisF4(31, variableName.size(), variableName, polynomialList, 1, 4);
            std::cout << basis.size() << std::endl;
            for (const auto& str : basis) {
                std::cout << str << std::endl;
            }
        #endif
    }
}

void benchmark_cyclic4_rational() {
    {
        std::vector<Monomial<Rational>> amon;
        amon.push_back(Monomial(TTerm({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(TTerm({0}), Rational(-1)));

        Polynomial<Rational, LexComp> a(std::move(amon));

        std::vector<Monomial<Rational>> bmon;
        bmon.push_back(Monomial(TTerm({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> b(std::move(bmon));

        std::vector<Monomial<Rational>> cmon;
        cmon.push_back(Monomial(TTerm({1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> c(std::move(cmon));

        std::vector<Monomial<Rational>> dmon;
        dmon.push_back(Monomial(TTerm({1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational, LexComp> d(std::move(dmon));

        TPolynomials<Rational, LexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasis, "buchberger_cyclic4 ").call(test);
    }
    {
        std::vector<Monomial<Rational>> amon;
        amon.push_back(Monomial(TTerm({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(TTerm({0}), Rational(-1)));

        Polynomial<Rational, LexComp> a(std::move(amon));

        std::vector<Monomial<Rational>> bmon;
        bmon.push_back(Monomial(TTerm({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(TTerm({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> b(std::move(bmon));

        std::vector<Monomial<Rational>> cmon;
        cmon.push_back(Monomial(TTerm({1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> c(std::move(cmon));

        std::vector<Monomial<Rational>> dmon;
        dmon.push_back(Monomial(TTerm({1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational, LexComp> d(std::move(dmon));

        TPolynomials<Rational, LexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisWithCriterias, "buchberger_with_criterion_cyclic4_rational ").call(test);
    }
}


void benchmark_cyclic4_prime_field() {
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(TTerm({1, 1, 1, 1}), PrimeField<31>(1)));
        amon.push_back(Monomial(TTerm({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(TTerm({1, 1, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(TTerm({1, 1, 0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(TTerm({1, 0, 1, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(TTerm({0, 1, 1, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(TTerm({1, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0, 1, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<31>>> dmon;
        dmon.push_back(Monomial(TTerm({1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisF4PrimeField, "f4_cyclic4 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<31>> i1({
            {{{1, 1, 1, 1}}, 1},
            {{{0}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i2({
            {{{1, 1, 1}}, 1},
            {{{1, 1, 0, 1}}, 1},
            {{{1, 0, 1, 1}}, 1},
            {{{0, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i3({
            {{{1, 1}}, 1},
            {{{0, 1, 1}}, 1},
            {{{1, 0, 0, 1}}, 1},
            {{{0, 0, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i4({
            {{{1}}, 1},
            {{{0, 1}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0, 0, 0, 1}}, 1},
        });

        gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp> ideal({i1, i2, i3, i4});
        test_time(FindGroebnerBasisF4LibModular, "GroebnerBasisLibF4_cyclic4 ").call(ideal);
    }
    {
        std::vector<std::string> polynomialList;
        polynomialList.emplace_back("x0+x1+x2+x3");
        polynomialList.emplace_back("x0*x1+x1*x2+x2*x3+x0*x3");
        polynomialList.emplace_back("x0*x1*x2+x1*x2*x3+x2*x3*x0+x0*x1*x3");
        polynomialList.emplace_back("x0*x1*x2*x3-1");
        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3"};
        test_time(FindGroebnerBasisOpenF4PrimeField, "openf4_cyclic4 ").call(variableName, polynomialList);
    }
}


void benchmark_katsura4() {
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(TTerm({2}), PrimeField<31>(1)));
        amon.push_back(Monomial(TTerm({0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(TTerm({0, 0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(TTerm({0, 0, 0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(TTerm({1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(TTerm({1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(TTerm({0, 1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(TTerm({0, 0, 1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(TTerm({0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(TTerm({0, 2}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(TTerm({0, 1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<31>>> dmon;
        dmon.push_back(Monomial(TTerm({1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(TTerm({0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(TTerm({0, 0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(TTerm({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisF4PrimeField, "f4_katsura4 ").call(test);
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

        gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp> ideal({i1, i2, i3, i4});
        test_time(FindGroebnerBasisF4LibModular, "GroebnerBasisLibF4_katsura4 ").call(ideal);
    }
    {
        std::vector<std::string> polynomialList;
        polynomialList.emplace_back("a^2-a+2*b^2+2*c^2+2*d^2");
        polynomialList.emplace_back("2*a*b+2*b*c-b+2*c*d");
        polynomialList.emplace_back("2*a*c+b^2+2*b*d-c");
        polynomialList.emplace_back("a+2*b+2*c+2*d-1");
        std::vector<std::string> variableName = {"a", "b", "c", "d"};
        test_time(FindGroebnerBasisOpenF4PrimeField, "openf4_katsura4 ").call(variableName, polynomialList);
    }
}

void benchmark_sym3_3() {
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(TTerm({0, 1, 3}), PrimeField<31>(1)));
        amon.push_back(Monomial(TTerm({1}), PrimeField<31>(1)));
        amon.push_back(Monomial(TTerm({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(TTerm({3, 0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(TTerm({0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(TTerm({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(TTerm({1, 3}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c};
        test_time(FindGroebnerBasisF4PrimeField, "f4_sym3-3 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<31>> i1({
            {{{0, 1, 3}}, 1},
            {{{1}}, 1},
            {{{0}}, -2},
        });

        gb::Polynomial<gb::fields::Modular<31>> i2({
            {{{3, 0, 1}}, 1},
            {{{0, 1}}, 1},
            {{{0}}, -2},
        });

        gb::Polynomial<gb::fields::Modular<31>> i3({
            {{{1, 3}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0}}, -2},
        });

        gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp> ideal({i1, i2, i3});
        test_time(FindGroebnerBasisF4LibModular, "GroebnerBasisLibF4_sym3-3 ").call(ideal);
    }
    {
        std::vector<std::string> polynomialList;
        polynomialList.emplace_back("a+b*c^3-2");
        polynomialList.emplace_back("a^3*c+b-2");
        polynomialList.emplace_back("a*b^3+c-2");
        std::vector<std::string> variableName = {"a", "b", "c"};
        test_time(FindGroebnerBasisOpenF4PrimeField, "openf4_sym3-3 ").call(variableName, polynomialList);
    }
    
}
