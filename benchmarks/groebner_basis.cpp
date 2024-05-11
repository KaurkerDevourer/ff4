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

    void FindGroebnerBasisF4LibModular(gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp>& ideal) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp> ideal2 = ideal;
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

    void FindGroebnerBasisWithCriteriasPrimeField(TPolynomials<PrimeField<1000000007>, GrevLexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<PrimeField<1000000007>, GrevLexComp> F2 = F;
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

    void FindGroebnerBasisF4PrimeField(TPolynomials<PrimeField<1000000007>, GrevLexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<PrimeField<1000000007>, GrevLexComp> F2 = F;
                NAlgo::F4::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::F4::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }


    void FindGroebnerBasisF4PrimeFieldBig(TPolynomials<PrimeField<1000000007>, GrevLexComp>& F, int timesToRun) {
        #ifdef NDEBUG
            for (int i = 0; i < timesToRun; i++) {
                TPolynomials<PrimeField<1000000007>, GrevLexComp> F2 = F;
                NAlgo::F4::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::F4::FindGroebnerBasis(F);
            std::cout << F.size() << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebnerBig(F));
        #endif
    }

    void FindGroebnerBasisOpenF4PrimeField(std::vector<std::string> variableName, std::vector<std::string> polynomialList) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                std::vector<std::string> basis = groebnerBasisF4(1000000007, variableName.size(), variableName, polynomialList, 1, 0);
            }
        #else
            std::vector<std::string> basis = groebnerBasisF4(1000000007, variableName.size(), variableName, polynomialList, 1, 4);
            std::cout << basis.size() << std::endl;
            for (const auto& str : basis) {
                std::cout << str << std::endl;
            }
        #endif
    }

    void FindGroebnerBasisOpenF4PrimeFieldBig(std::vector<std::string> variableName, std::vector<std::string> polynomialList, int timesToRun) {
        #ifdef NDEBUG
            for (int i = 0; i < timesToRun; i++) {
                std::vector<std::string> basis = groebnerBasisF4(1000000007, variableName.size(), variableName, polynomialList, 32, 0);
            }
        #else
            std::vector<std::string> basis = groebnerBasisF4(1000000007, variableName.size(), variableName, polynomialList, 1, 4);
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
        amon.push_back(Monomial(Term({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(Term({0}), Rational(-1)));

        Polynomial<Rational, LexComp> a(std::move(amon));

        std::vector<Monomial<Rational>> bmon;
        bmon.push_back(Monomial(Term({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(Term({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(Term({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(Term({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> b(std::move(bmon));

        std::vector<Monomial<Rational>> cmon;
        cmon.push_back(Monomial(Term({1, 1}), Rational(1)));
        cmon.push_back(Monomial(Term({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(Term({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(Term({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> c(std::move(cmon));

        std::vector<Monomial<Rational>> dmon;
        dmon.push_back(Monomial(Term({1}), Rational(1)));
        dmon.push_back(Monomial(Term({0, 1}), Rational(1)));
        dmon.push_back(Monomial(Term({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational, LexComp> d(std::move(dmon));

        TPolynomials<Rational, LexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasis, "buchberger_cyclic4 ").call(test);
    }
    {
        std::vector<Monomial<Rational>> amon;
        amon.push_back(Monomial(Term({1, 1, 1, 1}), Rational(1)));
        amon.push_back(Monomial(Term({0}), Rational(-1)));

        Polynomial<Rational, LexComp> a(std::move(amon));

        std::vector<Monomial<Rational>> bmon;
        bmon.push_back(Monomial(Term({1, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(Term({1, 1, 0, 1}), Rational(1)));
        bmon.push_back(Monomial(Term({1, 0, 1, 1}), Rational(1)));
        bmon.push_back(Monomial(Term({0, 1, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> b(std::move(bmon));

        std::vector<Monomial<Rational>> cmon;
        cmon.push_back(Monomial(Term({1, 1}), Rational(1)));
        cmon.push_back(Monomial(Term({1, 0, 0, 1}), Rational(1)));
        cmon.push_back(Monomial(Term({0, 1, 1}), Rational(1)));
        cmon.push_back(Monomial(Term({0, 0, 1, 1}), Rational(1)));

        Polynomial<Rational, LexComp> c(std::move(cmon));

        std::vector<Monomial<Rational>> dmon;
        dmon.push_back(Monomial(Term({1}), Rational(1)));
        dmon.push_back(Monomial(Term({0, 1}), Rational(1)));
        dmon.push_back(Monomial(Term({0, 0, 1}), Rational(1)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), Rational(1)));

        Polynomial<Rational, LexComp> d(std::move(dmon));

        TPolynomials<Rational, LexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisWithCriterias, "buchberger_with_criterion_cyclic4_rational ").call(test);
    }
}


void benchmark_cyclic4_prime_field() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> amon;
        amon.push_back(Monomial(Term({1, 1, 1, 1}), PrimeField<1000000007>(1)));
        amon.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<1000000007>>> bmon;
        bmon.push_back(Monomial(Term({1, 1, 1}), PrimeField<1000000007>(1)));
        bmon.push_back(Monomial(Term({1, 1, 0, 1}), PrimeField<1000000007>(1)));
        bmon.push_back(Monomial(Term({1, 0, 1, 1}), PrimeField<1000000007>(1)));
        bmon.push_back(Monomial(Term({0, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<1000000007>>> cmon;
        cmon.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(1)));
        cmon.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(1)));
        cmon.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<1000000007>(1)));
        cmon.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<1000000007>>> dmon;
        dmon.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        dmon.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(1)));
        dmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(1)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisF4PrimeField, "f4_cyclic4 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<1000000007>> i1({
            {{{1, 1, 1, 1}}, 1},
            {{{0}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i2({
            {{{1, 1, 1}}, 1},
            {{{1, 1, 0, 1}}, 1},
            {{{1, 0, 1, 1}}, 1},
            {{{0, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i3({
            {{{1, 1}}, 1},
            {{{0, 1, 1}}, 1},
            {{{1, 0, 0, 1}}, 1},
            {{{0, 0, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i4({
            {{{1}}, 1},
            {{{0, 1}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0, 0, 0, 1}}, 1},
        });

        gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp> ideal({i1, i2, i3, i4});
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
        std::vector<Monomial<PrimeField<1000000007>>> amon;
        amon.push_back(Monomial(Term({2}), PrimeField<1000000007>(1)));
        amon.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(2)));
        amon.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(2)));
        amon.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(2)));
        amon.push_back(Monomial(Term({1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<1000000007>>> bmon;
        bmon.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(2)));
        bmon.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        bmon.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        bmon.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<1000000007>>> cmon;
        cmon.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(1)));
        cmon.push_back(Monomial(Term({1, 0, 1}), PrimeField<1000000007>(2)));
        cmon.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        cmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<1000000007>>> dmon;
        dmon.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        dmon.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(2)));
        dmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(2)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(2)));
        dmon.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisF4PrimeField, "f4_katsura4 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<1000000007>> i1({
            {{{2}}, 1},
            {{{1}}, -1},
            {{{0, 2}}, 2},
            {{{0, 0, 2}}, 2},
            {{{0, 0, 0, 2}}, 2},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i2({
            {{{1, 1}}, 2},
            {{{0, 1, 1}}, 2},
            {{{0, 1}}, -1},
            {{{0, 0, 1, 1}}, 2},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i3({
            {{{1, 0, 1}}, 2},
            {{{0, 2}}, 1},
            {{{0, 1, 0, 1}}, 2},
            {{{0, 0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i4({
            {{{1}}, 1},
            {{{0, 1}}, 2},
            {{{0, 0, 1}}, 2},
            {{{0, 0, 0, 1}}, 2},
            {{{0}}, -1},
        });

        gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp> ideal({i1, i2, i3, i4});
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
        std::vector<Monomial<PrimeField<1000000007>>> amon;
        amon.push_back(Monomial(Term({0, 1, 3}), PrimeField<1000000007>(1)));
        amon.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        amon.push_back(Monomial(Term({0}), PrimeField<1000000007>(-2)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<1000000007>>> bmon;
        bmon.push_back(Monomial(Term({3, 0, 1}), PrimeField<1000000007>(1)));
        bmon.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(1)));
        bmon.push_back(Monomial(Term({0}), PrimeField<1000000007>(-2)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<1000000007>>> cmon;
        cmon.push_back(Monomial(Term({1, 3}), PrimeField<1000000007>(1)));
        cmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(1)));
        cmon.push_back(Monomial(Term({0}), PrimeField<1000000007>(-2)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> c(std::move(cmon));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {a, b, c};
        test_time(FindGroebnerBasisF4PrimeField, "f4_sym3-3 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<1000000007>> i1({
            {{{0, 1, 3}}, 1},
            {{{1}}, 1},
            {{{0}}, -2},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i2({
            {{{3, 0, 1}}, 1},
            {{{0, 1}}, 1},
            {{{0}}, -2},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i3({
            {{{1, 3}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0}}, -2},
        });

        gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp> ideal({i1, i2, i3});
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

void benchmark_cyclic5() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({1, 1, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        std::vector<Monomial<PrimeField<1000000007>>> mon4;

        mon4.push_back(Monomial(Term({1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 1, 1, 0, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 1, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p4(std::move(mon4));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {p0, p1, p2, p3, p4};

        test_time(FindGroebnerBasisF4PrimeField, "f4_cyclic-5 ").call(test);
    }
    {
        std::vector<std::string> polCyclic6;
        polCyclic6.emplace_back("x0+x1+x2+x3+x4");
        polCyclic6.emplace_back("x0*x1+x1*x2+x2*x3+x3*x4+x0*x4");
        polCyclic6.emplace_back("x0*x1*x2+x1*x2*x3+x2*x3*x4+x0*x1*x4+x0*x3*x4");
        polCyclic6.emplace_back("x0*x1*x2*x3+x1*x2*x3*x4+x0*x2*x3*x4+x0*x1*x3*x4+x0*x1*x2*x4");
        polCyclic6.emplace_back("x0*x1*x2*x3*x4-1");
        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3", "x4"};
        test_time(FindGroebnerBasisOpenF4PrimeField, "openf4_cyclic-5 ").call(variableName, polCyclic6);
    }
}

void benchmark_cyclic6() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({1, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({1, 1, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({1, 0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        std::vector<Monomial<PrimeField<1000000007>>> mon4;

        mon4.push_back(Monomial(Term({1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 1, 1, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 1, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p4(std::move(mon4));

        std::vector<Monomial<PrimeField<1000000007>>> mon5;

        mon5.push_back(Monomial(Term({1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 1, 1, 1, 0, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 1, 1, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 1, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({0, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p5(std::move(mon5));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {p0, p1, p2, p3, p4, p5};

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_cyclic-6 ").call(test, 1);
    }

    {
        std::vector<std::string> polCyclic6;
        polCyclic6.emplace_back("x0+x1+x2+x3+x4+x5");
        polCyclic6.emplace_back("x0*x1+x1*x2+x2*x3+x3*x4+x0*x5+x4*x5");
        polCyclic6.emplace_back("x0*x1*x2+x1*x2*x3+x2*x3*x4+x0*x1*x5+x0*x4*x5+x3*x4*x5");
        polCyclic6.emplace_back("x0*x1*x2*x3+x1*x2*x3*x4+x0*x1*x2*x5+x0*x1*x4*x5+x0*x3*x4*x5+x2*x3*x4*x5");
        polCyclic6.emplace_back("x0*x1*x2*x3*x4+x0*x1*x2*x3*x5+x0*x1*x2*x4*x5+x0*x1*x3*x4*x5+x0*x2*x3*x4*x5+x1*x2*x3*x4*x5");
        polCyclic6.emplace_back("x0*x1*x2*x3*x4*x5-1");
        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3", "x4", "x5"};
        test_time(FindGroebnerBasisOpenF4PrimeFieldBig, "openf4_cyclic-6 ").call(variableName, polCyclic6, 1);
    }
}

void benchmark_cyclic7() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({1, 1, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({1, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({1, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        std::vector<Monomial<PrimeField<1000000007>>> mon4;

        mon4.push_back(Monomial(Term({1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 1, 1, 0, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 1, 0, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({1, 0, 0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 0, 0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p4(std::move(mon4));

        std::vector<Monomial<PrimeField<1000000007>>> mon5;

        mon5.push_back(Monomial(Term({1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({0, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 1, 1, 1, 0, 0, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 1, 1, 0, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 1, 0, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({1, 0, 0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon5.push_back(Monomial(Term({0, 0, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p5(std::move(mon5));

        std::vector<Monomial<PrimeField<1000000007>>> mon6;

        mon6.push_back(Monomial(Term({1, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({1, 1, 1, 1, 1, 0, 1}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({1, 1, 1, 1, 0, 1, 1}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({1, 1, 1, 0, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({1, 1, 0, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({1, 0, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({0, 1, 1, 1, 1, 1, 1}), PrimeField<1000000007>(1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p6(std::move(mon6));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {p0, p1, p2, p3, p4, p5, p6};

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_cyclic-7 ").call(test, 1);
    }
    {
        std::vector<std::string> polCyclic7;
        polCyclic7.emplace_back("x0+x1+x2+x3+x4+x5+x6");
        polCyclic7.emplace_back("x0*x1+x1*x2+x2*x3+x3*x4+x4*x5+x0*x6+x5*x6");
        polCyclic7.emplace_back("x0*x1*x2+x1*x2*x3+x2*x3*x4+x3*x4*x5+x0*x1*x6+x0*x5*x6+x4*x5*x6");
        polCyclic7.emplace_back("x0*x1*x2*x3+x1*x2*x3*x4+x2*x3*x4*x5+x0*x1*x2*x6+x0*x1*x5*x6+x0*x4*x5*x6+x3*x4*x5*x6");
        polCyclic7.emplace_back("x0*x1*x2*x3*x4+x1*x2*x3*x4*x5+x0*x1*x2*x3*x6+x0*x1*x2*x5*x6+x0*x1*x4*x5*x6+x0*x3*x4*x5*x6+x2*x3*x4*x5*x6");
        polCyclic7.emplace_back("x0*x1*x2*x3*x4*x5+x0*x1*x2*x3*x4*x6+x0*x1*x2*x3*x5*x6+x0*x1*x2*x4*x5*x6+x0*x1*x3*x4*x5*x6+x0*x2*x3*x4*x5*x6+x1*x2*x3*x4*x5*x6");
        polCyclic7.emplace_back("x0*x1*x2*x3*x4*x5*x6-1");
        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3", "x4", "x5", "x6"};
        test_time(FindGroebnerBasisOpenF4PrimeFieldBig, "openf4_cyclic-7 ").call(variableName, polCyclic7, 1);
    }
}

void benchmark_katsura9() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> monk;

        monk.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        monk.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> pk(std::move(monk));

        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({2}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        std::vector<Monomial<PrimeField<1000000007>>> mon4;

        mon4.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p4(std::move(mon4));

        std::vector<Monomial<PrimeField<1000000007>>> mon5;

        mon5.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p5(std::move(mon5));

        std::vector<Monomial<PrimeField<1000000007>>> mon6;

        mon6.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p6(std::move(mon6));

        std::vector<Monomial<PrimeField<1000000007>>> mon7;

        mon7.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p7(std::move(mon7));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {pk, p0, p1, p2, p3, p4, p5, p6, p7};

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-9 ").call(test, 1);
    }

    {
        std::vector<std::string> polKatsura9;
        polKatsura9.emplace_back("x0+2*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8-1");
        polKatsura9.emplace_back("x0^2+2*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2-x0");
        polKatsura9.emplace_back("2*x0*x1+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8-x1");
        polKatsura9.emplace_back("x1^2+2*x0*x2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8-x2");
        polKatsura9.emplace_back("2*x1*x2+2*x0*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8-x3");
        polKatsura9.emplace_back("x2^2+2*x1*x3+2*x0*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8-x4");
        polKatsura9.emplace_back("2*x2*x3+2*x1*x4+2*x0*x5+2*x1*x6+2*x2*x7+2*x3*x8-x5");
        polKatsura9.emplace_back("x3^2+2*x2*x4+2*x1*x5+2*x0*x6+2*x1*x7+2*x2*x8-x6");
        polKatsura9.emplace_back("2*x3*x4+2*x2*x5+2*x1*x6+2*x0*x7+2*x1*x8-x7");
        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"};
        test_time(FindGroebnerBasisOpenF4PrimeFieldBig, "openf4_katsura9 ").call(variableName, polKatsura9, 1);
    }

}

void benchmark_katsura10() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> monk;

        monk.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        monk.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> pk(std::move(monk));

        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({2}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        std::vector<Monomial<PrimeField<1000000007>>> mon4;

        mon4.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p4(std::move(mon4));

        std::vector<Monomial<PrimeField<1000000007>>> mon5;

        mon5.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p5(std::move(mon5));

        std::vector<Monomial<PrimeField<1000000007>>> mon6;

        mon6.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p6(std::move(mon6));

        std::vector<Monomial<PrimeField<1000000007>>> mon7;

        mon7.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p7(std::move(mon7));

        std::vector<Monomial<PrimeField<1000000007>>> mon8;

        mon8.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon8.push_back(Monomial(Term({0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p8(std::move(mon8));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {pk, p0, p1, p2, p3, p4, p5, p6, p7, p8};

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-10 ").call(test, 1);
    }

    {
        std::vector<std::string> polKatsura10;
        polKatsura10.emplace_back("x0+2*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9-1");
        polKatsura10.emplace_back("x0^2+2*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2-x0");
        polKatsura10.emplace_back("2*x0*x1+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9-x1");
        polKatsura10.emplace_back("x1^2+2*x0*x2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9-x2");
        polKatsura10.emplace_back("2*x1*x2+2*x0*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9-x3");
        polKatsura10.emplace_back("x2^2+2*x1*x3+2*x0*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9-x4");
        polKatsura10.emplace_back("2*x2*x3+2*x1*x4+2*x0*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9-x5");
        polKatsura10.emplace_back("x3^2+2*x2*x4+2*x1*x5+2*x0*x6+2*x1*x7+2*x2*x8+2*x3*x9-x6");
        polKatsura10.emplace_back("2*x3*x4+2*x2*x5+2*x1*x6+2*x0*x7+2*x1*x8+2*x2*x9-x7");
        polKatsura10.emplace_back("x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x0*x8+2*x1*x9-x8");

        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"};
        test_time(FindGroebnerBasisOpenF4PrimeFieldBig, "openf4_katsura10 ").call(variableName, polKatsura10, 1);
    }

}


void benchmark_katsura11() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> monk;

        monk.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        monk.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> pk(std::move(monk));

        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({2}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        std::vector<Monomial<PrimeField<1000000007>>> mon4;

        mon4.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p4(std::move(mon4));

        std::vector<Monomial<PrimeField<1000000007>>> mon5;

        mon5.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p5(std::move(mon5));

        std::vector<Monomial<PrimeField<1000000007>>> mon6;

        mon6.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p6(std::move(mon6));

        std::vector<Monomial<PrimeField<1000000007>>> mon7;

        mon7.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p7(std::move(mon7));

        std::vector<Monomial<PrimeField<1000000007>>> mon8;

        mon8.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon8.push_back(Monomial(Term({0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p8(std::move(mon8));

        std::vector<Monomial<PrimeField<1000000007>>> mon9;

        mon9.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p9(std::move(mon9));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {pk, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9};

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-11 ").call(test, 1);
    }

    {
        std::vector<std::string> polKatsura11;
        polKatsura11.emplace_back("x0+2*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9+2*x10-1");
        polKatsura11.emplace_back("x0^2+2*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2+2*x10^2-x0");
        polKatsura11.emplace_back("2*x0*x1+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9+2*x9*x10-x1");
        polKatsura11.emplace_back("x1^2+2*x0*x2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9+2*x8*x10-x2");
        polKatsura11.emplace_back("2*x1*x2+2*x0*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9+2*x7*x10-x3");
        polKatsura11.emplace_back("x2^2+2*x1*x3+2*x0*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9+2*x6*x10-x4");
        polKatsura11.emplace_back("2*x2*x3+2*x1*x4+2*x0*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9+2*x5*x10-x5");
        polKatsura11.emplace_back("x3^2+2*x2*x4+2*x1*x5+2*x0*x6+2*x1*x7+2*x2*x8+2*x3*x9+2*x4*x10-x6");
        polKatsura11.emplace_back("2*x3*x4+2*x2*x5+2*x1*x6+2*x0*x7+2*x1*x8+2*x2*x9+2*x3*x10-x7");
        polKatsura11.emplace_back("x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x0*x8+2*x1*x9+2*x2*x10-x8");
        polKatsura11.emplace_back("2*x4*x5+2*x3*x6+2*x2*x7+2*x1*x8+2*x0*x9+2*x1*x10-x9");
        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};
        test_time(FindGroebnerBasisOpenF4PrimeFieldBig, "openf4_katsura11 ").call(variableName, polKatsura11, 1);
    }

}


void benchmark_katsura12() {
    {
        std::vector<Monomial<PrimeField<1000000007>>> monk;

        monk.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        monk.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> pk(std::move(monk));

        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({2}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        std::vector<Monomial<PrimeField<1000000007>>> mon4;

        mon4.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(1)));
        mon4.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon4.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p4(std::move(mon4));

        std::vector<Monomial<PrimeField<1000000007>>> mon5;

        mon5.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon5.push_back(Monomial(Term({0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p5(std::move(mon5));

        std::vector<Monomial<PrimeField<1000000007>>> mon6;

        mon6.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon6.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p6(std::move(mon6));

        std::vector<Monomial<PrimeField<1000000007>>> mon7;

        mon7.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon7.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p7(std::move(mon7));

        std::vector<Monomial<PrimeField<1000000007>>> mon8;

        mon8.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon8.push_back(Monomial(Term({0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon8.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p8(std::move(mon8));

        std::vector<Monomial<PrimeField<1000000007>>> mon9;

        mon9.push_back(Monomial(Term({0, 0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon9.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p9(std::move(mon9));

        std::vector<Monomial<PrimeField<1000000007>>> mon10;

        mon10.push_back(Monomial(Term({0, 0, 0, 0, 0, 2}), PrimeField<1000000007>(1)));
        mon10.push_back(Monomial(Term({0, 0, 0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon10.push_back(Monomial(Term({0, 0, 0, 1, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon10.push_back(Monomial(Term({0, 0, 1, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon10.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon10.push_back(Monomial(Term({1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon10.push_back(Monomial(Term({0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon10.push_back(Monomial(Term({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p10(std::move(mon10));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {pk, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10};

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-12 ").call(test, 1);
    }

    {
        std::vector<std::string> polKatsura12;
        polKatsura12.emplace_back("x0+2*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9+2*x10+2*x11-1");
        polKatsura12.emplace_back("x0^2+2*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2+2*x10^2+2*x11^2-x0");
        polKatsura12.emplace_back("2*x0*x1+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9+2*x9*x10+2*x10*x11-x1");
        polKatsura12.emplace_back("x1^2+2*x0*x2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9+2*x8*x10+2*x9*x11-x2");
        polKatsura12.emplace_back("2*x1*x2+2*x0*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9+2*x7*x10+2*x8*x11-x3");
        polKatsura12.emplace_back("x2^2+2*x1*x3+2*x0*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9+2*x6*x10+2*x7*x11-x4");
        polKatsura12.emplace_back("2*x2*x3+2*x1*x4+2*x0*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9+2*x5*x10+2*x6*x11-x5");
        polKatsura12.emplace_back("x3^2+2*x2*x4+2*x1*x5+2*x0*x6+2*x1*x7+2*x2*x8+2*x3*x9+2*x4*x10+2*x5*x11-x6");
        polKatsura12.emplace_back("2*x3*x4+2*x2*x5+2*x1*x6+2*x0*x7+2*x1*x8+2*x2*x9+2*x3*x10+2*x4*x11-x7");
        polKatsura12.emplace_back("x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x0*x8+2*x1*x9+2*x2*x10+2*x3*x11-x8");
        polKatsura12.emplace_back("2*x4*x5+2*x3*x6+2*x2*x7+2*x1*x8+2*x0*x9+2*x1*x10+2*x2*x11-x9");
        polKatsura12.emplace_back("x5^2+2*x4*x6+2*x3*x7+2*x2*x8+2*x1*x9+2*x0*x10+2*x1*x11-x10");
        std::vector<std::string> variableName = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11"};
        test_time(FindGroebnerBasisOpenF4PrimeFieldBig, "openf4_katsura12 ").call(variableName, polKatsura12, 1);
    }

}
