#include "../../external/GroebnerBasisFork/GroebnerLib/includes/PolynomialSet.hpp"
#include "../../external/GroebnerBasisFork/GroebnerLib/includes/Rational.hpp"
#include "../../external/GroebnerBasisFork/GroebnerLib/includes/Modular.hpp"
#include "../../external/GroebnerBasisFork/GroebnerLib/includes/F4GB.hpp"
#include "../benchmarking.h"
#include "../../lib/algo/buchberger.h"
#include "../../lib/algo/improved_buchberger.h"
#include "../../lib/algo/f4.h"
#include "../../lib/util/rational.h"
#include "../../lib/util/prime_field.h"

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

    void FindGroebnerBasisLibModular(gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp>& ideal) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp> ideal2 = ideal;
                ideal2.MakeGroebnerBasis();
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

    void FindGroebnerBasisPrimeField(TPolynomials<PrimeField<31>, GrevLexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<PrimeField<31>, GrevLexComp> F2 = F;
                NAlgo::Buchberger::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::Buchberger::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisImprovedBuchberger(TPolynomials<Rational, LexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<Rational, LexComp> F2 = F;
                NAlgo::ImprovedBuchberger::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::ImprovedBuchberger::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisImprovedBuchbergerPrimeField(TPolynomials<PrimeField<31>, GrevLexComp>& F) {
        #ifdef NDEBUG
            for (int i = 0; i < TimesToRun; i++) {
                TPolynomials<PrimeField<31>, GrevLexComp> F2 = F;
                NAlgo::ImprovedBuchberger::FindGroebnerBasis(F2);
            }
        #else
            NAlgo::ImprovedBuchberger::FindGroebnerBasis(F);
            std::cout << F << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebner(F));
        #endif
    }

    void FindGroebnerBasisImprovedBuchbergerPrimeFieldBig(TPolynomials<PrimeField<1000000007>, GrevLexComp>& F) {
        #ifdef NDEBUG
            NAlgo::ImprovedBuchberger::FindGroebnerBasis(F);
        #else
            NAlgo::ImprovedBuchberger::FindGroebnerBasis(F);
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


    void FindGroebnerBasisF4PrimeFieldBig(TPolynomials<PrimeField<1000000007>, GrevLexComp>& F) {
        #ifdef NDEBUG
            NAlgo::F4::FindGroebnerBasis(F, 32);
        #else
            NAlgo::F4::FindGroebnerBasis(F);
            std::cout << F.size() << std::endl;
            assert(NAlgo::NUtil::CheckBasisIsGroebnerBig(F));
        #endif
    }


    void FindGroebnerBasisF4LibModularBig(gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp>& ideal) {
        #ifdef NDEBUG
            gb::inplace_calculate_f4_gb(ideal);
        #else
            gb::inplace_calculate_f4_gb(ideal);
            std::cout << ideal << std::endl;
        #endif
    }
}

void benchmark_buchberger_cyclic4_rational() {
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
        test_time(FindGroebnerBasis, "buchberger_cyclic4_rational ").call(test);
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
        test_time(FindGroebnerBasisImprovedBuchberger, "improved_buchberger_cyclic4_rational ").call(test);
    }
}


void benchmark_buchberger_cyclic4() {
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(Term({1, 1, 1, 1}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(Term({1, 1, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({1, 1, 0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({1, 0, 1, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({0, 1, 1, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(Term({1, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<31>>> dmon;
        dmon.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisPrimeField, "buchberger_cyclic4 ").call(test);
    }
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(Term({1, 1, 1, 1}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(Term({1, 1, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({1, 1, 0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({1, 0, 1, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({0, 1, 1, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(Term({1, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<31>>> dmon;
        dmon.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(1)));

        Polynomial<PrimeField<31>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c, d};
        test_time(FindGroebnerBasisImprovedBuchbergerPrimeField, "improved_buchberger_cyclic4 ").call(test);
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
        test_time(FindGroebnerBasisLibModular, "GroebnerBasisLibBuchberger_cyclic4 ").call(ideal);
    }
}

void benchmark_buchberger_katsura4() {
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(Term({2}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(Term({0, 0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(Term({1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(Term({1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(Term({0, 2}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<31>>> dmon;
        dmon.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(Term({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c, d};

        test_time(FindGroebnerBasisImprovedBuchbergerPrimeField, "improved_buchberger_katsura4 ").call(test);
    }
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(Term({2}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(Term({0, 0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<31>(2)));
        amon.push_back(Monomial(Term({1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(Term({1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<31>(2)));
        bmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(Term({0, 2}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<31>(2)));
        cmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        std::vector<Monomial<PrimeField<31>>> dmon;
        dmon.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        dmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(2)));
        dmon.push_back(Monomial(Term({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> d(std::move(dmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c, d};

        test_time(FindGroebnerBasisPrimeField, "buchberger_katsura4 ").call(test);
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
        test_time(FindGroebnerBasisLibModular, "GroebnerBasisLibBuchberger_katsura4 ").call(ideal);
    }
}

void benchmark_buchberger_sym3_3() {
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(Term({0, 1, 3}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(Term({3, 0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(Term({1, 3}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c};

        test_time(FindGroebnerBasisImprovedBuchbergerPrimeField, "improved_buchberger_sym3_3 ").call(test);
    }
    {
        std::vector<Monomial<PrimeField<31>>> amon;
        amon.push_back(Monomial(Term({0, 1, 3}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        amon.push_back(Monomial(Term({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> a(std::move(amon));

        std::vector<Monomial<PrimeField<31>>> bmon;
        bmon.push_back(Monomial(Term({3, 0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({0, 1}), PrimeField<31>(1)));
        bmon.push_back(Monomial(Term({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> b(std::move(bmon));

        std::vector<Monomial<PrimeField<31>>> cmon;
        cmon.push_back(Monomial(Term({1, 3}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(Term({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c};

        test_time(FindGroebnerBasisPrimeField, "buchberger_sym3_3 ").call(test);
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
        test_time(FindGroebnerBasisLibModular, "GroebnerBasisLibBuchberger_sym3_3 ").call(ideal);
    }
}

void benchmark_buchberger_katsura5(){
    {
        std::vector<Monomial<PrimeField<31>>> monk;

        monk.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        monk.push_back(Monomial(Term({0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> pk(std::move(monk));

        std::vector<Monomial<PrimeField<31>>> mon0;

        mon0.push_back(Monomial(Term({2}), PrimeField<31>(1)));
        mon0.push_back(Monomial(Term({0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({0, 0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<31>>> mon1;

        mon1.push_back(Monomial(Term({1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<31>>> mon2;

        mon2.push_back(Monomial(Term({0, 2}), PrimeField<31>(1)));
        mon2.push_back(Monomial(Term({1, 0, 1}), PrimeField<31>(2)));
        mon2.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<31>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<31>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<31>>> mon3;

        mon3.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(2)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<31>(2)));
        mon3.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<31>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p3(std::move(mon3));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {pk, p0, p1, p2, p3};

        test_time(FindGroebnerBasisImprovedBuchbergerPrimeField, "improved_buchberger_katsura5 ").call(test);
    }
    {
        std::vector<Monomial<PrimeField<31>>> monk;

        monk.push_back(Monomial(Term({1}), PrimeField<31>(1)));
        monk.push_back(Monomial(Term({0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<31>(2)));
        monk.push_back(Monomial(Term({0}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> pk(std::move(monk));

        std::vector<Monomial<PrimeField<31>>> mon0;

        mon0.push_back(Monomial(Term({2}), PrimeField<31>(1)));
        mon0.push_back(Monomial(Term({0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({0, 0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<31>(2)));
        mon0.push_back(Monomial(Term({1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<31>>> mon1;

        mon1.push_back(Monomial(Term({1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<31>(2)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<31>>> mon2;

        mon2.push_back(Monomial(Term({0, 2}), PrimeField<31>(1)));
        mon2.push_back(Monomial(Term({1, 0, 1}), PrimeField<31>(2)));
        mon2.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<31>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<31>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<31>>> mon3;

        mon3.push_back(Monomial(Term({0, 1, 1}), PrimeField<31>(2)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<31>(2)));
        mon3.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<31>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<31>(-1)));

        Polynomial<PrimeField<31>, GrevLexComp> p3(std::move(mon3));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {pk, p0, p1, p2, p3};

        test_time(FindGroebnerBasisPrimeField, "buchberger_katsura5 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<31>> i1({
            {{{2}}, 1},
            {{{0, 2}}, 2},
            {{{0, 0, 2}}, 2},
            {{{0, 0, 0, 2}}, 2},
            {{{0, 0, 0, 0, 2}}, 2},
            {{{1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i2({
            {{{1, 1}}, 2},
            {{{0, 1, 1}}, 2},
            {{{0, 0, 1, 1}}, 2},
            {{{0, 0, 0, 1, 1}}, 2},
            {{{0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i3({
            {{{1, 0, 1}}, 2},
            {{{0, 2}}, 1},
            {{{0, 1, 0, 1}}, 2},
            {{{0, 0, 1, 0, 1}}, 2},
            {{{0, 0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i4({
            {{{0, 1, 1}}, 2},
            {{{1, 0, 0, 1}}, 2},
            {{{0, 1, 0, 0, 1}}, 2},
            {{{0, 0, 0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<31>> i5({
            {{{1}}, 1},
            {{{0, 1}}, 2},
            {{{0, 0, 1}}, 2},
            {{{0, 0, 0, 1}}, 2},
            {{{0, 0, 0, 0, 1}}, 2},
            {{{0}}, -1},
        });

        gb::PolynomialSet<gb::fields::Modular<31>, gb::DegReLexComp> ideal({i1, i2, i3, i4, i5});
        test_time(FindGroebnerBasisLibModular, "GroebnerBasisLibBuchberger_katsura5 ").call(ideal);
    }
}


void benchmark_cyclic4() {
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

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_cyclic-5 ").call(test);
    }

    {
        gb::Polynomial<gb::fields::Modular<1000000007>> i1({
            {{{1, 1, 1, 1, 1}}, 1},
            {{{0}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i2({
            {{{1}}, 1},
            {{{0, 1}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0, 0, 0, 1}}, 1},
            {{{0, 0, 0, 0, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i3({
            {{{1, 1}}, 1},
            {{{0, 1, 1}}, 1},
            {{{0, 0, 1, 1}}, 1},
            {{{1, 0, 0, 0, 1}}, 1},
            {{{0, 0, 0, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i4({
            {{{1, 1, 1}}, 1},
            {{{0, 1, 1, 1}}, 1},
            {{{1, 1, 0, 0, 1}}, 1},
            {{{1, 0, 0, 1, 1}}, 1},
            {{{0, 0, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i5({
            {{{1, 1, 1, 1}}, 1},
            {{{1, 1, 1, 0, 1}}, 1},
            {{{1, 1, 0, 1, 1}}, 1},
            {{{1, 0, 1, 1, 1}}, 1},
            {{{0, 1, 1, 1, 1}}, 1},
        });

        gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp> ideal({i1, i2, i3, i4, i5});
        test_time(FindGroebnerBasisF4LibModularBig, "GroebnerBasisLibF4_cyclic5 ").call(ideal);
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

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_cyclic-6 ").call(test);
    }

    {
        gb::Polynomial<gb::fields::Modular<1000000007>> i1({
            {{{1, 1, 1, 1, 1, 1}}, 1},
            {{{0}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i2({
            {{{1}}, 1},
            {{{0, 1}}, 1},
            {{{0, 0, 1}}, 1},
            {{{0, 0, 0, 1}}, 1},
            {{{0, 0, 0, 0, 1}}, 1},
            {{{0, 0, 0, 0, 0, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i3({
            {{{1, 1}}, 1},
            {{{0, 1, 1}}, 1},
            {{{0, 0, 1, 1}}, 1},
            {{{0, 0, 0, 1, 1}}, 1},
            {{{1, 0, 0, 0, 0, 1}}, 1},
            {{{0, 0, 0, 0, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i4({
            {{{1, 1, 1}}, 1},
            {{{0, 1, 1, 1}}, 1},
            {{{0, 0, 1, 1, 1}}, 1},
            {{{1, 1, 0, 0, 0, 1}}, 1},
            {{{1, 0, 0, 0, 1, 1}}, 1},
            {{{0, 0, 0, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i5({
            {{{1, 1, 1, 1}}, 1},
            {{{0, 1, 1, 1, 1}}, 1},
            {{{1, 1, 1, 0, 0, 1}}, 1},
            {{{1, 1, 0, 0, 1, 1}}, 1},
            {{{1, 0, 0, 1, 1, 1}}, 1},
            {{{0, 0, 1, 1, 1, 1}}, 1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i6({
            {{{1, 1, 1, 1, 1}}, 1},
            {{{1, 1, 1, 1, 0, 1}}, 1},
            {{{1, 1, 1, 0, 1, 1}}, 1},
            {{{1, 1, 0, 1, 1, 1}}, 1},
            {{{1, 0, 1, 1, 1, 1}}, 1},
            {{{0, 1, 1, 1, 1, 1}}, 1},
        });

        gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp> ideal({i1, i2, i3, i4, i5, i6});
        test_time(FindGroebnerBasisF4LibModularBig, "GroebnerBasisLibF4_cyclic6 ").call(ideal);
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

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_cyclic-7 ").call(test);
    }
}

void benchmark_katsura5(){
    {
        std::vector<Monomial<PrimeField<1000000007>>> monk;

        monk.push_back(Monomial(Term({1}), PrimeField<1000000007>(1)));
        monk.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0, 0, 0, 0, 1}), PrimeField<1000000007>(2)));
        monk.push_back(Monomial(Term({0}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> pk(std::move(monk));

        std::vector<Monomial<PrimeField<1000000007>>> mon0;

        mon0.push_back(Monomial(Term({2}), PrimeField<1000000007>(1)));
        mon0.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({0, 0, 0, 0, 2}), PrimeField<1000000007>(2)));
        mon0.push_back(Monomial(Term({1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p0(std::move(mon0));

        std::vector<Monomial<PrimeField<1000000007>>> mon1;

        mon1.push_back(Monomial(Term({1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 0, 0, 1, 1}), PrimeField<1000000007>(2)));
        mon1.push_back(Monomial(Term({0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p1(std::move(mon1));

        std::vector<Monomial<PrimeField<1000000007>>> mon2;

        mon2.push_back(Monomial(Term({0, 2}), PrimeField<1000000007>(1)));
        mon2.push_back(Monomial(Term({1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1, 0, 1}), PrimeField<1000000007>(2)));
        mon2.push_back(Monomial(Term({0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p2(std::move(mon2));

        std::vector<Monomial<PrimeField<1000000007>>> mon3;

        mon3.push_back(Monomial(Term({0, 1, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 1, 0, 0, 1}), PrimeField<1000000007>(2)));
        mon3.push_back(Monomial(Term({0, 0, 0, 1}), PrimeField<1000000007>(-1)));

        Polynomial<PrimeField<1000000007>, GrevLexComp> p3(std::move(mon3));

        TPolynomials<PrimeField<1000000007>, GrevLexComp> test = {pk, p0, p1, p2, p3};

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura5 ").call(test);
    }
    {
        gb::Polynomial<gb::fields::Modular<1000000007>> i1({
            {{{2}}, 1},
            {{{0, 2}}, 2},
            {{{0, 0, 2}}, 2},
            {{{0, 0, 0, 2}}, 2},
            {{{0, 0, 0, 0, 2}}, 2},
            {{{1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i2({
            {{{1, 1}}, 2},
            {{{0, 1, 1}}, 2},
            {{{0, 0, 1, 1}}, 2},
            {{{0, 0, 0, 1, 1}}, 2},
            {{{0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i3({
            {{{1, 0, 1}}, 2},
            {{{0, 2}}, 1},
            {{{0, 1, 0, 1}}, 2},
            {{{0, 0, 1, 0, 1}}, 2},
            {{{0, 0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i4({
            {{{0, 1, 1}}, 2},
            {{{1, 0, 0, 1}}, 2},
            {{{0, 1, 0, 0, 1}}, 2},
            {{{0, 0, 0, 1}}, -1},
        });

        gb::Polynomial<gb::fields::Modular<1000000007>> i5({
            {{{1}}, 1},
            {{{0, 1}}, 2},
            {{{0, 0, 1}}, 2},
            {{{0, 0, 0, 1}}, 2},
            {{{0, 0, 0, 0, 1}}, 2},
            {{{0}}, -1},
        });

        gb::PolynomialSet<gb::fields::Modular<1000000007>, gb::DegReLexComp> ideal({i1, i2, i3, i4, i5});
        test_time(FindGroebnerBasisF4LibModularBig, "GroebnerBasisLibF4_katsura5 ").call(ideal);
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

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-9 ").call(test);
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

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-10 ").call(test);
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

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-11 ").call(test);
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

        test_time(FindGroebnerBasisF4PrimeFieldBig, "f4_katsura-12 ").call(test);
    }
}
