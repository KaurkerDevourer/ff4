#include <iostream>

#include "../lib/algo/f4.h"
#include "../lib/util/rational.h"
#include "../lib/util/prime_field.hpp"
#include "../lib/algo/util/groebner_basis_util.h"

void test_f4() {
    using namespace FF4::NUtils;
    // cyclic-4
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
        std::cout << "F4: " << test << std::endl;
        FF4::NAlgo::F4::FindGroebnerBasis(test);
        std::cout << "Size of Groebner basis by F4: " << test.size() << std::endl;
        std::cout << test << std::endl;
        assert(FF4::NAlgo::NUtil::CheckBasisIsGroebner(test));
    }
    // cyclic-4-prime-field
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
        std::cout << "F4: " << test << std::endl;
        FF4::NAlgo::F4::FindGroebnerBasis(test);
        std::cout << "Size of Groebner basis by F4: " << test.size() << std::endl;
        std::cout << test << std::endl;
        assert(FF4::NAlgo::NUtil::CheckBasisIsGroebner(test));
    }

    // katsura4
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
        std::cout << "F4: " << test << std::endl;
        FF4::NAlgo::F4::FindGroebnerBasis(test);
        std::cout << "Size of Groebner basis by F4: " << test.size() << std::endl;
        std::cout << test << std::endl;
        assert(FF4::NAlgo::NUtil::CheckBasisIsGroebner(test));
    }

    // sym3-3
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
        cmon.push_back(Monomial(TTerm({1, 3, 0}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0, 0, 1}), PrimeField<31>(1)));
        cmon.push_back(Monomial(TTerm({0}), PrimeField<31>(-2)));

        Polynomial<PrimeField<31>, GrevLexComp> c(std::move(cmon));

        TPolynomials<PrimeField<31>, GrevLexComp> test = {a, b, c};
        std::cout << "F4: " << test << std::endl;
        FF4::NAlgo::F4::FindGroebnerBasis(test);
        std::cout << "Size of Groebner basis by F4: " << test.size() << std::endl;
        std::cout << test << std::endl;
        assert(FF4::NAlgo::NUtil::CheckBasisIsGroebner(test));
    }
}
