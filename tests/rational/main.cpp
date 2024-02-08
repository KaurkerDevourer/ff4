#include "../../lib/util/rational.h"
#include "../testing.h"
#include <iostream>
#include <cassert>

int main() {
    NUtils::Rational r(3, 6);
    ASSERT_EQUAL(r.GetNumerator(), (int64_t)1);
    ASSERT_EQUAL(r.GetDenominator(), (int64_t)2);

    r *= 2;
    ASSERT_EQUAL(r.GetNumerator(), (int64_t)1);
    ASSERT_EQUAL(r.GetDenominator(), (int64_t)1);

    NUtils::Rational r2(3, 4);

    r *= r2;
    ASSERT_EQUAL(r, r2);

    r += r2;
    ASSERT_EQUAL(r, 2 * r2);

    NUtils::Rational k(1, 2);
    NUtils::Rational t(1, 3);
    NUtils::Rational p(1, 6);
    NUtils::Rational p2(5, 6);
    ASSERT_EQUAL(k - t, p);
    ASSERT_EQUAL(k + t, p2);
    ASSERT_EQUAL(k * t, p);
    ASSERT_EQUAL(k / t, NUtils::Rational(3, 2));

    std::cout << "Successfully tested Rational" << std::endl;
}