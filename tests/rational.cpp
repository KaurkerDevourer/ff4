#include "../lib/util/rational.h"
#include "testing.h"
#include <iostream>
#include <cassert>

void test_rational() {
    using namespace FF4::NUtils;
    Rational r(3, 6);
    ASSERT_EQUAL(r.GetNumerator(), (int64_t)1);
    ASSERT_EQUAL(r.GetDenominator(), (int64_t)2);

    r *= 2;
    ASSERT_EQUAL(r.GetNumerator(), (int64_t)1);
    ASSERT_EQUAL(r.GetDenominator(), (int64_t)1);

    Rational r2(3, 4);

    r *= r2;
    ASSERT_EQUAL(r, r2);

    r += r2;
    ASSERT_EQUAL(r, 2 * r2);

    Rational k(1, 2);
    Rational t(1, 3);
    Rational p(1, 6);
    Rational p2(5, 6);
    ASSERT_EQUAL(k - t, p);
    ASSERT_EQUAL(k + t, p2);
    ASSERT_EQUAL(k * t, p);
    ASSERT_EQUAL(k / t, Rational(3, 2));

    ASSERT_EQUAL(-k - t, -p2);

    std::cout << "Successfully tested Rational" << std::endl;
}
