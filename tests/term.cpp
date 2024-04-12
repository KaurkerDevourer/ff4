#include "../lib/util/term.h"
#include "testing.h"
#include <iostream>
#include <cassert>

void test_term() {
    using namespace NUtils;
    TTerm a({1, 2, 3});
    TTerm b({3, 2, 1});
    ASSERT_EQUAL(lcm(a, b), TTerm({3, 2, 3}));
    ASSERT_EQUAL(gcd(a, b), TTerm({1, 2, 1}));
    ASSERT_EQUAL(a.IsDivisibleBy(b), false);
    ASSERT_EQUAL(b.IsDivisibleBy(a), false);
    a = TTerm({2, 4, 6});
    b = TTerm({0, 3, 5});
    ASSERT_EQUAL(lcm(a, b), a);
    ASSERT_EQUAL(gcd(a, b), b);
    ASSERT_EQUAL(a.IsDivisibleBy(b), true);
    ASSERT_EQUAL(b.IsDivisibleBy(a), false);

    ASSERT_EQUAL(a * b, TTerm({2, 7, 11}));
    ASSERT_EQUAL(a / b, TTerm({2, 1, 1}));


    std::cout << "Successfully tested Term" << std::endl;
}
