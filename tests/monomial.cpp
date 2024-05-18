#include "../lib/util/monomial.h"
#include "../lib/util/comp.h"
#include "testing.h"
#include <iostream>
#include <cassert>

void test_monomial() {
    using namespace FF4::NUtils;
    Monomial a(Term({0, 2}), 2);
    Monomial b(Term({1}), 3);

    ASSERT_EQUAL(a * b, Monomial(Term({1, 2}), 6));
    ASSERT_EQUAL((a * b) / b, a);
    ASSERT_EQUAL((a * b) / a, b);
    ASSERT_EQUAL((a * b), (b * a));
    ASSERT_EQUAL(LexComp()(a, b), true);
    ASSERT_EQUAL(LexComp()(b, a), false);

    std::cout << "Successfully tested Monomial" << std::endl;
}
