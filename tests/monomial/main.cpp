#include "../../lib/util/monomial.h"
#include "../testing.h"
#include <iostream>
#include <cassert>

int main() {
    using namespace NUtils;
    Monomial a(TTerm({0, 2}), 2);
    Monomial b(TTerm({1, 0}), 3);

    ASSERT_EQUAL(a * b, Monomial(TTerm({1, 2}), 6));
    ASSERT_EQUAL((a * b) / b, a);
    ASSERT_EQUAL((a * b) / a, b);
    ASSERT_EQUAL((a * b), (b * a));
    ASSERT_EQUAL((a < b), true);
    ASSERT_EQUAL((b > a), true);
    ASSERT_EQUAL((b >= a), true);
    ASSERT_EQUAL((a > b), false);

    std::cout << "Successfully tested Monomial" << std::endl;
}
