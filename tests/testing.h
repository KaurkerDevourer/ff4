#pragma once
#include <iostream>
#include <cassert>

#define ASSERT_EQUAL(A, B) \
    if (A != B) { \
        std::cout << A << " != " << B << std::endl; \
        assert(false); \
    }
