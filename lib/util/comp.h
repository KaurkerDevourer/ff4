#pragma once

#include "monomial.h"

namespace NUtils {
    class LexComp {
    public:
        template <typename T>
        bool operator()(const Monomial<T>& left, const Monomial<T>& right) const noexcept {
            assert(left.GetCoef() != 0);
            assert(right.GetCoef() != 0);
            return left.GetTerm() < right.GetTerm();
        }
    };
}