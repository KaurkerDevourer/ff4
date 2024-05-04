#pragma once

#include "monomial.hpp"
#include "critical_pair.hpp"

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

    class DexComp {
    public:
        template <typename T>
        bool operator()(const CriticalPair<T>& left, const CriticalPair<T>& right) const noexcept {
            assert(left.GetGlcm().GetCoef() != 0);
            assert(right.GetGlcm().GetCoef() != 0);
            if (left.GetGlcm().GetTerm().GetDegree() != right.GetGlcm().GetTerm().GetDegree()) {
                return left.GetGlcm().GetTerm().GetDegree() < right.GetGlcm().GetTerm().GetDegree();
            }
            return left.GetGlcm().GetTerm() < right.GetGlcm().GetTerm();
        }
    };

    class TTermReverseComp {
    public:
        bool operator()(const TTerm& left, const TTerm& right) const noexcept {
            return right < left;
        }
    };
}
