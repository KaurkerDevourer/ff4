#pragma once

#include "monomial.hpp"
#include "critical_pair.hpp"

namespace NUtils {
    class LexComp {
    public:
        template <typename T>
        bool operator()(const CriticalPair<T, LexComp>& left, const CriticalPair<T, LexComp>& right) const noexcept {
            assert(left.GetGlcm().GetCoef() != 0);
            assert(right.GetGlcm().GetCoef() != 0);
            return left.GetGlcm().GetTerm() < right.GetGlcm().GetTerm();
        }

        template <typename T>
        bool operator()(const Monomial<T>& left, const Monomial<T>& right) const noexcept {
            assert(left.GetCoef() != 0);
            assert(right.GetCoef() != 0);
            return left.GetTerm() < right.GetTerm();
        }
    };

    class RevLexComp {
    public:
        template <typename T>
        bool operator()(const Monomial<T>& left, const Monomial<T>& right) const noexcept {
            assert(left.GetCoef() != 0);
            assert(right.GetCoef() != 0);
            const TTerm& l = left.GetTerm();
            const TTerm& r = right.GetTerm();
            if (l.size() != r.size()) {
                return l.size() > r.size();
            }
            for (int i = l.size() - 1; i >= 0; i--) {
                if (l[i] != r[i]) {
                    return l[i] > r[i];
                }
            }
            return false;
        }
    };

    class GrevLexComp {
    public:
        template <typename T>
        bool operator()(const Monomial<T>& left, const Monomial<T>& right) const noexcept {
            assert(left.GetCoef() != 0);
            assert(right.GetCoef() != 0);
            if (left.GetTerm().GetDegree() != right.GetTerm().GetDegree()) {
                return left.GetTerm().GetDegree() < right.GetTerm().GetDegree();
            }
            return RevLexComp()(left, right);
        }

        template <typename T>
        bool operator()(const CriticalPair<T, GrevLexComp>& left, const CriticalPair<T, GrevLexComp>& right) const noexcept {
            return GrevLexComp()(left.GetGlcm(), right.GetGlcm());
        }
    };

    class TTermReverseComp {
    public:
        bool operator()(const TTerm& left, const TTerm& right) const noexcept {
            return right < left;
        }
    };
}
