#pragma once

#include "monomial.h"
#include "critical_pair.h"

namespace FF4 {
    namespace NUtils {
        class LexComp {
        public:
            bool operator()(const TTerm& left, const TTerm right) const noexcept {
                return left < right;
            }

            template <typename T>
            bool operator()(const Monomial<T>& left, const Monomial<T>& right) const noexcept {
                assert(left.GetCoef() != 0);
                assert(right.GetCoef() != 0);
                return LexComp()(left.GetTerm(), right.GetTerm());
            }

            template <typename T>
            bool operator()(const CriticalPair<T, LexComp>& left, const CriticalPair<T, LexComp>& right) const noexcept {
                return LexComp()(left.GetGlcm(), right.GetGlcm());
            }
        };

        class RevLexComp {
        public:
            bool operator()(const TTerm& left, const TTerm& right) const noexcept {
                if (left.size() != right.size()) {
                    return left.size() > right.size();
                }
                for (int i = left.size() - 1; i >= 0; i--) {
                    if (left[i] != right[i]) {
                        return left[i] > right[i];
                    }
                }
                return false;
            }

            template <typename T>
            bool operator()(const Monomial<T>& left, const Monomial<T>& right) const noexcept {
                assert(left.GetCoef() != 0);
                assert(right.GetCoef() != 0);
                const TTerm& l = left.GetTerm();
                const TTerm& r = right.GetTerm();
                return RevLexComp()(l, r);
            }
        };

        class GrevLexComp {
        public:

            bool operator()(const TTerm& left, const TTerm& right) const noexcept {
                if (left.TotalDegree() != right.TotalDegree()) {
                    return left.TotalDegree() < right.TotalDegree();
                }
                return RevLexComp()(left, right);
            }

            template <typename T>
            bool operator()(const Monomial<T>& left, const Monomial<T>& right) const noexcept {
                assert(left.GetCoef() != 0);
                assert(right.GetCoef() != 0);
                return GrevLexComp()(left.GetTerm(), right.GetTerm());
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
}
