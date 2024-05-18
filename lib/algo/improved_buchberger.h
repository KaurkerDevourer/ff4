#pragma once
#include "util/groebner_basis_util.h"
#include <algorithm>
#include <cassert>

namespace FF4 {
    namespace NAlgo {
        namespace ImprovedBuchberger {
            template <typename TCoef, typename TComp>
            void FindGroebnerBasis(NUtils::TPolynomials<TCoef, TComp>& F) {
                NUtil::TPolynomialSet<TCoef, TComp> polynomials;
                NUtil::TPairsSet<TCoef, TComp> pairs_to_check;
                for (auto& f : F) {
                    NUtil::UpdateCriticalPairs(polynomials, pairs_to_check, f);
                }

                while(!pairs_to_check.empty()) {
                    NUtils::CriticalPair<TCoef, TComp> cp = (*pairs_to_check.begin());
                    pairs_to_check.erase(pairs_to_check.begin());

                    NUtils::Polynomial<TCoef, TComp> S = cp.GetGlcm() / cp.GetLeftTerm() * cp.GetLeft() - cp.GetGlcm() / cp.GetRightTerm() * cp.GetRight();

                    if (!NUtil::InplaceReduceToZero(S, polynomials)) {
                        NUtil::UpdateCriticalPairs(polynomials, pairs_to_check, S);
                    }
                }
                NUtil::UpdateBasis(polynomials, F);
            }
        }
    }
}
