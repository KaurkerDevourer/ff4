#pragma once
#include "util/groebner_basis_util.h"

namespace FF4 {
    namespace NAlgo {
        namespace Buchberger {
            template <typename TCoef, typename TComp>
            void FindGroebnerBasis(NUtils::TPolynomials<TCoef, TComp>& F) {
                std::queue<std::pair<size_t, size_t> > pairs_to_check = NUtil::GetPairsToCheck(F.size());

                while(!pairs_to_check.empty()) {
                    const NUtils::Polynomial<TCoef, TComp>& fi = F[pairs_to_check.front().first];
                    const NUtils::Polynomial<TCoef, TComp>& fj = F[pairs_to_check.front().second];
                    pairs_to_check.pop();
                    const NUtils::Monomial<TCoef>& gi = fi.GetLeadingMonomial();
                    const NUtils::Monomial<TCoef>& gj = fj.GetLeadingMonomial();
                    NUtils::Monomial<TCoef> glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                    NUtils::Polynomial<TCoef, TComp> S = fi * (glcm / gi) - fj * (glcm / gj);
                    if (!NUtil::ReduceToZero(S, F)) {
                        for (size_t i = 0; i < F.size(); i++) {
                            pairs_to_check.push({i, F.size()});
                        }
                        F.push_back(S);
                    }
                }
            }
        }
    }
}
