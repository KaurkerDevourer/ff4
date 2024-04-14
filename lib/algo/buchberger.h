#pragma once
#include "../util/polynomial.h"

namespace NAlgo {
    namespace Buchberger {
        template <typename TCoef>
        bool ReduceToZero(NUtils::Polynomial<TCoef>& F, NUtils::TPolynomials<TCoef>& polynomialsSet) {
            if (F.IsZero()) {
                return true;
            }
            bool changed = true;
            while(changed && !F.IsZero()) {
                changed = false;
                for (const auto& f : polynomialsSet) {
                    while (!F.IsZero() && F.GetHeadMonomial().GetTerm().IsDivisibleBy(f.GetHeadMonomial().GetTerm())) {
                        F -= f * (F.GetHeadMonomial() / f.GetHeadMonomial());
                        changed = true;
                    }
                }
            }
            return F.IsZero();
        }

        template <typename TCoef>
        void FindGroebnerBasis(NUtils::TPolynomials<TCoef>& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = NUtils::GetPairsToCheck(F.size());

            while(!pairs_to_check.empty()) {
                const NUtils::Polynomial<TCoef>& fi = F[pairs_to_check.front().first];
                const NUtils::Polynomial<TCoef>& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                const NUtils::Monomial<TCoef>& gi = fi.GetHeadMonomial();
                const NUtils::Monomial<TCoef>& gj = fj.GetHeadMonomial();
                NUtils::Monomial<TCoef> glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                NUtils::Polynomial<TCoef> S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!ReduceToZero(S, F)) {
                    for (size_t i = 0; i < F.size(); i++) {
                        pairs_to_check.push({i, F.size()});
                    }
                    F.push_back(S);
                }
            }
        }
    }
}
