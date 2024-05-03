#pragma once
#include "util/groebner_basis_util.h"

namespace NAlgo {
    namespace Buchberger {
        using namespace NUtils;
        template <typename TCoef>
        void FindGroebnerBasis(TPolynomials<TCoef>& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheck(F.size());

            while(!pairs_to_check.empty()) {
                const Polynomial<TCoef>& fi = F[pairs_to_check.front().first];
                const Polynomial<TCoef>& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                const Monomial<TCoef>& gi = fi.GetHeadMonomial();
                const Monomial<TCoef>& gj = fj.GetHeadMonomial();
                Monomial<TCoef> glcm = Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                Polynomial<TCoef> S = fi * (glcm/gi) - fj * (glcm/gj);
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
