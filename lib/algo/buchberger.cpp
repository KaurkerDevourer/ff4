#include "buchberger.h"
#include <queue>

namespace NAlgo {
    namespace Buchberger {

        void FindGroebnerBasis(NUtils::TPolynomials& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check;
            for (size_t i = 0; i < F.size(); i++) {
                for (size_t j = i + 1; j < F.size(); j++) {
                    pairs_to_check.push({i, j});
                }
            }
            while(!pairs_to_check.empty()) {
                const NUtils::Polynomial& fi = F[pairs_to_check.front().first];
                const NUtils::Polynomial& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                const NUtils::Monomial& gi = fi.GetHeadMonomial();
                const NUtils::Monomial& gj = fj.GetHeadMonomial();
                NUtils::Monomial glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), 1);
                NUtils::Polynomial S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!ReduceToZero(S, F)) {
                    for (size_t i = 0; i < F.size(); i++) {
                        pairs_to_check.push({i, F.size()});
                    }
                    F.push_back(S);
                }
            }
        }

        bool ReduceToZero(NUtils::Polynomial& F, NUtils::TPolynomials& polynomialsSet) {
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
    }
}
