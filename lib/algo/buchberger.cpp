#include "buchberger.h"
#include <queue>

namespace NAlgo {
    namespace Buchberger {

        void DoProcess(NUtils::TPolynomials& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check;
            for (size_t i = 0; i < F.size(); i++) {
                for (size_t j = i + 1; j < F.size(); j++) {
                    pairs_to_check.push({i, j});
                }
            }
            while(pairs_to_check.size()) {
                const NUtils::Polynomial& fi = F[pairs_to_check.front().first];
                const NUtils::Polynomial& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                const NUtils::Monomial& gi = fi.GetHeadMonomial();
                const NUtils::Monomial& gj = fj.GetHeadMonomial();
                NUtils::Monomial glcm = lcm(gi, gj);
                NUtils::Polynomial S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!S.ReduceBy(F)) {
                    for (size_t i = 0; i < F.size(); i++) {
                        pairs_to_check.push({i, F.size()});
                    }
                    F.push_back(S);
                }
            }
        }
    }
}