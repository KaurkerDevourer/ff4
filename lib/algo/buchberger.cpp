#include "buchberger.h"

namespace NAlgo {
    namespace Buchberger {

        void DoProcess(NUtils::TPolynomials& F) {
            bool ready = false;
            while(!ready) {
                ready = true;
                for (size_t i = 0; i < F.size(); i++) {
                    for (size_t j = i + 1; j < F.size(); j++) {
                        const NUtils::Monomial& gi = F[i].GetHeadMonomial();
                        const NUtils::Monomial& gj = F[j].GetHeadMonomial();
                        NUtils::Monomial glcm = lcm(gi, gj);
                        NUtils::Polynomial S = F[i] * (glcm/gi) - F[j] * (glcm/gj);
                        if (!S.ReduceBy(F)) {
                            ready = false;
                            F.push_back(S);
                        }
                    }
                }
            }
        }
    }
}