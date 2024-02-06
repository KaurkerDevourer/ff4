#include "buchberger.h"

namespace NAlgo {
    namespace Buchberger {

        void DoProcess(TPolynomials& F) {
            bool ready = false;
            while(!ready) {
                ready = true;
                for (size_t i = 0; i < F.size(); i++) {
                    for (size_t j = i + 1; j < F.size(); j++) {
                        Monomial& gi = F[i].GetHeadMonomial();
                        Monomial& gj = F[j].GetHeadMonomial();
                        Monomial lcm = Monomial::lcm(gi, gj);
                        Polynomial S = lcm/gi * F[i].GetHeadMonomial() - lcm/gj * F[j].GetHeadMonomial();
                        if (S.ReduceBy(F)) {
                            ready = false;
                            F.push_back(S);
                        }
                    }
                }
            }
        }
    }
}