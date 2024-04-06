#pragma once
#include "../util/polynomial.h"

namespace NAlgo {
    namespace Buchberger {
        void FindGroebnerBasis(NUtils::TPolynomials& F);

        bool ReduceToZero(NUtils::Polynomial& F, NUtils::TPolynomials& polynomialsSet);
    }
}
