#pragma once

namespace NAlgo {
    namespace F4 {
        template <typename TCoef>
        void FindGroebnerBasis(NUtils::TPolynomials<TCoef>& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = NUtils::GetPairsToCheck(F.size());
        }
    }
}