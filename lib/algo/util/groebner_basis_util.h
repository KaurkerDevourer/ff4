#pragma once

#include "../../util/polynomial.hpp"

namespace NAlgo {
    namespace NUtil {
        using namespace NUtils;
        template <typename TCoef>
        bool ReduceToZero(Polynomial<TCoef>& F, const TPolynomials<TCoef>& polynomialsSet) {
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
        bool CheckProductCriteria(const Polynomial<TCoef>& a, const Polynomial<TCoef>& b) {
            const Monomial<TCoef>& am = a.GetHeadMonomial();
            const Monomial<TCoef>& bm = b.GetHeadMonomial();
            const TTerm t = lcm(am.GetTerm(), bm.GetTerm());
            return (t == am.GetTerm() * bm.GetTerm());
        }

        template <typename TCoef>
        bool CheckBasisIsGroebner(const TPolynomials<TCoef>& basis) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheck(basis.size());
            while(!pairs_to_check.empty()) {
                const Polynomial<TCoef>& fi = basis[pairs_to_check.front().first];
                const Polynomial<TCoef>& fj = basis[pairs_to_check.front().second];
                pairs_to_check.pop();
                const Monomial<TCoef>& gi = fi.GetHeadMonomial();
                const Monomial<TCoef>& gj = fj.GetHeadMonomial();
                Monomial<TCoef> glcm = Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                Polynomial<TCoef> S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!ReduceToZero(S, basis)) {
                    return false;
                }
            }
            return true;
        }
    }
}
