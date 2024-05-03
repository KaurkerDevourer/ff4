#pragma once

#include "../../util/polynomial.hpp"

namespace NAlgo {
    namespace NUtil {
        using namespace NUtils;
        template <typename TCoef>
        bool ReduceToZero(Polynomial<TCoef>& F, TPolynomials<TCoef>& polynomialsSet) {
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
    }
}
