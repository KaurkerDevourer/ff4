#pragma once

#include "../../util/polynomial.hpp"

namespace NAlgo {
    namespace NUtil {
        using namespace NUtils;

        std::queue<std::pair<size_t, size_t>> GetPairsToCheck(size_t sz) {
            std::queue<std::pair<size_t, size_t>> pairs_to_check;
            for (size_t i = 0; i < sz; i++) {
                for (size_t j = i + 1; j < sz; j++) {
                    pairs_to_check.push({i, j});
                }
            }
            return pairs_to_check;
        }

        template <typename TCoef, typename TComp>
        bool ReduceToZero(Polynomial<TCoef, TComp>& F, const TPolynomials<TCoef, TComp>& polynomialsSet) {
            if (F.IsZero()) {
                return true;
            }
            bool changed = true;
            while(changed && !F.IsZero()) {
                changed = false;
                for (const auto& f : polynomialsSet) {
                    while (!F.IsZero() && F.GetLeadingTerm().IsDivisibleBy(f.GetLeadingTerm())) {
                        F -= f * (F.GetLeadingMonomial() / f.GetLeadingMonomial());
                        changed = true;
                    }
                }
            }
            return F.IsZero();
        }

        template <typename TCoef, typename TComp>
        bool CheckProductCriteria(const Polynomial<TCoef, TComp>& a, const Polynomial<TCoef, TComp>& b) {
            const Monomial<TCoef>& am = a.GetLeadingMonomial();
            const Monomial<TCoef>& bm = b.GetLeadingMonomial();
            const TTerm t = gcd(am.GetTerm(), bm.GetTerm());
            return t.GetDegree() == 0;
        }

        template <typename TCoef, typename TComp>
        bool CheckBasisIsGroebner(const TPolynomials<TCoef, TComp>& basis) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheck(basis.size());
            while(!pairs_to_check.empty()) {
                const Polynomial<TCoef, TComp>& fi = basis[pairs_to_check.front().first];
                const Polynomial<TCoef, TComp>& fj = basis[pairs_to_check.front().second];
                pairs_to_check.pop();
                const Monomial<TCoef>& gi = fi.GetLeadingMonomial();
                const Monomial<TCoef>& gj = fj.GetLeadingMonomial();
                Monomial<TCoef> glcm = Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                Polynomial<TCoef, TComp> S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!ReduceToZero(S, basis)) {
                    return false;
                }
            }
            return true;
        }
    }
}
