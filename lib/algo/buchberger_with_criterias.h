#pragma once
#include "../util/polynomial.h"

namespace NAlgo {
    namespace BuchbergerWithCreterias {

        template <typename TCoef>
        bool CheckProductCreteria(NUtils::Polynomial<TCoef>& a, NUtils::Polynomial<TCoef>& b) {
                const NUtils::Monomial<TCoef>& am = a.GetHeadMonomial();
                const NUtils::Monomial<TCoef>& bm = b.GetHeadMonomial();
                const NUtils::TTerm t = lcm(am.GetTerm(), bm.GetTerm());
                return (t == am.GetTerm() * bm.GetTerm());
        }

        // template <typename TCoef>
        // bool CheckSecondCreteria(TPolynomial<TCoef>& a, TPolynomial<TCoef>& b) {
        //         const NUtils::Monomial<TCoef>& am = a.GetHeadMonomial();
        //         const NUtils::Monomial<TCoef>& bm = b.GetHeadMonomial();
        //         const NUtils::TTerm t = lcm(am.GetTerm(), bm.GetTerm());
        //         return t != am.GetTerm() && t != bm.GetTerm();
        // }
        // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
        template <typename TCoef>
        std::queue<std::pair<size_t, size_t>> GetPairsToCheckWithCriterias(NUtils::TPolynomials<TCoef> Polynomials) {
            std::queue<std::pair<size_t, size_t>> pairs_to_check;
            for (size_t i = 0; i < Polynomials.size(); i++) {
                for (size_t j = i + 1; j < Polynomials.size(); j++) {
                    if (CheckProductCreteria(Polynomials[i], Polynomials[j])) {
                        continue;
                    }
                    // if (!CheckSecondCreteria(Polynomials[i], Polynomials[j])) {
                    //     continue;
                    // }
                    pairs_to_check.push({i, j});
                }
            }
            return pairs_to_check;
        }

        template <typename TCoef>
        bool ReduceToZero(NUtils::Polynomial<TCoef>& F, NUtils::TPolynomials<TCoef>& polynomialsSet) {
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
        void FindGroebnerBasis(NUtils::TPolynomials<TCoef>& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheckWithCriterias(F);

            while(!pairs_to_check.empty()) {
                const NUtils::Polynomial<TCoef>& fi = F[pairs_to_check.front().first];
                const NUtils::Polynomial<TCoef>& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                const NUtils::Monomial<TCoef>& gi = fi.GetHeadMonomial();
                const NUtils::Monomial<TCoef>& gj = fj.GetHeadMonomial();
                NUtils::Monomial<TCoef> glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                NUtils::Polynomial<TCoef> S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!ReduceToZero(S, F)) {
                    for (size_t i = 0; i < F.size(); i++) {
                        if (!CheckProductCreteria(F[i], S)) {
                            pairs_to_check.push({i, F.size()});
                        }
                    }
                    F.push_back(S);
                }
            }
        }
    }
}
