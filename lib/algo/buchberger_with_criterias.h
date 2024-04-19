#pragma once
#include "../util/polynomial.h"
#include <algorithm>
#include <cassert>

namespace NAlgo {
    namespace BuchbergerWithCreterias {
        using namespace NUtils;

        template <typename TCoef>
        bool CheckProductCreteria(Polynomial<TCoef>& a, Polynomial<TCoef>& b) {
                const Monomial<TCoef>& am = a.GetHeadMonomial();
                const Monomial<TCoef>& bm = b.GetHeadMonomial();
                const TTerm t = lcm(am.GetTerm(), bm.GetTerm());
                return (t == am.GetTerm() * bm.GetTerm());
        }

        // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
        template <typename TCoef>
        std::queue<std::pair<size_t, size_t>> GetPairsToCheckWithCriterias(TPolynomials<TCoef> polynomials) {
            std::queue<std::pair<size_t, size_t>> pairs_to_check;
            for (size_t i = 0; i < polynomials.size(); i++) {
                for (size_t j = i + 1; j < polynomials.size(); j++) {
                    if (CheckProductCreteria(polynomials[i], polynomials[j])) {
                        continue;
                    }
                    pairs_to_check.push({i, j});
                }
            }
            return pairs_to_check;
        }

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
        void FindGroebnerBasis(TPolynomials<TCoef>& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheckWithCriterias(F);

            while(!pairs_to_check.empty()) {
                const Polynomial<TCoef>& fi = F[pairs_to_check.front().first];
                const Polynomial<TCoef>& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                const Monomial<TCoef>& gi = fi.GetHeadMonomial();
                const Monomial<TCoef>& gj = fj.GetHeadMonomial();
                Monomial<TCoef> glcm = Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                Polynomial<TCoef> S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!ReduceToZero(S, F)) {
                    size_t idx = F.size();
                    // std::cout << S << std::endl;
                    // std::cout << pairs_to_check.size() << ' ' << F.size() << std::endl;
                    F.push_back(std::move(S));
                    for (size_t i = 0; i < idx; i++) {
                        if (CheckProductCreteria(F[i], F[idx])) {
                            continue;
                        }
                        pairs_to_check.push({i, idx});
                    }
                }
            }
        }
    }
}
