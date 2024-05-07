#pragma once
#include "util/groebner_basis_util.h"
#include <algorithm>
#include <cassert>

namespace FF4 {
    namespace NAlgo {
        namespace BuchbergerWithCreteria {
            using TPairsQueue = std::queue<std::pair<size_t, size_t>>;

            // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
            template <typename TCoef, typename TComp>
            TPairsQueue GetPairsToCheckWithCriteria(const NUtils::TPolynomials<TCoef, TComp>& polynomials) {
                TPairsQueue pairs_to_check;
                for (size_t i = 0; i < polynomials.size(); i++) {
                    for (size_t j = i + 1; j < polynomials.size(); j++) {
                        if (NUtil::CheckProductCriteria(polynomials[i], polynomials[j])) {
                            continue;
                        }
                        pairs_to_check.push({i, j});
                    }
                }
                return pairs_to_check;
            }

            template <typename TCoef, typename TComp>
            bool InplaceReduceToZero(NUtils::Polynomial<TCoef, TComp>& F, NUtils::TPolynomials<TCoef, TComp>& polynomialsSet) {
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
            void FindGroebnerBasis(NUtils::TPolynomials<TCoef, TComp>& F) {
                TPairsQueue pairs_to_check = GetPairsToCheckWithCriteria(F);

                while(!pairs_to_check.empty()) {
                    const NUtils::Polynomial<TCoef, TComp>& fi = F[pairs_to_check.front().first];
                    const NUtils::Polynomial<TCoef, TComp>& fj = F[pairs_to_check.front().second];
                    pairs_to_check.pop();
                    const NUtils::Monomial<TCoef>& gi = fi.GetLeadingMonomial();
                    const NUtils::Monomial<TCoef>& gj = fj.GetLeadingMonomial();
                    NUtils::Monomial<TCoef> glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                    NUtils::Polynomial<TCoef, TComp> S = fi * (glcm/gi) - fj * (glcm/gj);
                    if (!NUtil::InplaceReduceToZero(S, F)) {
                        size_t idx = F.size();
                        // std::cout << S << std::endl;
                        // std::cout << pairs_to_check.size() << ' ' << F.size() << std::endl;
                        F.push_back(std::move(S));
                        for (size_t i = 0; i < idx; i++) {
                            if (NUtil::CheckProductCriteria(F[i], F[idx])) {
                                continue;
                            }
                            pairs_to_check.push({i, idx});
                        }
                    }
                }
            }
        }
    }
}
