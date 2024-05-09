#pragma once

#include "../../util/polynomial.h"

namespace FF4 {
    namespace NAlgo {
        namespace NUtil {

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
            bool InplaceReduceToZero(NUtils::Polynomial<TCoef, TComp>& F, const NUtils::TPolynomials<TCoef, TComp>& polynomialsSet) {
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
            bool CheckProductCriteria(const NUtils::Polynomial<TCoef, TComp>& a, const NUtils::Polynomial<TCoef, TComp>& b) {
                const NUtils::Monomial<TCoef>& am = a.GetLeadingMonomial();
                const NUtils::Monomial<TCoef>& bm = b.GetLeadingMonomial();
                const NUtils::Term t = gcd(am.GetTerm(), bm.GetTerm());
                return t.TotalDegree() == 0;
            }

            template <typename TCoef, typename TComp>
            bool CheckBasisIsGroebner(const NUtils::TPolynomials<TCoef, TComp>& basis) {
                std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheck(basis.size());
                while(!pairs_to_check.empty()) {
                    const NUtils::Polynomial<TCoef, TComp>& fi = basis[pairs_to_check.front().first];
                    const NUtils::Polynomial<TCoef, TComp>& fj = basis[pairs_to_check.front().second];
                    pairs_to_check.pop();
                    const NUtils::Monomial<TCoef>& gi = fi.GetLeadingMonomial();
                    const NUtils::Monomial<TCoef>& gj = fj.GetLeadingMonomial();
                    NUtils::Monomial<TCoef> glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                    NUtils::Polynomial<TCoef, TComp> S = fi * (glcm / gi) - fj * (glcm / gj);
                    if (!InplaceReduceToZero(S, basis)) {
                        return false;
                    }
                }
                return true;
            }

            template <typename TCoef, typename TComp>
            std::queue<std::pair<size_t, size_t> > GetPairsToCheckWithCriterias(const NUtils::TPolynomials<TCoef, TComp>& basis) {
                std::queue<std::pair<size_t, size_t> > pairs_to_check;
                std::vector<std::pair<NUtils::Term, std::pair<size_t, size_t>>> terms;
                for (size_t i = 0; i < basis.size(); i++) {
                    for (size_t j = i + 1; j < basis.size(); j++) {
                        if (NUtil::CheckProductCriteria(basis[i], basis[j])) {
                            continue;
                        }
                        bool isLcmCriteria = false;
                        for (size_t k = 0; k < terms.size(); k++) {
                            if (terms[k].first.IsDivisibleBy(basis[j].GetLeadingTerm()) && terms[k].second.first != j && terms[k].second.second != j) {
                                isLcmCriteria = true;
                                break;
                            }
                            if (terms[k].first.IsDivisibleBy(basis[i].GetLeadingTerm()) && terms[k].second.first != i && terms[k].second.second != i) {
                                isLcmCriteria = true;
                                break;
                            }
                        }
                        if (isLcmCriteria) {
                            continue;
                        }
                        pairs_to_check.push({i, j});
                        terms.push_back({lcm(basis[i].GetLeadingTerm(), basis[j].GetLeadingTerm()), {i, j}});
                    }
                }
                return pairs_to_check;
            }


            template <typename TCoef, typename TComp>
            bool CheckBasisIsGroebnerBig(const NUtils::TPolynomials<TCoef, TComp>& basis) {
                std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheckWithCriterias(basis);
                std::cout << pairs_to_check.size() << std::endl;
                while(!pairs_to_check.empty()) {
                    const NUtils::Polynomial<TCoef, TComp>& fi = basis[pairs_to_check.front().first];
                    const NUtils::Polynomial<TCoef, TComp>& fj = basis[pairs_to_check.front().second];
                    pairs_to_check.pop();
                    const NUtils::Monomial<TCoef>& gi = fi.GetLeadingMonomial();
                    const NUtils::Monomial<TCoef>& gj = fj.GetLeadingMonomial();
                    NUtils::Monomial<TCoef> glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                    NUtils::Polynomial<TCoef, TComp> S = fi * (glcm / gi) - fj * (glcm / gj);
                    if (!InplaceReduceToZero(S, basis)) {
                        return false;
                    }
                }
                return true;
            }
        }
    }
}
