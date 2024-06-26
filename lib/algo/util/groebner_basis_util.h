#pragma once

#include "../../util/polynomial.h"
#include "../../util/critical_pair.h"
#include <set>
#include <unordered_set>

namespace FF4 {
    namespace NAlgo {
        namespace NUtil {
            template <typename TCoef, typename TComp>
            using TPairsSet = std::set<NUtils::CriticalPair<TCoef, TComp>, TComp>;

            template <typename TCoef, typename TComp>
            using TPolynomialSet = std::set<NUtils::Polynomial<TCoef, TComp>, TComp>;

            std::queue<std::pair<size_t, size_t>> GetPairsToCheck(size_t sz) {
                std::queue<std::pair<size_t, size_t>> pairs_to_check;
                for (size_t i = 0; i < sz; i++) {
                    for (size_t j = i + 1; j < sz; j++) {
                        pairs_to_check.push({i, j});
                    }
                }
                return pairs_to_check;
            }

            template <typename TCoef, typename TComp, typename TContainter>
            bool InplaceReduceToZero(NUtils::Polynomial<TCoef, TComp>& F, const TContainter& polynomialsSet) {
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
            void EraseByLcm(TPairsSet<TCoef, TComp>& pairs_to_check, const NUtils::Polynomial<TCoef, TComp>& f) {
                for (auto it = pairs_to_check.begin(); it != pairs_to_check.end();) {
                    if (gcd(f.GetLeadingTerm(), it->GetRightTerm()).IsOne()) {
                        ++it;
                        continue;
                    }

                    bool deleted = false;
                    const NUtils::Term& left = it->GetGlcm();
                    for (auto jt = pairs_to_check.begin(); jt != pairs_to_check.end(); ++jt) {
                        if (it == jt) {
                            continue;
                        }
                        if (left.IsDivisibleBy(jt->GetGlcm())) {
                            it = pairs_to_check.erase(it);
                            deleted = true;
                            break;
                        }
                    }

                    if (!deleted) {
                        ++it;
                    }
                }
            }

            template <typename TCoef, typename TComp>
            bool VGBL(const NUtils::CriticalPair<TCoef, TComp>& cp, const NUtils::Term term) {
                if (!term.IsDivisibleBy(gcd(cp.GetLeftTerm(), cp.GetRightTerm()))) {
                    return false;
                }
                std::vector<uint16_t> p1 = cp.GetLeftTerm().GetData();
                std::vector<uint16_t> p2 = term.GetData();
                std::vector<uint16_t> p3 = cp.GetRightTerm().GetData();

                if (p2.size() > p1.size() && p2.size() > p3.size()) {
                    return false;
                }

                size_t sz = std::min(p2.size(), std::min(p1.size(), p3.size()));
                for (size_t i = 0; i < sz; i++) {
                    if (std::min(p1[i], p3[i]) == 0) {
                        continue;
                    }
                    if (std::max(p1[i], p3[i]) >= p2[i]) {
                        continue;
                    }
                    return false;
                }
                if (p2.size() > p1.size()) {
                    for (size_t i = sz; i < std::min(p2.size(), p3.size()); i++) {
                        if (p3[i] >= p2[i]) {
                            continue;
                        }
                        return false;
                    }
                }
                if (p2.size() > p3.size()) {
                    for (size_t i = sz; i < std::min(p2.size(), p1.size()); i++) {
                        if (p1[i] >= p2[i]) {
                            continue;
                        }
                        return false;
                    }
                }
                return true;
            }

            template <typename TCoef, typename TComp>
            void InsertByLcm(TPairsSet<TCoef, TComp>& old_crit_pairs, TPairsSet<TCoef, TComp>& new_crit_pairs, const NUtils::Polynomial<TCoef, TComp>& f) {
                for (const auto& cp : old_crit_pairs) {
                    if (cp.GetGlcm().IsDivisibleBy(f.GetLeadingTerm())) {
                        continue;
                    }
                    if (VGBL(cp, f.GetLeadingTerm())) {
                        continue;
                    }
                    new_crit_pairs.insert(cp);
                }
            }

            template <typename TCoef, typename TComp>
            void InsertByGcd(TPairsSet<TCoef, TComp>& all_crit, TPairsSet<TCoef, TComp>& new_crit_pairs, const NUtils::Polynomial<TCoef, TComp>& f) {
                for (const auto& cp : all_crit) {
                    if (!gcd(f.GetLeadingTerm(), cp.GetRightTerm()).IsOne()) {
                        new_crit_pairs.insert(cp);
                    }
                }
            }

            template <typename TCoef, typename TComp>
            void UpdateCriticalPairs(TPolynomialSet<TCoef, TComp>& polynomials, TPairsSet<TCoef, TComp>& old_crit_pairs, NUtils::Polynomial<TCoef, TComp>& g) {
                g.Normalize();
                TPairsSet<TCoef, TComp> all_crit, new_crit_pairs;
                auto [fit, _] = polynomials.insert(g);
                for (auto it = polynomials.begin(); it != polynomials.end(); ++it) {
                    if (it != fit) {
                        all_crit.insert(NUtils::CriticalPair(*fit, *it));
                    }
                }

                EraseByLcm(all_crit, *fit);
                InsertByGcd(all_crit, new_crit_pairs, *fit);
                InsertByLcm(old_crit_pairs, new_crit_pairs, *fit);
                old_crit_pairs = std::move(new_crit_pairs);
            }

            template <typename TCoef, typename TComp>
            void UpdateBasis(const NUtil::TPolynomialSet<TCoef, TComp>& polynomials, NUtils::TPolynomials<TCoef, TComp>& F) {
                F.clear();
                F.reserve(polynomials.size());
                for (const auto& x : polynomials) {
                    F.push_back(x);
                }
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
            bool CheckProductCriteria(const NUtils::Polynomial<TCoef, TComp>& a, const NUtils::Polynomial<TCoef, TComp>& b) {
                const NUtils::Monomial<TCoef>& am = a.GetLeadingMonomial();
                const NUtils::Monomial<TCoef>& bm = b.GetLeadingMonomial();
                const NUtils::Term t = gcd(am.GetTerm(), bm.GetTerm());
                return t.IsOne();
            }

            template <typename TCoef, typename TComp>
            std::queue<std::pair<size_t, size_t> > GetPairsToCheckWithCriteria(const NUtils::TPolynomials<TCoef, TComp>& basis) {
                std::queue<std::pair<size_t, size_t> > pairs_to_check;
                std::vector<std::pair<NUtils::Term, std::pair<size_t, size_t>>> terms;
                for (size_t i = 0; i < basis.size(); i++) {
                    for (size_t j = i + 1; j < basis.size(); j++) {
                        if (CheckProductCriteria(basis[i], basis[j])) {
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
                std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheckWithCriteria(basis);
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
