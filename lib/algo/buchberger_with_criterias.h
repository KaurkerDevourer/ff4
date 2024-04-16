#pragma once
#include "../util/polynomial.h"
#include <algorithm>
#include <cassert>

namespace NAlgo {
    namespace BuchbergerWithCreterias {

        template <typename TCoef>
        bool CheckProductCreteria(NUtils::Polynomial<TCoef>& a, NUtils::Polynomial<TCoef>& b) {
                const NUtils::Monomial<TCoef>& am = a.GetHeadMonomial();
                const NUtils::Monomial<TCoef>& bm = b.GetHeadMonomial();
                const NUtils::TTerm t = lcm(am.GetTerm(), bm.GetTerm());
                return (t == am.GetTerm() * bm.GetTerm());
        }

        template <typename TCoef>
        bool CheckChainCreteria(NUtils::TPolynomials<TCoef> polynomials, size_t a, size_t b, std::vector<std::vector<size_t> >& pairs_on_each) {
            std::vector<size_t>& pairs_a = pairs_on_each[a];
            std::vector<size_t>& pairs_b = pairs_on_each[b];
            assert(std::is_sorted(pairs_a.begin(), pairs_a.end()));
            assert(std::is_sorted(pairs_b.begin(), pairs_b.end()));
            size_t i = 0, j = 0;

            const NUtils::Monomial<TCoef>& am = polynomials[a].GetHeadMonomial();
            const NUtils::Monomial<TCoef>& bm = polynomials[b].GetHeadMonomial();

            const NUtils::TTerm lcmterm = lcm(am.GetTerm(), bm.GetTerm());
            while(i < pairs_a.size() && j < pairs_b.size()) {
                if (pairs_a[i] == pairs_b[j]) {
                    const NUtils::TTerm k = polynomials[pairs_a[i]].GetHeadMonomial().GetTerm();
                    if (lcmterm.IsDivisibleBy(k)) {
                        return true;
                    }
                    i++;
                    j++;
                } else if (pairs_a[i] < pairs_b[j]) {
                    i++;
                } else {
                    j++;
                }
            }
            return false;
        }
        // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
        template <typename TCoef>
        std::queue<std::pair<size_t, size_t>> GetPairsToCheckWithCriterias(NUtils::TPolynomials<TCoef> polynomials/*, std::vector<std::vector<size_t> >& pairs_on_each*/) {
            std::queue<std::pair<size_t, size_t>> pairs_to_check;
            for (size_t i = 0; i < polynomials.size(); i++) {
                for (size_t j = i + 1; j < polynomials.size(); j++) {
                    if (CheckProductCreteria(polynomials[i], polynomials[j])) {
                        continue;
                    }
                    // if (CheckChainCreteria(polynomials, i, j, pairs_on_each)) {
                    //     continue;
                    // }
                    pairs_to_check.push({i, j});
                    // pairs_on_each[i].push_back(j);
                    // pairs_on_each[j].push_back(i);
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
            // std::vector<std::vector<size_t> > pairs_on_each(F.size());
            std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheckWithCriterias(F/*, pairs_on_each*/);

            while(!pairs_to_check.empty()) {
                const NUtils::Polynomial<TCoef>& fi = F[pairs_to_check.front().first];
                const NUtils::Polynomial<TCoef>& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                const NUtils::Monomial<TCoef>& gi = fi.GetHeadMonomial();
                const NUtils::Monomial<TCoef>& gj = fj.GetHeadMonomial();
                NUtils::Monomial<TCoef> glcm = NUtils::Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                NUtils::Polynomial<TCoef> S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!ReduceToZero(S, F)) {
                    size_t idx = F.size();
                    // std::cout << S << std::endl;
                    // std::cout << pairs_to_check.size() << ' ' << F.size() << std::endl;
                    F.push_back(std::move(S));
                    // pairs_on_each.push_back({});
                    for (size_t i = 0; i < idx; i++) {
                        if (CheckProductCreteria(F[i], F[idx])) {
                            continue;
                        }
                        // if (CheckChainCreteria(F, i, idx, pairs_on_each)) {
                        //     continue;
                        // }
                        pairs_to_check.push({i, idx});
                        // pairs_on_each[i].push_back(idx);
                        // pairs_on_each[idx].push_back(i);
                    }
                }
            }
        }
    }
}
