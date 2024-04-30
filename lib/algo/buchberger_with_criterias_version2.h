#pragma once
#include "util/groebner_basis_util.h"
#include <algorithm>
#include <cassert>

namespace NAlgo {
    namespace BuchbergerWithCreteriasVersion2 {
        using namespace NUtils;
        using TPairsQueue = std::queue<std::pair<size_t, size_t>>;

        template <typename TCoef>
        bool CheckProductCreteria(const Polynomial<TCoef>& a,const  Polynomial<TCoef>& b) {
                const Monomial<TCoef>& am = a.GetHeadMonomial();
                const Monomial<TCoef>& bm = b.GetHeadMonomial();
                const TTerm t = lcm(am.GetTerm(), bm.GetTerm());
                return (t == am.GetTerm() * bm.GetTerm());
        }

        // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
        template <typename TCoef>
        TPairsQueue GetPairsToCheckWithCriterias(const TPolynomials<TCoef>& polynomials) {
            TPairsQueue pairs_to_check;
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

        // THIS IS THE DIFFERENCE BETWEEN no-number version. Problem is - in cyclic and katsura, this version is 10x times faster. But for sym-sym for some reason this version is 10x slower.
        template <typename TCoef>
        void TailReduce(TPolynomials<TCoef>& F, Polynomial<TCoef>& S) {
            const TTerm& t = S.GetHeadMonomial().GetTerm();
            for (size_t i = 0; i < F.size(); i++) {
                if (!F[i].IsRemoved() && F[i].GetHeadMonomial().GetTerm().IsDivisibleBy(t)) {
                    F[i].MakeRemoved();
                }
            }
        }

        template <typename TCoef>
        void FindGroebnerBasis(TPolynomials<TCoef>& F) {
            std::queue<std::pair<size_t, size_t> > pairs_to_check = GetPairsToCheckWithCriterias(F);

            while(!pairs_to_check.empty()) {
                const Polynomial<TCoef>& fi = F[pairs_to_check.front().first];
                const Polynomial<TCoef>& fj = F[pairs_to_check.front().second];
                pairs_to_check.pop();
                if (fi.IsRemoved() || fj.IsRemoved()) {
                    continue;
                }
                const Monomial<TCoef>& gi = fi.GetHeadMonomial();
                const Monomial<TCoef>& gj = fj.GetHeadMonomial();
                Monomial<TCoef> glcm = Monomial(lcm(gi.GetTerm(), gj.GetTerm()), TCoef(1));
                Polynomial<TCoef> S = fi * (glcm/gi) - fj * (glcm/gj);
                if (!NUtil::ReduceToZero(S, F)) {
                    TailReduce(F, S);
                    size_t idx = F.size();
                    // std::cout << S << std::endl;
                    // std::cout << pairs_to_check.size() << ' ' << F.size() << std::endl;
                    F.push_back(std::move(S));
                    for (size_t i = 0; i < idx; i++) {
                        if (F[i].IsRemoved()) {
                            continue;
                        }
                        if (CheckProductCreteria(F[i], F[idx])) {
                            continue;
                        }
                        pairs_to_check.push({i, idx});
                    }
                }
            }
            size_t j = 0;
            for (size_t i = 0; i + j < F.size();) {
                if (F[i + j].IsRemoved()) {
                    j++;
                } else {
                    F[i] = F[i + j];
                    i++;
                }
            }
            F.erase(F.end() - j, F.end());
        }
    }
}
