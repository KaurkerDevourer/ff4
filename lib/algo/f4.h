#pragma once
#include "../util/critical_pair.hpp"
#include "util/groebner_basis_util.h"
#include "util/matrix_reduction.h"
#include <cassert>

namespace NAlgo {
    namespace F4 {
        using namespace NUtils;

        template <typename TCoef, typename TComp>
        using TPairsSet = std::multiset<CriticalPair<TCoef, TComp>, TComp>;

        template <typename TCoef, typename TComp>
        using TPairsVector = std::vector<CriticalPair<TCoef, TComp>>;

        // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
        template <typename TCoef, typename TComp>
        TPairsSet<TCoef, TComp> GetPairsToCheckWithCriterias(const TPolynomials<TCoef, TComp>& polynomials) {
            TPairsSet<TCoef, TComp> pairs_to_check;
            for (size_t i = 0; i < polynomials.size(); i++) {
                for (size_t j = i + 1; j < polynomials.size(); j++) {
                    if (NUtil::CheckProductCriteria(polynomials[i], polynomials[j])) {
                        continue;
                    }
                    pairs_to_check.insert(CriticalPair(&polynomials, i, j));
                }
            }
            return pairs_to_check;
        }

        template <typename TCoef, typename TComp>
        TPairsVector<TCoef, TComp> Select(TPairsSet<TCoef, TComp>& pairs_to_check) {
            TPairsVector<TCoef, TComp> selectionGroup;
            uint64_t value = pairs_to_check.begin()->GetDegree();
            while(pairs_to_check.size() && pairs_to_check.begin()->GetDegree() == value) {
                selectionGroup.push_back(*pairs_to_check.begin());
                pairs_to_check.erase(pairs_to_check.begin());
            }
            return selectionGroup;
        }

        template <typename TCoef, typename TComp>
        void UpdateL(TPolynomials<TCoef, TComp>& L, const TTerm& term, const TPolynomials<TCoef, TComp>& F, NUtil::TDiffSet& diff, NUtil::TDiffSet& done) {
            for (const auto& polynomial : F) {
                const auto& t = polynomial.GetHeadMonomial().GetTerm();
                if (term.IsDivisibleBy(t)) {
                    Polynomial<TCoef, TComp> reducer = (term / t) * polynomial;
                    for (const auto& m : reducer.GetMonomials()) {
                        if (!done.contains(m.GetTerm())) {
                            diff.insert(m.GetTerm());
                        }
                    }
                    L.push_back(std::move(reducer));
                    break;
                }
            }
        }

        template <typename TCoef, typename TComp>
        NUtil::TSymbolicPreprocessingResult<TCoef, TComp> SymbolicPreprocessing(TPairsVector<TCoef, TComp>& selected, const TPolynomials<TCoef, TComp>& F) {
            TPolynomials<TCoef, TComp> L;
            L.reserve(selected.size() * 2);
            for (const auto& pair : selected) {
                L.push_back(pair.GetGlcmTerm() / pair.GetLeftTerm() * pair.GetLeft()); // add left
                L.push_back(pair.GetGlcmTerm() / pair.GetRightTerm() * pair.GetRight()); // add right
            }

            NUtil::TDiffSet diff;
            for (const auto& l : L) {
                auto it = diff.begin();
                for (const auto& m : l.GetMonomials()) {
                    it = diff.insert(it, m.GetTerm());
                }
            }

            NUtil::TDiffSet done;
            for (const auto& l : L) {
                diff.erase(l.GetHeadMonomial().GetTerm());
                done.insert(l.GetHeadMonomial().GetTerm());
            }

            while(!diff.empty()) {
                TTerm term = *diff.begin();
                diff.erase(diff.begin());
                done.insert(done.begin(), term);
                UpdateL(L, term, F, diff, done);
            }

            return {L, done};
        }

        template <typename TCoef, typename TComp>
        TPolynomials<TCoef, TComp> Reduce(TPairsVector<TCoef, TComp>& selected, TPolynomials<TCoef, TComp>& F) {
            NUtil::TSymbolicPreprocessingResult<TCoef, TComp> L = SymbolicPreprocessing(selected, F);
            return NUtil::MatrixReduction(L);
        }

        template <typename TCoef, typename TComp>
        void FindGroebnerBasis(TPolynomials<TCoef, TComp>& F) {
            TPairsSet<TCoef, TComp> pairs_to_check = GetPairsToCheckWithCriterias(F);

            while(!pairs_to_check.empty()) {
                TPairsVector<TCoef, TComp> selection_group = Select(pairs_to_check);
                TPolynomials<TCoef, TComp> G = Reduce(selection_group, F);
                for (size_t i = 0; i < G.size(); i++) {
                    size_t idx = F.size();
                    F.push_back(G[i]);
                    for (size_t j = 0; j < idx; j++) {
                        if (NUtil::CheckProductCriteria(F[j], F[idx])) {
                            continue;
                        }
                        pairs_to_check.insert(CriticalPair(&F, j, idx));
                    }
                }
            };
        }
    }
}
