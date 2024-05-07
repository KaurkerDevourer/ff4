#pragma once
#include "../util/critical_pair.h"
#include "util/groebner_basis_util.h"
#include "util/matrix_reduction.h"
#include <cassert>

namespace FF4 {
    namespace NAlgo {
        namespace F4 {

            template <typename TCoef, typename TComp>
            using TPairsSet = std::set<NUtils::CriticalPair<TCoef, TComp>, TComp>;

            template <typename TCoef, typename TComp>
            using TPairsVector = std::vector<NUtils::CriticalPair<TCoef, TComp>>;

            // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
            template <typename TCoef, typename TComp>
            TPairsSet<TCoef, TComp> GetPairsToCheckWithCriteria(const NUtils::TPolynomials<TCoef, TComp>& polynomials) {
                TPairsSet<TCoef, TComp> pairs_to_check;
                for (size_t i = 0; i < polynomials.size(); i++) {
                    for (size_t j = i + 1; j < polynomials.size(); j++) {
                        if (NUtil::CheckProductCriteria(polynomials[i], polynomials[j])) {
                            continue;
                        }
                        pairs_to_check.insert(NUtils::CriticalPair(polynomials, i, j));
                    }
                }
                return pairs_to_check;
            }

            template <typename TCoef, typename TComp>
            TPairsVector<TCoef, TComp> Select(TPairsSet<TCoef, TComp>& pairs_to_check) {
                TPairsVector<TCoef, TComp> selectionGroup;
                uint64_t value = pairs_to_check.begin()->TotalDegree();
                while(pairs_to_check.size() && pairs_to_check.begin()->TotalDegree() == value) {
                    selectionGroup.push_back(*pairs_to_check.begin());
                    pairs_to_check.erase(pairs_to_check.begin());
                }
                return selectionGroup;
            }

            template <typename TCoef, typename TComp>
            void UpdateL(NUtils::TPolynomials<TCoef, TComp>& L, const NUtils::Term& term, const NUtils::TPolynomials<TCoef, TComp>& F, NUtil::TTermSet<TComp>& diff, NUtil::TTermSet<TComp>& done) {
                for (const auto& polynomial : F) {
                    const auto& t = polynomial.GetLeadingTerm();
                    if (term.IsDivisibleBy(t)) {
                        NUtils::Polynomial<TCoef, TComp> reducer = (term / t) * polynomial;
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
            NUtil::TSymbolicPreprocessingResult<TCoef, TComp> SymbolicPreprocessing(TPairsVector<TCoef, TComp>& selected, const NUtils::TPolynomials<TCoef, TComp>& F) {
                NUtils::TPolynomials<TCoef, TComp> L;
                L.reserve(selected.size() * 2);
                for (const auto& pair : selected) {
                    L.push_back(pair.GetGlcmTerm() / pair.GetLeftTerm() * pair.GetLeft());
                    L.push_back(pair.GetGlcmTerm() / pair.GetRightTerm() * pair.GetRight());
                }

                NUtil::TTermSet<TComp> diff;
                for (const auto& l : L) {
                    auto it = diff.begin();
                    for (const auto& m : l.GetMonomials()) {
                        it = diff.insert(it, m.GetTerm());
                    }
                }

                NUtil::TTermSet<TComp> done;
                for (const auto& l : L) {
                    diff.erase(l.GetLeadingTerm());
                    done.insert(l.GetLeadingTerm());
                }

                while(!diff.empty()) {
                    NUtils::Term term = *diff.begin();
                    diff.erase(diff.begin());
                    done.insert(done.begin(), term);
                    UpdateL(L, term, F, diff, done);
                }

                return {L, done};
            }

            template <typename TCoef, typename TComp>
            NUtils::TPolynomials<TCoef, TComp> Reduce(TPairsVector<TCoef, TComp>& selected, NUtils::TPolynomials<TCoef, TComp>& F) {
                NUtil::TSymbolicPreprocessingResult<TCoef, TComp> L = SymbolicPreprocessing(selected, F);
                return NUtil::MatrixReduction(L);
            }

            template <typename TCoef, typename TComp>
            void FindGroebnerBasis(NUtils::TPolynomials<TCoef, TComp>& F) {
                for (auto& f : F) {
                    f.Normalize();
                }
                TPairsSet<TCoef, TComp> pairs_to_check = GetPairsToCheckWithCriteria(F);

                while(!pairs_to_check.empty()) {
                    TPairsVector<TCoef, TComp> selection_group = Select(pairs_to_check);
                    NUtils::TPolynomials<TCoef, TComp> G = Reduce(selection_group, F);
                    for (size_t i = 0; i < G.size(); i++) {
                        G[i].Normalize();
                        const NUtils::Polynomial<TCoef, TComp>& h = G[i];
                        size_t idx = F.size();
                        F.push_back(h);
                        for (size_t j = 0; j < idx; j++) {
                            // KIND OF GMI INSTALLATION
                            if (NUtil::CheckProductCriteria(F[j], h)) {
                                continue;
                            }
                            for (auto it = pairs_to_check.begin(); it != pairs_to_check.end();) {
                                const NUtils::CriticalPair<TCoef, TComp>& cp = (*it);
                                if (cp.GetGlcmTerm().IsDivisibleBy(h.GetLeadingTerm()) &&
                                lcm(cp.GetLeftTerm(), h.GetLeadingTerm()) != cp.GetGlcmTerm() &&
                                lcm(cp.GetRightTerm(), h.GetLeadingTerm()) != cp.GetGlcmTerm()) {
                                    it = pairs_to_check.erase(it);
                                } else {
                                    ++it;
                                }
                            }
                            pairs_to_check.insert(NUtils::CriticalPair(F, j, idx));
                        }
                    }
                };
            }
        }
    }
}
