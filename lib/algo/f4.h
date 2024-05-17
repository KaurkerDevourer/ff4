#pragma once
#include "../util/critical_pair.h"
#include "util/groebner_basis_util.h"
#include "util/matrix_reduction.h"
#include <cassert>

namespace FF4 {
    namespace NAlgo {
        namespace F4 {

            template <typename TCoef, typename TComp>
            using TPairsVector = std::vector<NUtils::CriticalPair<TCoef, TComp>>;

            template <typename TCoef, typename TComp>
            TPairsVector<TCoef, TComp> Select(NUtil::TPairsSet<TCoef, TComp>& pairs_to_check) {
                TPairsVector<TCoef, TComp> selectionGroup;
                NUtils::Term::Degree value = pairs_to_check.begin()->TotalDegree();
                while(pairs_to_check.size() && pairs_to_check.begin()->TotalDegree() == value) {
                    selectionGroup.push_back(*pairs_to_check.begin());
                    pairs_to_check.erase(pairs_to_check.begin());
                }
                return selectionGroup;
            }

            template <typename TCoef, typename TComp>
            void UpdateL(NUtils::TPolynomials<TCoef, TComp>& L, const NUtils::Term& term, const NUtil::TPolynomialSet<TCoef, TComp>& polynomials, NUtil::TTermHashSet& diff, NUtil::TTermHashSet& done) {
                for (const auto& polynomial : polynomials) {
                    const auto& t = polynomial.GetLeadingTerm();
                    if (term.IsDivisibleBy(t)) {
                        NUtils::Polynomial<TCoef, TComp> reducer = (term / t) * polynomial;
                        L.push_back(std::move(reducer));
                        for (const auto& m : L.back().GetMonomials()) {
                            if (!done.contains(m.GetTerm())) {
                                diff.insert(m.GetTerm());
                            }
                        }
                        break;
                    }
                }
            }

            template <typename TCoef, typename TComp>
            NUtil::TSymbolicPreprocessingResult<TCoef, TComp> SymbolicPreprocessing(TPairsVector<TCoef, TComp>& selected, const NUtil::TPolynomialSet<TCoef, TComp>& polynomials) {
                NUtils::TPolynomials<TCoef, TComp> L;
                L.reserve(selected.size() * 3);
                for (const auto& pair : selected) {
                    L.push_back(pair.GetGlcm() / pair.GetLeftTerm() * pair.GetLeft());
                    L.push_back(pair.GetGlcm() / pair.GetRightTerm() * pair.GetRight());
                }

                NUtil::TTermHashSet diff;
                for (const auto& l : L) {
                    auto it = diff.begin();
                    for (const auto& m : l.GetMonomials()) {
                        it = diff.insert(it, m.GetTerm());
                    }
                }

                NUtil::TTermHashSet done;
                for (const auto& l : L) {
                    diff.erase(l.GetLeadingTerm());
                    done.insert(l.GetLeadingTerm());
                }

                while(!diff.empty()) {
                    const NUtils::Term& term = *diff.begin();
                    auto extracted = diff.extract(diff.begin());
                    done.insert(std::move(extracted));
                    UpdateL(L, term, polynomials, diff, done);
                }
                std::vector<NUtils::Term> done_sorted;
                done_sorted.reserve(done.size());
                for (auto& x : done) {
                    done_sorted.push_back(x);
                }
                std::sort(done_sorted.begin(), done_sorted.end(), TComp());

                return {std::move(L), std::move(done_sorted)};
            }

            template <typename TCoef, typename TComp>
            NUtils::TPolynomials<TCoef, TComp> Reduce(TPairsVector<TCoef, TComp>& selected, NUtil::TPolynomialSet<TCoef, TComp>& polynomials) {
                NUtil::TSymbolicPreprocessingResult<TCoef, TComp> L = SymbolicPreprocessing(selected, polynomials);
                return NUtil::MatrixReduction(L);
            }

            template <typename TCoef, typename TComp>
            void FindGroebnerBasis(NUtils::TPolynomials<TCoef, TComp>& F) {
                NUtil::TPolynomialSet<TCoef, TComp> polynomials;
                NUtil::TPairsSet<TCoef, TComp> pairs_to_check;
                for (auto& f : F) {
                    NUtil::UpdateCriticalPairs(polynomials, pairs_to_check, f);
                }

                while(!pairs_to_check.empty()) {
                    TPairsVector<TCoef, TComp> selection_group = Select(pairs_to_check);
                    NUtils::TPolynomials<TCoef, TComp> G = Reduce(selection_group, polynomials);
                    for (auto& g : G) {
                        NUtil::UpdateCriticalPairs(polynomials, pairs_to_check, g);
                    }
                }
                NUtil::UpdateBasis(polynomials, F);
            }
        }
    }
}
