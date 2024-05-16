#pragma once
#include "../util/critical_pair.h"
#include "util/groebner_basis_util.h"
#include "util/matrix_reduction.h"
#include <unordered_map>
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
            void UpdateL(NUtils::TPolynomials<TCoef, TComp>& L, size_t idx, const NUtil::TPolynomialSet<TCoef, TComp>& polynomials, NUtil::TTermHashSet& termsMap, std::vector<NUtils::Term>& terms) {
                for (const auto& polynomial : polynomials) {
                    const auto& t = polynomial.GetLeadingTerm();
                    if (terms[idx].IsDivisibleBy(t)) {
                        NUtils::Polynomial<TCoef, TComp> reducer = (terms[idx] / t) * polynomial;
                        for (const auto& m : reducer.GetMonomials()) {
                            auto [_, inserted] = termsMap.insert(m.GetTerm());
                            if (inserted) {
                                terms.push_back(m.GetTerm());
                            }
                        }
                        L.push_back(std::move(reducer));
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

                NUtil::TTermHashSet termsMap;
                std::vector<NUtils::Term> terms;
                for (const auto& l : L) {
                    const std::vector<NUtils::Monomial<TCoef>>& monomials = l.GetMonomials();
                    for (size_t i = 0; i < monomials.size(); i++) {
                        auto [_, inserted] = termsMap.insert(monomials[i].GetTerm());
                        if (inserted) {
                            terms.push_back(monomials[i].GetTerm());
                        }
                    }
                }

                size_t i = 0;
                while(i < terms.size()) {
                    UpdateL(L, i, polynomials, termsMap, terms);
                    i++;
                }

                std::sort(terms.begin(), terms.end(), TComp());

                return {L, terms};
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
