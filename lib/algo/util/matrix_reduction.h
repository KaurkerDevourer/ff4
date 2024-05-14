#pragma once
#include "../../util/polynomial.h"
#include "../../util/comp.h"
#include "../../util/matrix.h"
#include <set>
#include <map>
#include <cstring>

namespace FF4 {
    namespace NAlgo {
        namespace NUtil {
            template <typename TComp>
            using TTermSet = std::set<NUtils::Term, TComp>;

            template<typename TCoef, typename TComp>
            using TSymbolicPreprocessingResult = std::pair<NUtils::TPolynomials<TCoef, TComp>, TTermSet<TComp>>;

            template <typename TCoef, typename TComp>
            NUtils::Matrix<TCoef> FillMatrixAndLeadingTerms(NUtils::TPolynomials<TCoef, TComp>& F, std::map<NUtils::Term, size_t>& Mp, TTermSet<TComp>& leadingTerms, std::vector<NUtils::Term>& vTerms, const TTermSet<TComp>& diffSet) {
                size_t cnt = 0;
                size_t swp = 0;
                auto it = Mp.begin();
                std::vector<bool> not_pivot(F.size());
                for (size_t i = 0; i < F.size(); i++) {
                    auto [_, inserted] = leadingTerms.insert(F[i].GetLeadingTerm());
                    if (!inserted) {
                        not_pivot[i] = true;
                        swp++;
                        continue;
                    }
                    it = Mp.insert(it, {F[i].GetLeadingTerm(), cnt});
                    vTerms[cnt] = F[i].GetLeadingTerm();
                    cnt++;
                }

                NUtils::Matrix<TCoef> matrix(F.size(), diffSet.size(), F.size() - swp);
                cnt = diffSet.size() - 1;
                for (const auto& term : diffSet) {
                    if (Mp.find(term) == Mp.end()) {
                        it = Mp.insert(it, {term, cnt});
                        vTerms[cnt] = term;
                        cnt--;
                    }
                }

                for (size_t i = 0, j = 0; i < F.size(); i++) {
                    if (not_pivot[i]) {
                        j++;
                        continue;
                    }
                    for (const auto& m : F[i].GetMonomials()) {
                        const auto& term = m.GetTerm();
                        matrix(i - j, Mp[term]) = m.GetCoef();
                    }
                }

                for (size_t i = 0, j = 0; i < F.size(); i++) {
                    if (!not_pivot[i]) {
                        continue;
                    }
                    for (const auto& m : F[i].GetMonomials()) {
                        const auto& term = m.GetTerm();
                        matrix(F.size() - 1 - j, Mp[term]) = m.GetCoef();
                    }
                    j++;
                }

                return matrix;
            }

            template <typename TCoef>
            void GaussElimination(NUtils::Matrix<TCoef>& matrix, size_t pivots) {
                std::vector<bool> used(matrix.N_);
                for (size_t j = pivots; j < matrix.M_; j++) {
                    for (size_t i = pivots; i < matrix.N_; i++) {
                        if (used[i] || matrix(i, j) == 0) {
                            continue;
                        }
                        used[i] = true;
                        TCoef factor = matrix(i, j);
                        if (factor != 1) {
                            for (size_t k = j; k < matrix.M_; k++) {
                                matrix(i, k) /= factor;
                            }
                        }
                        
                        for (size_t k = pivots; k < matrix.N_; k++) {
                            if (k == i) {
                                continue;
                            }
                            TCoef factor = matrix(k, j);
                            if (factor != 0) {
                                if (factor == 1) {
                                    for (size_t q = j; q < matrix.M_; q++) {
                                        matrix(k, q) -= matrix(i, q);
                                    }
                                } else {
                                    for (size_t q = j; q < matrix.M_; q++) {
                                        matrix(k, q) -= matrix(i, q) * factor;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }

            template <typename TCoef, typename TComp>
            NUtils::TPolynomials<TCoef, TComp> GetReducedPolynomials(const NUtils::Matrix<TCoef>& matrix, const std::vector<NUtils::Term>& vTerms, size_t pivots) {
                NUtils::TPolynomials<TCoef, TComp> reduced;
                reduced.reserve(matrix.N_ - pivots);
                for (size_t i = pivots; i < matrix.N_; i++) {
                    std::vector<NUtils::Monomial<TCoef>> mons;
                    for (int j = 0; j < matrix.M_; j++) {
                        if (matrix(i, j) == 0) {
                            continue;
                        }
                        mons.emplace_back(vTerms[j], matrix(i, j));
                    }
                    if (!mons.empty()) {
                        reduced.emplace_back(std::move(mons));
                    }
                }
                return reduced;
            }

            template <typename TCoef>
            void TRSM(NUtils::Matrix<TCoef>& matrix, size_t pivots) {
                for (size_t j = pivots - 1; j > 0; j--) {
                    std::vector<size_t> next;
                    next.reserve((matrix.M_ - pivots) / 2);
                    for (size_t k = pivots; k < matrix.M_; k++) {
                        if (matrix(j, k) != 0) {
                            next.push_back(k);
                        }
                    }
                    for (size_t i = 0; i < j; i++) {
                        if (matrix(i, j) == 0) {
                            continue;
                        }
                        TCoef factor = matrix(i, j);

                        for (size_t k = 0; k < next.size(); k++) {
                            matrix(i, next[k]) -= factor * matrix(j, next[k]);
                        }

                        matrix(i, j) = 0;
                    }
                }
            }

            template <typename TCoef>
            void AXPY(NUtils::Matrix<TCoef>& matrix, size_t pivots) {
                for (size_t i = 0; i < pivots; i++) {
                    std::vector<size_t> next;
                    next.reserve((matrix.M_ - pivots) / 2);
                    for (size_t j = pivots; j < matrix.M_; j++) {
                        if (matrix(i, j) != 0) {
                            next.push_back(j);
                        }
                    }
                    for (size_t k = pivots; k < matrix.N_; k++) {
                        if (matrix(k, i) == 0) {
                            continue;
                        }
                        for (size_t j = 0; j < next.size(); j++) {
                            matrix(k, next[j]) -= matrix(k, i) * matrix(i, next[j]);
                        }
                        matrix(k, i) = 0;
                    }
                }
            }

            template <typename TCoef, typename TComp>
            NUtils::TPolynomials<TCoef, TComp> MatrixReduction(TSymbolicPreprocessingResult<TCoef, TComp>& L) {
                TTermSet<TComp>& diffSet = L.second;
                NUtils::TPolynomials<TCoef, TComp>& F = L.first;
                std::sort(F.begin(), F.end(), [](const NUtils::Polynomial<TCoef, TComp>& a, const NUtils::Polynomial<TCoef, TComp>& b){
                    return TComp()(b, a);
                });

                std::map<NUtils::Term, size_t> Mp;
                std::vector<NUtils::Term> vTerms(diffSet.size());

                TTermSet<TComp> leadingTerms;
                NUtils::Matrix<TCoef> matrix = FillMatrixAndLeadingTerms(F, Mp, leadingTerms, vTerms, diffSet);
                size_t pivots = matrix.pivots_;
                TRSM(matrix, pivots);
                AXPY(matrix, pivots);

                GaussElimination(matrix, pivots);

                return GetReducedPolynomials<TCoef, TComp>(matrix, vTerms, pivots);
            }
        }
    }
}
