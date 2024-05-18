#pragma once
#include "../../util/polynomial.h"
#include "../../util/comp.h"
#include "../../util/matrix.h"
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <thread>

namespace FF4 {
    namespace NAlgo {
        namespace NUtil {
            using TTermHashSet = std::unordered_set<NUtils::Term, NUtils::TermHasher>;

            template<typename TCoef, typename TComp>
            using TSymbolicPreprocessingResult = std::pair<NUtils::TPolynomials<TCoef, TComp>, std::vector<NUtils::Term>>;

            template <typename TCoef, typename TComp>
            size_t FillMatrix(NUtils::TPolynomials<TCoef, TComp>& F, NUtils::Matrix<TCoef>& matrix, std::vector<NUtils::Term>& vTerms, const std::vector<NUtils::Term>& diffSet, std::vector<std::vector<size_t> >& nnext) {
                size_t cnt = 0;
                size_t swp = 0;
                std::vector<bool> not_pivot(F.size());
                TTermHashSet leadingTerms;
                std::unordered_map<NUtils::Term, size_t, NUtils::TermHasher> Mp;
                for (size_t i = 0; i < F.size(); i++) {
                    auto [_, inserted] = leadingTerms.insert(F[i].GetLeadingTerm());
                    if (!inserted) {
                        not_pivot[i] = true;
                        swp++;
                        continue;
                    }
                    Mp[F[i].GetLeadingTerm()] = cnt;
                    vTerms[cnt] = F[i].GetLeadingTerm();
                    cnt++;
                }

                cnt = diffSet.size() - 1;
                for (auto& term : diffSet) {
                    if (Mp.find(term) == Mp.end()) {
                        Mp[term] = cnt;
                        vTerms[cnt] = term;
                        cnt--;
                    }
                }

                nnext.reserve(F.size() - swp);
                for (size_t i = 0, j = 0; i < F.size(); i++) {
                    if (not_pivot[i]) {
                        j++;
                        continue;
                    }
                    std::vector<size_t> next;
                    next.reserve(F[i].GetMonomials().size());
                    for (const auto& m : F[i].GetMonomials()) {
                        const auto& term = m.GetTerm();
                        size_t column = Mp[term];
                        matrix(i - j, column) = m.GetCoef();
                        next.push_back(column);
                    }
                    nnext.push_back(std::move(next));
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

                return F.size() - swp;
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
                                for (size_t q = j; q < matrix.M_; q++) {
                                    matrix(k, q) -= matrix(i, q) * factor;
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
            void NOTRSM(NUtils::Matrix<TCoef>& matrix, size_t pivots, const std::vector<std::vector<size_t> >& nnext, size_t numThreads) {
                if (numThreads == 1) {
                    for (size_t i = 0; i < pivots; i++) {
                        const auto& next = nnext[i];
                        for (size_t j = pivots; j < matrix.N_; j++) {
                            if (matrix(j, i) != 0) {
                                TCoef factor = matrix(j, i);
                                for (size_t k = 0; k < next.size(); k++) {
                                    matrix(j, next[k]) -= factor * matrix(i, next[k]);
                                }
                            }
                        }
                    }
                } else {
                    for (size_t i = 0; i < pivots; i++) {
                        const auto& next = nnext[i];
                        std::vector<std::thread> p;
                        for (size_t j = pivots; j < matrix.N_; j++) {
                            if (matrix(j, i) != 0) {
                                p.push_back(std::thread([&]() {
                                    const auto local_next = next;
                                    TCoef factor = matrix(j, i);
                                    for (size_t k = 0; k < local_next.size(); k++) {
                                        matrix(j, local_next[k]) -= factor * matrix(i, local_next[k]);
                                    }
                                }));
                            }
                        }
                        for (size_t i = 0; i < p.size(); i++) {
                            p[i].join();
                        }
                    }
                }
            }

            template <typename TCoef>
            void TRSM(NUtils::Matrix<TCoef>& matrix, size_t pivots) {
                for (size_t j = pivots - 1; j > 0; j--) {
                    for (size_t i = 0; i < j; i++) {
                        if (matrix(i, j) == 0) {
                            continue;
                        }
                        TCoef factor = matrix(i, j);

                        for (size_t k = pivots; k < matrix.M_; k++) {
                            matrix(i, k) -= factor * matrix(j, k);
                        }

                        matrix(i, j) = 0;
                    }
                }
            }

            template <typename TCoef>
            void AXPY(NUtils::Matrix<TCoef>& matrix, size_t pivots) {
                for (size_t i = 0; i < pivots; i++) {
                    for (size_t k = pivots; k < matrix.N_; k++) {
                        if (matrix(k, i) == 0) {
                            continue;
                        }
                        for (size_t j = pivots; j < matrix.M_; j++) {
                            matrix(k, j) -= matrix(k, i) * matrix(i, j);
                        }
                        matrix(k, i) = 0;
                    }
                }
            }

            template <typename TCoef, typename TComp>
            NUtils::TPolynomials<TCoef, TComp> MatrixReduction(TSymbolicPreprocessingResult<TCoef, TComp>& L, size_t numThreads) {
                std::vector<NUtils::Term>& diffSet = L.second;
                NUtils::TPolynomials<TCoef, TComp>& F = L.first;
                std::sort(F.begin(), F.end(), [](const NUtils::Polynomial<TCoef, TComp>& a, const NUtils::Polynomial<TCoef, TComp>& b){
                    return TComp()(b, a);
                });

                std::vector<NUtils::Term> vTerms(diffSet.size());

                NUtils::Matrix<TCoef> matrix(F.size(), diffSet.size());
                std::vector<std::vector<size_t>> nnext;
                size_t pivots = FillMatrix(F, matrix, vTerms, diffSet, nnext);
                NOTRSM(matrix, pivots, nnext, numThreads); // switch comments, to test gbla.
                // TRSM(matrix, pivots);
                // AXPY(matrix, pivots);

                GaussElimination(matrix, pivots);

                return GetReducedPolynomials<TCoef, TComp>(matrix, vTerms, pivots);
            }
        }
    }
}
