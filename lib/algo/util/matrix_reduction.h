#pragma once
#include "../../util/polynomial.hpp"
#include "../../util/comp.hpp"
#include <set>
#include <map>
#include <cstring>

namespace FF4 {
    namespace NAlgo {
        namespace NUtil {
            using namespace NUtils;
            template <typename TComp>
            using TDiffSet = std::set<TTerm, TComp>;

            template<typename TCoef, typename TComp>
            using TSymbolicPreprocessingResult = std::pair<TPolynomials<TCoef, TComp>, TDiffSet<TComp>>;

            template <typename TCoef, typename TComp>
            TPolynomials<TCoef, TComp> MatrixReduction(NUtil::TSymbolicPreprocessingResult<TCoef, TComp>& L) {
                TDiffSet<TComp>& diffSet = L.second;
                TPolynomials<TCoef, TComp>& F = L.first;
                std::map<TTerm, size_t> Mp;
                std::vector<TTerm> vTerms(diffSet.size());
                size_t cnt = 0;
                auto it = Mp.begin();
                for (const auto& term : diffSet) {
                    it = Mp.insert(it, {term, cnt++});
                    vTerms[vTerms.size() - cnt] = term;
                }
                TCoef matrix[F.size()][diffSet.size()];
                for (size_t i = 0; i < F.size(); i++) {
                    for (size_t j = 0; j < diffSet.size(); j++) {
                        matrix[i][j] = 0;
                    }
                }
                std::set<TTerm> leadingTerms;
                for (size_t i = 0; i < F.size(); i++) {
                    leadingTerms.insert(F[i].GetLeadingTerm());
                    for (const auto& m : F[i].GetMonomials()) {
                        const auto& term = m.GetTerm();
                        matrix[i][Mp[term]] = m.GetCoef();
                    }
                }

                bool used[F.size()];
                for (size_t i = 0; i < F.size(); i++) {
                    used[i] = false;
                }

                for (size_t j = 0; j < diffSet.size(); j++) {
                    for (size_t i = 0; i < F.size(); i++) {
                        if (used[i] || matrix[i][j] == 0) {
                            continue;
                        }
                        used[i] = true;
                        TCoef factor = matrix[i][j];
                        if (factor != 1) {
                            for (size_t k = j; k < diffSet.size(); k++) {
                                matrix[i][k] /= factor;
                            }
                        }
                        
                        for (size_t k = 0; k < F.size(); k++) {
                            if (k == i) {
                                continue;
                            }
                            TCoef factor = matrix[k][j];
                            if (factor != 0) {
                                for (size_t q = j; q < diffSet.size(); q++) {
                                    matrix[k][q] -= matrix[i][q] * factor;
                                }
                            }
                        }
                        break;
                    }
                }
                TPolynomials<TCoef, TComp> reduced;
                reduced.reserve(F.size());
                for (size_t i = 0; i < F.size(); i++) {
                    std::vector<Monomial<TCoef>> mons;
                    bool shouldAdd = false;
                    for (int j = 0; j < diffSet.size(); j++) {
                        if (matrix[i][j] == 0) {
                            continue;
                        }
                        if (mons.empty()) {
                            if (leadingTerms.contains(vTerms[j])) {
                                break;
                            } else {
                                shouldAdd = true;
                            }
                        }
                        mons.emplace_back(vTerms[j], matrix[i][j]);
                    }
                    if (shouldAdd) {
                        reduced.emplace_back(std::move(mons));
                    }
                }
                return std::move(reduced);
            }
        }
    }
}
