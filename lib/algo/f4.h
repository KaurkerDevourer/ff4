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

            template <typename TCoef, typename TComp>
            using TPolynomialSet = std::set<NUtils::Polynomial<TCoef, TComp>, TComp>;

            template <typename TCoef, typename TComp>
            using TEvaluatedResults = std::vector<std::pair<TPolynomialSet<TCoef, TComp>, NUtils::TPolynomials<TCoef, TComp>>>;

            template <typename TCoef, typename TComp>
            TPairsVector<TCoef, TComp> Select(TPairsSet<TCoef, TComp>& pairs_to_check) {
                TPairsVector<TCoef, TComp> selectionGroup;
                NUtils::Term::Degree value = pairs_to_check.begin()->TotalDegree();
                while(pairs_to_check.size() && pairs_to_check.begin()->TotalDegree() == value) {
                    selectionGroup.push_back(*pairs_to_check.begin());
                    pairs_to_check.erase(pairs_to_check.begin());
                }
                return selectionGroup;
            }

            template <typename TCoef, typename TComp>
            NUtils::Polynomial<TCoef, TComp> Simplify(const NUtils::Term& term, const NUtils::Polynomial<TCoef, TComp>& polynomial, const TEvaluatedResults<TCoef, TComp>& evaluatedResults) {
                for (const auto& div : term.GetAllDivisors()) {
                    const NUtils::Polynomial<TCoef, TComp> mul = term * polynomial;
                    for (const auto& [start, evaluated] : evaluatedResults) {
                        if (!start.contains(mul)) {
                            continue;
                        }
                        for (const auto& ev_pol : evaluated) {
                            if (ev_pol.GetLeadingTerm() == mul.GetLeadingTerm()) {
                                if (div.IsOne()) {
                                    return term * ev_pol;
                                } else if (div != term) {
                                    return Simplify(term / div, ev_pol, evaluatedResults);
                                } else {
                                    return ev_pol;
                                }
                            }
                        }
                    }
                }
                return term * polynomial;
            }

            template <typename TCoef, typename TComp>
            void UpdateL(NUtils::TPolynomials<TCoef, TComp>& L, const NUtils::Term& term, const TPolynomialSet<TCoef, TComp>& polynomials, NUtil::TTermSet<TComp>& diff, NUtil::TTermSet<TComp>& done, const TEvaluatedResults<TCoef, TComp>& evaluatedResults) {
                for (const auto& polynomial : polynomials) {
                    const auto& t = polynomial.GetLeadingTerm();
                    if (term.IsDivisibleBy(t)) {
                        NUtils::Polynomial<TCoef, TComp> reducer = Simplify((term / t), polynomial, evaluatedResults);
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
            NUtil::TSymbolicPreprocessingResult<TCoef, TComp> SymbolicPreprocessing(TPairsVector<TCoef, TComp>& selected, const TPolynomialSet<TCoef, TComp>& polynomials, const TEvaluatedResults<TCoef, TComp>& evaluatedResults) {
                NUtils::TPolynomials<TCoef, TComp> L;
                L.reserve(selected.size() * 3);
                for (const auto& pair : selected) {
                    L.push_back(Simplify(pair.GetGlcmTerm() / pair.GetLeftTerm(), pair.GetLeft(), evaluatedResults));
                    L.push_back(Simplify(pair.GetGlcmTerm() / pair.GetRightTerm(), pair.GetRight(), evaluatedResults));
                }

                NUtil::TTermSet<TComp> diff;
                NUtil::TTermSet<TComp> done;
                for (const auto& l : L) {
                    auto it = diff.begin();
                    for (const auto& m : l.GetMonomials()) {
                        it = diff.insert(it, m.GetTerm());
                    }
                    diff.erase(l.GetLeadingTerm());
                    done.insert(l.GetLeadingTerm());
                }

                while(!diff.empty()) {
                    NUtils::Term term = *diff.begin();
                    diff.erase(diff.begin());
                    done.insert(done.begin(), term);
                    UpdateL(L, term, polynomials, diff, done, evaluatedResults);
                }

                return {L, done};
            }

            template <typename TCoef, typename TComp>
            void EraseByLcm(TPairsSet<TCoef, TComp>& pairs_to_check, const NUtils::Polynomial<TCoef, TComp>& f) {
                for (auto it = pairs_to_check.begin(); it != pairs_to_check.end();) {
                    if (gcd(f.GetLeadingTerm(), it->GetRightTerm()).IsOne()) {
                        ++it;
                        continue;
                    }

                    bool deleted = false;
                    for (auto jt = pairs_to_check.begin(); jt != pairs_to_check.end(); ++jt) {
                        if (it == jt) {
                            continue;
                        }
                        if (lcm(f.GetLeadingTerm(), it->GetRightTerm()).IsDivisibleBy(lcm(f.GetLeadingTerm(), jt->GetRightTerm()))) {
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
            void EraseByLead(TPolynomialSet<TCoef, TComp>& polynomials, auto fit) {
                for (auto it = polynomials.begin(); it != polynomials.end(); ) {
                    if (it == fit) {
                        ++it;
                        continue;
                    }
                    if (it->GetLeadingTerm().IsDivisibleBy(fit->GetLeadingTerm())) {
                        it = polynomials.erase(it);
                    } else {
                        ++it;
                    }
                }
            }

            template <typename TCoef, typename TComp>
            bool CheckLcm(const NUtils::CriticalPair<TCoef, TComp>& cp, const NUtils::Term& h) {
                return cp.GetGlcmTerm().IsDivisibleBy(h) &&
                    lcm(cp.GetLeftTerm(), h) != cp.GetGlcmTerm() &&
                    lcm(cp.GetRightTerm(), h) != cp.GetGlcmTerm();
            }

            template <typename TCoef, typename TComp>
            void InsertByLcm(TPairsSet<TCoef, TComp>& old_crit_pairs, TPairsSet<TCoef, TComp>& new_crit_pairs, const NUtils::Polynomial<TCoef, TComp>& f) {
                for (const auto& cp : old_crit_pairs) {
                    if (cp.GetGlcmTerm().IsDivisibleBy(f.GetLeadingTerm())) {
                        continue;
                    }
                    if (CheckLcm(cp, f.GetLeadingTerm())) {
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
                    if (it != fit && !it->GetLeadingTerm().IsDivisibleBy(fit->GetLeadingTerm())) {
                        all_crit.insert(NUtils::CriticalPair(*fit, *it));
                    }
                }

                EraseByLcm(all_crit, *fit);
                InsertByGcd(all_crit, new_crit_pairs, *fit);
                InsertByLcm(old_crit_pairs, new_crit_pairs, *fit);
                old_crit_pairs = std::move(new_crit_pairs);

                EraseByLead(polynomials, fit);
            }

            template <typename TCoef, typename TComp>
            NUtils::TPolynomials<TCoef, TComp> Reduce(TPairsVector<TCoef, TComp>& selected, TPolynomialSet<TCoef, TComp>& polynomials, TEvaluatedResults<TCoef, TComp>& evaluatedResults) {
                NUtil::TSymbolicPreprocessingResult<TCoef, TComp> L = SymbolicPreprocessing(selected, polynomials, evaluatedResults);
                NUtil::PolynomialsPair<TCoef, TComp> result = NUtil::MatrixReduction(L);
                TPolynomialSet<TCoef, TComp> start;
                for (const auto& pol : L.first) {
                    start.insert(pol);
                }
                evaluatedResults.emplace_back(start, result.second);
                return result.first;
            }

            template <typename TCoef, typename TComp>
            void UpdateBasis(const TPolynomialSet<TCoef, TComp>& polynomials, NUtils::TPolynomials<TCoef, TComp>& F) {
                F.clear();
                F.reserve(polynomials.size());
                for (const auto& x : polynomials) {
                    F.push_back(x);
                }
            }

            template <typename TCoef, typename TComp>
            void FindGroebnerBasis(NUtils::TPolynomials<TCoef, TComp>& F) {
                TPolynomialSet<TCoef, TComp> polynomials;
                TPairsSet<TCoef, TComp> pairs_to_check;
                TEvaluatedResults<TCoef, TComp> evaluatedResults;
                for (auto& f : F) {
                    UpdateCriticalPairs(polynomials, pairs_to_check, f);
                }

                while(!pairs_to_check.empty()) {
                    TPairsVector<TCoef, TComp> selection_group = Select(pairs_to_check);
                    NUtils::TPolynomials<TCoef, TComp> G = Reduce(selection_group, polynomials, evaluatedResults);
                    for (auto& g : G) {
                        UpdateCriticalPairs(polynomials, pairs_to_check, g);
                    }
                }
                UpdateBasis(polynomials, F);
            }
        }
    }
}
