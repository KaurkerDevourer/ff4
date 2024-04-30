#pragma once
#include "../util/critical_pair.hpp"
#include "../util/comp.hpp"
#include <cassert>
#include <set>

namespace NAlgo {
    namespace F4 {
        using namespace NUtils;

        template<typename TCoef>
        using TPairsSet = std::set<CriticalPair<TCoef>, DexComp>;

        using TDiffSet = std::set<TTerm, TTermReverseComp>;

        template<typename TCoef>
        using TSymbolicPreprocessingResult = std::pair<Polynomials<TCoef>, TDiffSet>

        template <typename TCoef>
        bool CheckProductCreteria(const Polynomial<TCoef>& a, const Polynomial<TCoef>& b) {
            const Monomial<TCoef>& am = a.GetHeadMonomial();
            const Monomial<TCoef>& bm = b.GetHeadMonomial();
            const TTerm t = lcm(am.GetTerm(), bm.GetTerm());
            return (t == am.GetTerm() * bm.GetTerm());
        }

        // https://apmi.bsu.by/assets/files/agievich/em-atk.pdf
        template <typename TCoef>
        TPairsSet<TCoef> GetPairsToCheckWithCriterias(const TPolynomials<TCoef>& polynomials) {
            TPairsSet<TCoef> pairs_to_check;
            for (size_t i = 0; i < polynomials.size(); i++) {
                for (size_t j = i + 1; j < polynomials.size(); j++) {
                    if (CheckProductCreteria(polynomials[i], polynomials[j])) {
                        continue;
                    }
                    pairs_to_check.insert(CriticalPair(&polynomials, i, j));
                }
            }
            return pairs_to_check;
        }

        template <typename TCoef>
        TPairsSet<TCoef> Select(TPairsSet<TCoef>& pairs_to_check) {
            TPairsSet<TCoef> selectionGroup;
            uint64_t value = pairs_to_check.begin()->GetDegree();
            while(pairs_to_check.size() && pairs_to_check.begin()->GetDegree() == value) {
                selectionGroup.insert(selectionGroup.begin(), *pairs_to_check.begin());
                pairs_to_check.erase(pairs_to_check.begin());
            }
            return selectionGroup;
        }

        template <typename TCoef>
        void UpdateL(TPolynomials<TCoef>& L, const TTerm& term, const TPolynomials<TCoef>& F, TDiffSet& diff, TDiffSet& done) {
            for (const auto& polynomial : F) {
                const auto& t = polynomial.GetHeadMonomial().GetTerm();
                if (term.IsDivisibleBy(t)) {
                    Polynomial<TCoef> reducer = (term / t) * polynomial;
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

        template <typename TCoef>
        TSymbolicPreprocessingResult<TCoef> SymbolicPreprocessing(TPairsSet<TCoef>& selected, const TPolynomials<TCoef>& F) {
            TPolynomials<TCoef> L;
            L.reserve(selected.size() * 2);
            for (const auto& pair : selected) {
                L.push_back(pair.GetGlcmTerm() / pair.GetLeftTerm() * pair.GetLeft()); // add left
                L.push_back(pair.GetGlcmTerm() / pair.GetRightTerm() * pair.GetRight()); // add right
            }
            TDiffSet diff;
            for (const auto& l : L) {
                auto it = diff.begin();
                for (const auto& m : l.GetMonomials()) {
                    it = diff.insert(it, m.GetTerm());
                }
            }

            TDiffSet done;
            for (const auto& l : L) {
                diff.erase(l.GetHeadMonomial().GetTerm());
                done.insert(l.GetHeadMonomial().GetTerm());
            }

            while(!diff.empty()) {
                const TTerm& term = *diff.begin();
                diff.erase(diff.begin());
                done.insert(done.begin(), term);
                UpdateL(L, term, F, diff, done);
            }

            return {L, done};
        }

        template <typename TCoef>
        void Reduce(TPairsSet<TCoef>& selected, TPolynomials<TCoef>& F) {
            TSymbolicPreprocessingResult<TCoef> L = SymbolicPreprocessing(selected, F);

        }

        template <typename TCoef>
        void FindGroebnerBasis(TPolynomials<TCoef>& F) {
            TPairsSet<TCoef> pairs_to_check = GetPairsToCheckWithCriterias(F);

            while(!pairs_to_check.empty()) {
                TPairsSet<TCoef> selection_group = Select(pairs_to_check);
                Reduce(selection_group, F);
            };
        }
    }
}
