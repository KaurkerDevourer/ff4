#pragma once

#include "polynomial.hpp"

namespace NUtils {

    template<typename TCoef, typename TComp>
    class CriticalPair {
        public:
            CriticalPair() = default;
            CriticalPair(const TPolynomials<TCoef, TComp>* F, size_t i, size_t j)
                : F_(F)
                , left_idx_(i)
                , right_idx_(j)
                , Glcm_(Monomial(lcm((*F)[i].GetHeadMonomial().GetTerm(), (*F)[j].GetHeadMonomial().GetTerm()), TCoef(1)))
                , degree_(Glcm_.GetTerm().GetDegree())
            {
            }

            uint64_t GetDegree() const noexcept {
                return degree_;
            }

            const Monomial<TCoef>& GetGlcm() const noexcept {
                return Glcm_;
            }

            const TTerm& GetGlcmTerm() const noexcept {
                return Glcm_.GetTerm();
            }

            const Polynomial<TCoef, TComp>& GetLeft() const noexcept {
                return (*F_)[left_idx_];
            }

            const Polynomial<TCoef, TComp>& GetRight() const noexcept {
                return (*F_)[right_idx_];
            }

            const TTerm& GetLeftTerm() const noexcept {
                return (*F_)[left_idx_].GetHeadMonomial().GetTerm();
            }

            const TTerm& GetRightTerm() const noexcept {
                return (*F_)[right_idx_].GetHeadMonomial().GetTerm();
            }

        private:
            const TPolynomials<TCoef, TComp>* F_;
            size_t left_idx_;
            size_t right_idx_;
            Monomial<TCoef> Glcm_;
            uint64_t degree_;
    };
}
