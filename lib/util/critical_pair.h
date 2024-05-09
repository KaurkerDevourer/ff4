#pragma once

#include "polynomial.h"

namespace FF4 {
    namespace NUtils {

        template<typename TCoef, typename TComp>
        class CriticalPair {
            public:
                CriticalPair(const NUtils::Polynomial<TCoef, TComp>& left, const NUtils::Polynomial<TCoef, TComp>& right)
                    : left_(left)
                    , right_(right)
                    , Glcm_(Monomial(lcm(left.GetLeadingTerm(), right.GetLeadingTerm()), TCoef(1)))
                    , degree_(Glcm_.GetTerm().TotalDegree())
                {
                }

                uint64_t TotalDegree() const noexcept {
                    return degree_;
                }

                const Monomial<TCoef>& GetGlcm() const noexcept {
                    return Glcm_;
                }

                const Term& GetGlcmTerm() const noexcept {
                    return Glcm_.GetTerm();
                }

                const Polynomial<TCoef, TComp>& GetLeft() const noexcept {
                    return left_;
                }

                const Polynomial<TCoef, TComp>& GetRight() const noexcept {
                    return right_;
                }

                const Term& GetLeftTerm() const noexcept {
                    return left_.GetLeadingTerm();
                }

                const Term& GetRightTerm() const noexcept {
                    return right_.GetLeadingTerm();
                }

                friend std::ostream& operator<<(std::ostream& out, const CriticalPair& cp) noexcept {
                    return out << cp.degree_ << " | " << cp.GetLeftTerm() << " " << cp.GetRightTerm() << " : " << cp.GetGlcmTerm();
                }

            private:
                const NUtils::Polynomial<TCoef, TComp>& left_;
                const NUtils::Polynomial<TCoef, TComp>& right_;
                Monomial<TCoef> Glcm_;
                uint64_t degree_;
        };
    }
}
