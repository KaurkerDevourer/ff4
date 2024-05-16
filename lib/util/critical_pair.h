#pragma once

#include "polynomial.h"

namespace FF4 {
    namespace NUtils {

        template<typename TCoef, typename TComp>
        class CriticalPair {
            public:
                CriticalPair(const Polynomial<TCoef, TComp>& left, const Polynomial<TCoef, TComp>& right)
                    : left_(left)
                    , right_(right)
                    , Glcm_(lcm(left.GetLeadingTerm(), right.GetLeadingTerm()))
                    , degree_(Glcm_.TotalDegree())
                {
                }

                Term::Degree TotalDegree() const noexcept {
                    return degree_;
                }

                const Term& GetGlcm() const noexcept {
                    return Glcm_;
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
                    return out << cp.degree_ << " | " << cp.GetLeftTerm() << " " << cp.GetRightTerm() << " : " << cp.GetGlcm();
                }

            private:
                const Polynomial<TCoef, TComp>& left_;
                const Polynomial<TCoef, TComp>& right_;
                Term Glcm_;
                Term::Degree degree_;
        };
    }
}
