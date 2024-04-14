#pragma once
#include "term.h"

namespace NUtils {
    template <typename TCoef>
    class Monomial {
    public:
        Monomial() = default;

        Monomial(TTerm&& term, TCoef coef = 1)
        : term_(std::move(term))
        , coef_(coef)
        {
        }

        const TTerm& GetTerm() const noexcept {
            return term_;
        }

        const TCoef& GetCoef() const noexcept {
            return coef_;
        }

        int64_t GetNumerator() const noexcept {
            return coef_.GetNumerator();
        }

        int64_t GetDenominator() const noexcept {
            return coef_.GetDenominator();
        }

        void AddCoef(const Monomial& other) noexcept {
            coef_ += other.GetCoef();
        }

        void SubCoef(const Monomial& other) noexcept {
            coef_ -= other.GetCoef();
        }

        Monomial operator+() const noexcept {
            return *this;
        }

        Monomial operator-() const noexcept {
            TTerm tmpTerm = term_;
            return Monomial(std::move(tmpTerm), -coef_);
        }

        friend bool operator==(const Monomial& left, const Monomial& right) noexcept {
            if (left.coef_ != right.coef_) {
                return false;
            }
            if (left.coef_ == 0) {
                return true;
            }
            return left.term_ == right.term_;
        }

        friend bool operator!=(const Monomial& left, const Monomial& right) noexcept {
            return !(left == right);
        }

        Monomial& operator*=(const Monomial& other) noexcept {
            coef_ *= other.GetCoef();
            term_ *= other.GetTerm();
            return *this;
        }

        friend Monomial operator*(Monomial left, const Monomial& right) noexcept {
            left *= right;
            return left;
        }

        Monomial& operator*=(const TCoef& coef) noexcept {
            coef_ *= coef;
            return *this;
        }

        friend Monomial operator*(Monomial left, const TCoef& coef) noexcept {
            left *= coef;
            return left;
        }

        Monomial& operator/=(const TCoef& coef) {
            assert(coef != 0);
            coef_ /= coef;
            return *this;
        }

        friend Monomial operator/(Monomial left, const TCoef& coef) {
            left /= coef;
            return left;
        }

        Monomial& operator/=(const Monomial& other) {
            assert(other.GetCoef() != 0);
            coef_ /= other.GetCoef();
            term_ /= other.GetTerm();
            return *this;
        }

        friend Monomial operator/(Monomial left, const Monomial& right) {
            left /= right;
            return left;
        }

        friend bool operator<(const Monomial& left, const Monomial& right) noexcept {
            return std::lexicographical_compare(left.term_.begin(), left.term_.end(), right.term_.begin(), right.term_.end());
        }

        friend bool operator>(const Monomial& left, const Monomial& right) noexcept {
            return right < left;
        }

        friend bool operator<=(const Monomial& left, const Monomial& right) noexcept {
            return !(left > right);
        }

        friend bool operator>=(const Monomial& left, const Monomial& right) noexcept {
            return !(left < right);
        }

        friend std::ostream& operator<<(std::ostream& out, const Monomial& monomial) noexcept {
            if (monomial.GetCoef() == 1) {
                out << monomial.GetTerm();
            } else if (monomial.GetCoef() == -1) {
                out << " - " << monomial.GetTerm();
            } else if (monomial.GetCoef() != 0) {
                out << monomial.GetCoef() << monomial.GetTerm();
            }
            return out;
        }

    private:
        TTerm term_;
        TCoef coef_;
    };
}
