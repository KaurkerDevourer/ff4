#include "monomial.h"
#include <cassert>

namespace NUtils {
    Monomial::Monomial(TTerm term, Rational coef)
    : term_(std::move(term))
    , coef_(coef)
    {
    }

    const TTerm& Monomial::GetTerm() const noexcept {
        return term_;
    }

    const Rational& Monomial::GetCoef() const noexcept {
        return coef_;
    }

    int64_t Monomial::GetNumerator() const noexcept {
        return coef_.GetNumerator();
    }

    int64_t Monomial::GetDenominator() const noexcept {
        return coef_.GetDenominator();
    }

    void Monomial::AddCoef(const Monomial& other) noexcept {
        coef_ += other.GetCoef();
    }

    void Monomial::SubCoef(const Monomial& other) noexcept {
        coef_ -= other.GetCoef();
    }

    bool operator<(const Monomial& left, const Monomial& right) noexcept {
        return std::lexicographical_compare(left.term_.begin(), left.term_.end(), right.term_.begin(), right.term_.end());
    }

    bool operator>(const Monomial& left, const Monomial& right) noexcept {
        return right < left;
    }

    bool operator<=(const Monomial& left, const Monomial& right) noexcept {
        return !(left > right);
    }

    bool operator>=(const Monomial& left, const Monomial& right) noexcept {
        return !(left < right);
    }

    bool operator==(const Monomial& left, const Monomial& right) noexcept {
        if (left.coef_ != right.coef_) {
            return false;
        }
        return left.term_ == right.term_;
    }

    bool operator!=(const Monomial& left, const Monomial& right) noexcept {
        return !(left == right);
    }

    Monomial Monomial::operator+() const noexcept {
        return *this;
    }

    Monomial Monomial::operator-() const noexcept {
        return Monomial(term_, -coef_);
    }

    Monomial& Monomial::operator*=(const Monomial& other) noexcept {
        coef_ *= other.GetCoef();
        term_ *= other.GetTerm();
        return *this;
    }

    Monomial operator*(Monomial left, const Monomial& right) noexcept {
        left *= right;
        return left;
    }

    Monomial& Monomial::operator*=(const Rational& coef) noexcept {
        coef_ *= coef;
        return *this;
    }

    Monomial operator*(Monomial left, const Rational& coef) noexcept {
        left *= coef;
        return left;
    }

    Monomial& Monomial::operator/=(const Rational& coef) {
        assert(coef != 0);
        coef_ /= coef;
        return *this;
    }

    Monomial operator/(Monomial left, const Rational& coef) {
        left /= coef;
        return left;
    }

    Monomial& Monomial::operator/=(const Monomial& other) {
        assert(other.GetCoef() != 0);
        coef_ /= other.GetCoef();
        term_ /= other.GetTerm();
        return *this;
    }

    Monomial operator/(Monomial left, const Monomial& right) {
        left /= right;
        return left;
    }

    std::ostream& operator<<(std::ostream& out, const Monomial& monomial) noexcept {
        out << monomial.GetCoef() << monomial.GetTerm();
        return out;
    }
}
