#include "monomial.h"

namespace NUtils {
    Monomial::Monomial(std::vector<int64_t> term, Rational coef)
    : term_(term)
    , coef_(coef)
    {
        Normalize(term_);
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

    bool Monomial::IsDivisibleBy(const Monomial& other) const noexcept {
        if (other.GetCoef() == 0) {
            return false;
        }
        const TTerm& other_term = other.GetTerm();
        if (other_term.size() > term_.size()) {
            return false;
        }
        for (size_t i = 0; i < other_term.size(); i++) {
            if (other_term[i] > term_[i]) {
                return false;
            }
        }
        return true;
    };

    bool operator<(const Monomial& left, const Monomial& right) noexcept {
        std::lexicographical_compare(left.term_.begin(), left.term_.end(), right.term_.begin(), right.term_.end());
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
        if (other.term_.size() > term_.size()) {
            term_.resize(other.term_.size());
        }
        for (size_t i = 0; i < other.term_.size(); i++) {
            term_[i] += other.term_[i];
        }
        Normalize(term_);
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

    Monomial& Monomial::operator/=(const Monomial& other) {
        //assert(other != 0);
        coef_ /= other.GetCoef();
        if (other.term_.size() > term_.size()) {
            term_.resize(other.term_.size());
        }
        for (size_t i = 0; i < other.term_.size(); i++) {
            term_[i] -= other.term_[i];
        }
        Normalize(term_);
        return *this;
    }

    Monomial operator/(Monomial left, const Monomial& right) {
        left /= right;
        return left;
    }

    Monomial gcd(const Monomial& left, const Monomial& right) noexcept {
        const TTerm& leftTerm = left.GetTerm();
        const TTerm& rightTerm = right.GetTerm();

        size_t sz = std::min(leftTerm.size(), rightTerm.size());
        TTerm term(sz);
        for (size_t i = 0; i < sz; i++) {
            term[i] = std::min(leftTerm[i], rightTerm[i]);
        }
        return Monomial(std::move(term), 1);
    }

    Monomial lcm(const Monomial& left, const Monomial& right) noexcept {
        const TTerm& leftTerm = left.GetTerm();
        const TTerm& rightTerm = right.GetTerm();
    
        size_t sz = std::min(leftTerm.size(), rightTerm.size());
        size_t lcm_sz = std::max(leftTerm.size(), rightTerm.size());
        TTerm term(lcm_sz);
        for (size_t i = 0; i < sz; i++) {
            term[i] = std::max(leftTerm[i], rightTerm[i]);
        }
        for (size_t i = leftTerm.size(); i < lcm_sz; i++) {
            term[i] = rightTerm[i];
        }
        for (size_t i = rightTerm.size(); i < lcm_sz; i++) {
            term[i] = leftTerm[i];
        }
        return Monomial(std::move(term), 1);
    }

    // private
    void Monomial::Normalize(TTerm& term) {
        while(term.size() && term.back() == 0) {
            term.pop_back();
        }
    }
}
