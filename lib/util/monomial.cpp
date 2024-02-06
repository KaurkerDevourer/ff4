#include "monomial.h"

namespace NUtils {
    Monomial::Monomial(std::vector<int64_t>&& term, Rational coef_)
    : term_(std::move(term))
    , coef_(coef)
    {
        Normalize(term_);
    }

    TTerm& Monomial::GetTerm() const noexcept {
        return term_;
    }

    Rational& Monomial::GetCoef() const noexcept {
        return coef_;
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
            if (other_term[i] > term[i]) {
                return false;
            }
        }
        return true;
    };

    bool operator<(const Monomial& left, const Monomial& right) noexcept {
        for (size_t i = 0; i < std::min(left.term_.size(), right.term_.size()); i++) {
            if (left.term_[i] != right.term_[i]) {
                return left.term[i] < right.term[i];
            }
        }
        return left.term_.size() < right.term_.size();
    }

    bool operator>(const Monomial& left, const Monomial& right) noexcept {
        for (size_t i = 0; i < std::min(left.term_.size(), right.term_.size()); i++) {
            if (left.term_[i] != right.term_[i]) {
                return left.term[i] > right.term[i];
            }
        }
        return left.term_.size() > right.term_.size();
    }

    bool operator<=(const Monomial& left, const Monomial& right) noexcept {
        for (size_t i = 0; i < std::min(left.term_.size(), right.term_.size()); i++) {
            if (left.term_[i] != right.term_[i]) {
                return left.term[i] < right.term[i];
            }
        }
        return left.term_.size() <= right.term_.size();
    }

    bool operator>=(const Monomial& left, const Monomial& right) noexcept {
        for (size_t i = 0; i < std::min(left.term_.size(), right.term_.size()); i++) {
            if (left.term_[i] != right.term_[i]) {
                return left.term[i] > right.term[i];
            }
        }
        return left.term_.size() >= right.term_.size();
    }

    bool operator==(const Monomial& left, const Monomial& right) noexcept {
        if (left.coef_ != right.coef_) {
            return false;
        }
        if (left.term_.size() != right.term_.size()) {
            return false;
        }
        for (size_t i = 0; i < std::min(left.term_.size(), right.term_.size()); i++) {
            if (left.term_[i] != right.term_[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Monomial& left, const Monomial& right) noexcept {
        if (left.coef_ != right.coef_) {
            return true;
        }
        if (left.term_.size() != right.term_.size()) {
            return true;
        }
        for (size_t i = 0; i < std::min(left.term_.size(), right.term_.size()); i++) {
            if (left.term_[i] != right.term_[i]) {
                return true;
            }
        }
        return false;
    }

    Monomial& Monomial::operator+() const noexcept {
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

    Monomial gcd(Monomial left, const Monomial& right) noexcept {
        size_t sz = std::min(left.term_.size(), right.term_.size());
        TTerm term(sz);
        for (size_t i = 0; i < sz; i++) {
            term[i] = std::min(left.term_[i], right.term_[i]);
        }
        Normalize(term);
        return Monomial(std::move(term), 1);
    }

    Monomial lcm(Monomial left, const Monomial& right) noexcept {
        size_t sz = std::min(left.term_.size(), right.term_.size());
        size_t lcm_sz = std::max(left.term_.size(), right.term_.size());
        TTerm term(lcm_sz);
        for (size_t i = 0; i < sz; i++) {
            term[i] = std::max(left.term_[i], right.term_[i]);
        }
        for (size_t i = left.term_.size(); i < lcm_sz; i++) {
            term[i] = right.term_[i];
        }
        for (size_t i = right.term_.size(); i < lcm_sz; i++) {
            term[i] = left.term_[i];
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