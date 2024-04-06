#include "polynomial.h"
#include <algorithm>
#include <vector>

namespace NUtils {
    Polynomial::Polynomial(TMonomials&& monomials)
    : monomials_(Normalize(std::move(monomials)))
    {
    }

    const TMonomials& Polynomial::GetMonomials() const noexcept {
        return monomials_;
    }

    const Monomial& Polynomial::GetHeadMonomial() const noexcept {
        return monomials_[0];
    }

    bool Polynomial::IsZero() const noexcept {
        return monomials_.empty();
    }

    Polynomial Polynomial::operator+() const noexcept {
        return *this;
    }

    Polynomial Polynomial::operator-() const noexcept {
        TMonomials monomials(monomials_.size());
        for (size_t i = 0; i < monomials_.size(); i++) {
            monomials[i] = -monomials_[i];
        }
        return Polynomial(std::move(monomials));
    }

    bool operator==(const Polynomial& left, const Polynomial& right) noexcept {
        const TMonomials& leftMonomials = left.GetMonomials();
        const TMonomials& rightMonomials = right.GetMonomials();
        if (leftMonomials.size() != rightMonomials.size()) {
            return false;
        }
        for (size_t i = 0; i < leftMonomials.size(); i++) {
            if (leftMonomials[i] != rightMonomials[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Polynomial& left, const Polynomial& right) noexcept {
        return !(left == right);
    }

    Polynomial& Polynomial::operator+=(const Polynomial& other) noexcept {
        TMonomials monomials;
        monomials.reserve(monomials_.size() + other.monomials_.size());
        size_t i = 0;
        size_t j = 0;
        while(i != monomials_.size() && j != other.monomials_.size()) {
            if (monomials_[i] > other.monomials_[j]) {
                monomials.push_back(monomials_[i]);
                i++;
            } else if (monomials_[i] < other.monomials_[j]) {
                monomials.push_back(other.monomials_[j]);
                j++;
            } else {
                monomials_[i].AddCoef(other.monomials_[j]);
                if (monomials_[i].GetCoef() != 0) {
                    monomials.push_back(monomials_[i]);
                }
                i++;
                j++;
            }
        }
        while(i != monomials_.size()) {
            monomials.push_back(monomials_[i]);
            i++;
        }
        while(j != other.monomials_.size()) {
            monomials.push_back(other.monomials_[j]);
            j++;
        }
        monomials_ = std::move(monomials);
        return *this;
    }

    Polynomial& Polynomial::operator-=(const Polynomial& other) noexcept {
        *this += -other;
        return *this;
    }

    Polynomial& Polynomial::operator*=(const Monomial& monomial) noexcept {
        for (size_t i = 0; i < monomials_.size(); i++) {
            monomials_[i] *= monomial;
        }
        return *this;
    }

    Polynomial operator*(Polynomial left, const Monomial& right) noexcept {
        left *= right;
        return left;
    }

    Polynomial operator*(const Monomial& left, Polynomial right) noexcept {
        right *= left;
        return right;
    }

    Polynomial operator-(Polynomial left, const Polynomial& right) noexcept {
        left -= right;
        return left;
    }

    Polynomial operator+(Polynomial left, const Polynomial& right) noexcept {
        left += right;
        return left;
    }

    Polynomial& Polynomial::operator*=(const Rational& coef) noexcept {
        for (size_t i = 0; i < monomials_.size(); i++) {
            monomials_[i] *= coef;
        }
        return *this;
    }

    Polynomial operator*(Polynomial left, const Rational& right) noexcept {
        left *= right;
        return left;
    }

    Polynomial operator*(const Rational& left, Polynomial right) noexcept {
        right *= left;
        return right;
    }

    Polynomial& Polynomial::operator/=(const Rational& coef) noexcept {
        for (size_t i = 0; i < monomials_.size(); i++) {
            monomials_[i] /= coef;
        }
        return *this;
    }

    Polynomial operator/(Polynomial left, const Rational& right) noexcept {
        left *= right;
        return left;
    }

    Polynomial& Polynomial::operator*=(const Polynomial& polynomial) noexcept {
        const TMonomials& monomials = polynomial.GetMonomials();
        Polynomial newPolynomial;
        for (size_t i = 0; i < monomials.size(); i++) {
            Polynomial tmpMonomials = (*this);
            tmpMonomials *= monomials[i];
            newPolynomial += tmpMonomials;
        }
        *this = std::move(newPolynomial);
        return *this;
    }

    Polynomial operator*(Polynomial left, const Polynomial& right) noexcept {
        left *= right;
        return left;
    }

    std::ostream& operator<<(std::ostream& out, const Polynomial& polynomial) noexcept {
        const TMonomials& monomials = polynomial.GetMonomials();
        for (size_t i = 0; i < monomials.size(); i++) {
            if (i > 0 || monomials[i].GetCoef() > 0) {
                out << " + ";
            }
            out << monomials[i];
        }
        return out;
    }

    //private
    TMonomials&& Polynomial::Normalize(TMonomials&& monomials) {
        bool isSorted = true;
        size_t cnt = 0;
        for (size_t i = 0; i < monomials.size();) {
            if (monomials[i].GetCoef().GetNumerator() == 0) {
                cnt++;
                continue;
            }
            monomials[i] = monomials[i + cnt];
            if (i != 0 && monomials[i] < monomials[i - 1]) {
                isSorted = false;
            }
            i++;
        }
        if (cnt) {
            monomials.erase(monomials.end() - cnt, monomials.end());
        }

        if (!isSorted) {
            std::sort(monomials.rbegin(), monomials.rend());
        }
        return std::move(monomials);
    }
}
