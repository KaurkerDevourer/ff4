#include "polynomial.h"

namespace NUtils {
    Polynomial::Polynomial(TMonomials&& monomials)
    : monomials_(monomials)
    {
        Normalize(monomials_);
    }

    TMonomials& Polynomial::GetMonomials() const noexcept {
        return monomials_;
    }

    TMonomial& Polynomial::GetHeadMonomial() const noexcept {
        return monomials_[0];
    }

    Polynomial& Polynomial::operator+() const noexcept {
        return *this;
    }

    bool Polynomial::IsZero() const noexcept {
        return monomials_.empty();
    }

    Polynomial Polynomial::operator-() const noexcept {
        TMonomials monomials(monomials_.size());
        for (size_t i = 0; i < monomials.size(); i++) {
            monomials[i] = -monomials_[i];
        }
        return Polynomial(std::move(monomials));
    }

    Polynomial& Polynomial::operator+=(const Polynomial& other) noexcept {
        TMonomials monomials;
        monomials.reserve(monomials_.size() + other.monomials_.size());
        size_t i = 0;
        size_t j = 0;
        while(i != monomials_.size() && j != other.monomials_.size()) {
            if (monomials_[i] < other.monomials_[j]) {
                monomials.push_back(monomials_[i]);
                i++;
            } else if (monomials_[i] > other.monomials_[j]) {
                monomials.push_back(other.monomials_[j]);
                j++;
            } else {
                monomials_[i].coef_ += other.monomials_[j].coef_;
                if (monomials_[i].coef != 0) {
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
        while(j != monomials_.size()) {
            monomials.push_back(other.monomials_[j]);
            j++;
        }
        monomials.shrink_to_fit();
        monomials_ = std::move(monomials);
        return *this;
    }

    Polynomial& Polynomial::operator-=(const Polynomial& other) noexcept {
        TMonomials monomials;
        monomials.reserve(monomials_.size() + other.monomials_.size());
        size_t i = 0;
        size_t j = 0;
        while(i != monomials_.size() && j != other.monomials_.size()) {
            if (monomials_[i] < other.monomials_[j]) {
                monomials.push_back(monomials_[i]);
                i++;
            } else if (monomials_[i] > other.monomials_[j]) {
                monomials.push_back(-other.monomials_[j]);
                j++;
            } else {
                monomials_[i].coef_ -= other.monomials_[j].coef_;
                if (monomials_[i].coef != 0) {
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
        while(j != monomials_.size()) {
            monomials.push_back(-other.monomials_[j]);
            j++;
        }
        monomials.shrink_to_fit();
        monomials_ = std::move(monomials);
        return *this;
    }

    Polynomial& Polynomial::operator*=(const Monomial& monomial) noexcept {
        for (size_t i = 0; i < monomials_.size(); i++) {
            monomials_[i] *= monomial;
        }
        return *this;
    }

    bool Polynomial::ReduceBy(const std::vector<Polynomial>& F) noexcept {
        bool changed = true;
        while(changed) {
            changed = false;
            for (size_t i = 0; i < F.size(); i++) {
                while (GetHeadMonomial().IsDivisibleBy(F[i].GetHeadMonomial())) {
                    *this -= (GetHeadMonomial() / F[i].GetHeadMonomial()) * F[i];
                    changed = true;
                }
            }
        }
        return IsZero();
    }

    //private
    void Polynomial::Normalize(TMonomials& monomials) {
        bool isSorted = true;
        size_t cnt = 0;
        for (size_t i = 0; i < monomials.size();) {
            if (monomials[i] == {}) {
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
            std::erase(monomials.end() - cnt, monomials.end());
        }

        if (!isSorted) {
            std::sort(monomials.rbegin(), monomials.rend());
        }
    }
}