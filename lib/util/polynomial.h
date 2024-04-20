#pragma once
#include "monomial.h"
#include "comp.h"
#include <vector>
#include <queue>

namespace NUtils {

    template <typename TCoef>
    using TMonomials = std::vector<Monomial<TCoef>>;

    template <typename TCoef, typename TComp>
    void DebugCheckPolynomialIsCorrect(TMonomials<TCoef>& monomials) {
        bool isSorted = true;
        size_t cnt = 0;
        for (size_t i = 0; i < monomials.size(); i++) {
            if (monomials[i].GetCoef() == 0) {
                cnt++;
                continue;
            }
            monomials[i] = monomials[i + cnt];
            if (i != 0 && TComp()(monomials[i - 1], monomials[i])) {
                isSorted = false;
            }
            i++;
        }
        if (cnt) {
            std::cout << "Zero coef" << std::endl;
            monomials.erase(monomials.end() - cnt, monomials.end());
        }

        if (!isSorted) {
            std::cout << "!IsSorted" << std::endl;
            std::sort(monomials.rbegin(), monomials.rend(), TComp());
        }
    }

    template <typename TCoef, typename TComp = LexComp>
    class Polynomial {
    public:
        Polynomial() = default;

        Polynomial(TMonomials<TCoef>&& monomials)
        : monomials_(std::move(monomials))
        {
            #ifdef NDEBUG
                // nondebug
            #else
                DebugCheckPolynomialIsCorrect<TCoef, TComp>(monomials_);
            #endif
        }

        const TMonomials<TCoef>& GetMonomials() const noexcept {
            return monomials_;
        }

        const Monomial<TCoef>& GetHeadMonomial() const noexcept {
            return monomials_[0];
        }

        bool IsZero() const noexcept {
            return monomials_.empty();
        }

        Polynomial operator+() const noexcept {
            return *this;
        }

        Polynomial operator-() const noexcept {
            TMonomials<TCoef> monomials(monomials_.size());
            for (size_t i = 0; i < monomials_.size(); i++) {
                monomials[i] = -monomials_[i];
            }
            return Polynomial(std::move(monomials));
        }

        friend bool operator==(const Polynomial& left, const Polynomial& right) noexcept {
            const TMonomials<TCoef>& lefTMonomials = left.GetMonomials();
            const TMonomials<TCoef>& righTMonomials = right.GetMonomials();
            if (lefTMonomials.size() != righTMonomials.size()) {
                return false;
            }
            for (size_t i = 0; i < lefTMonomials.size(); i++) {
                if (lefTMonomials[i] != righTMonomials[i]) {
                    return false;
                }
            }
            return true;
        }

        friend bool operator!=(const Polynomial& left, const Polynomial& right) noexcept {
            return !(left == right);
        }

        Polynomial& operator+=(const Polynomial& other) noexcept {
            TMonomials<TCoef> monomials;
            monomials.reserve(monomials_.size() + other.monomials_.size());
            size_t i = 0;
            size_t j = 0;
            while(i != monomials_.size() && j != other.monomials_.size()) {
                if (TComp()(other.monomials_[j], monomials_[i])) {
                    monomials.push_back(monomials_[i]);
                    i++;
                } else if (TComp()(monomials_[i], other.monomials_[j])) {
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

        friend Polynomial operator+(Polynomial left, const Polynomial& right) noexcept {
            left += right;
            return left;
        }

        Polynomial& operator-=(const Polynomial& other) noexcept {
            *this += -other;
            return *this;
        }

        friend Polynomial operator-(Polynomial left, const Polynomial& right) noexcept {
            left -= right;
            return left;
        }

        Polynomial& operator*=(const Monomial<TCoef>& monomial) noexcept {
            for (size_t i = 0; i < monomials_.size(); i++) {
                monomials_[i] *= monomial;
            }
            return *this;
        }

        friend Polynomial operator*(Polynomial left, const Monomial<TCoef>& right) noexcept {
            left *= right;
            return left;
        }

        friend Polynomial operator*(const Monomial<TCoef>& left, Polynomial right) noexcept {
            right *= left;
            return right;
        }

        Polynomial& operator*=(const TCoef& coef) noexcept {
            for (size_t i = 0; i < monomials_.size(); i++) {
                monomials_[i] *= coef;
            }
            return *this;
        }

        friend Polynomial operator*(Polynomial left, const TCoef& right) noexcept {
            left *= right;
            return left;
        }

        friend Polynomial operator*(const TCoef& left, Polynomial right) noexcept {
            right *= left;
            return right;
        }

        Polynomial& operator*=(const Polynomial& polynomial) noexcept {
            const TMonomials<TCoef>& monomials = polynomial.GetMonomials();
            Polynomial newPolynomial;
            for (size_t i = 0; i < monomials.size(); i++) {
                Polynomial tmpMonomials = (*this);
                tmpMonomials *= monomials[i];
                newPolynomial += tmpMonomials;
            }
            *this = std::move(newPolynomial);
            return *this;
        }

        friend Polynomial operator*(Polynomial left, const Polynomial& right) noexcept {
            left *= right;
            return left;
        }

        Polynomial& operator/=(const TCoef& coef) noexcept {
            for (size_t i = 0; i < monomials_.size(); i++) {
                monomials_[i] /= coef;
            }
            return *this;
        }

        friend Polynomial operator/(Polynomial left, const TCoef& right) noexcept {
            left /= right;
            return left;
        }

        friend std::ostream& operator<<(std::ostream& out, const Polynomial& polynomial) noexcept {
            const TMonomials<TCoef>& monomials = polynomial.GetMonomials();
            for (size_t i = 0; i < monomials.size(); i++) {
                if (i > 0 && monomials[i].GetCoef().MoreThanZero()) {
                    out << " + ";
                }
                out << monomials[i];
            }
            return out;
        }

    private:
        TMonomials<TCoef> monomials_;
    };

    template <typename TCoef>
    using TPolynomials = std::vector<Polynomial<TCoef>>;

    std::queue<std::pair<size_t, size_t>> GetPairsToCheck(size_t sz) {
        std::queue<std::pair<size_t, size_t>> pairs_to_check;
        for (size_t i = 0; i < sz; i++) {
            for (size_t j = i + 1; j < sz; j++) {
                pairs_to_check.push({i, j});
            }
        }
        return pairs_to_check;
    }

    template <typename TCoef>
    std::ostream& operator<<(std::ostream& out, const TPolynomials<TCoef>& polynomials) noexcept {
        out << "{\n";
        for (const Polynomial<TCoef>& polynomial : polynomials) {
            out << polynomial << '\n';
        }
        return out << "}\n";
    }
}
