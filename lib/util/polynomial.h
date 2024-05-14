#pragma once
#include "monomial.h"
#include <vector>
#include <queue>

namespace FF4 {
    namespace NUtils {

        template <typename TCoef, typename TComp>
        class Polynomial {
            using TMonomials = std::vector<Monomial<TCoef>>;

            bool isPolynomialCorrect(TMonomials& monomials) {
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
                    for (size_t i = 0; i < monomials.size(); i++) {
                        std::cout << monomials[i] << ' ';
                    }
                    monomials.erase(monomials.end() - cnt, monomials.end());
                    return false;
                }

                if (!isSorted) {
                    std::cout << "!IsSorted" << std::endl;
                    for (size_t i = 0; i < monomials.size(); i++) {
                        std::cout << monomials[i] << ' ';
                    }
                    std::cout << std::endl;
                    std::sort(monomials.rbegin(), monomials.rend(), TComp());
                    return false;
                }
                return true;
            }

        public:
            Polynomial() = default;

            Polynomial(TMonomials&& monomials)
            : monomials_(std::move(monomials))
            {
                assert(isPolynomialCorrect(monomials_));
            }

            const TMonomials& GetMonomials() const noexcept {
                return monomials_;
            }

            const Monomial<TCoef>& GetLeadingMonomial() const noexcept {
                return monomials_[0];
            }

            const Term& GetLeadingTerm() const noexcept {
                return monomials_[0].GetTerm();
            }

            bool IsZero() const noexcept {
                return monomials_.empty();
            }

            Polynomial operator+() const noexcept {
                return *this;
            }

            void Normalize() noexcept {
                const TCoef leadingCoef = monomials_[0].GetCoef();
                if (leadingCoef == 1) {
                    return;
                }
                assert(leadingCoef != 0);
                for (size_t i = 0; i < monomials_.size(); i++) {
                    monomials_[i] /= leadingCoef;
                }
            }

            Polynomial operator-() const noexcept {
                TMonomials monomials = monomials_;
                for (size_t i = 0; i < monomials.size(); i++) {
                    monomials[i] *= -1;
                }
                return Polynomial(std::move(monomials));
            }

            friend bool operator==(const Polynomial& left, const Polynomial& right) noexcept {
                const TMonomials& leftMonomials = left.GetMonomials();
                const TMonomials& rightMonomials = right.GetMonomials();
                if (leftMonomials.size() != rightMonomials.size()) {
                    return false;
                }
                return leftMonomials == rightMonomials;
            }

            friend bool operator!=(const Polynomial& left, const Polynomial& right) noexcept {
                return !(left == right);
            }

            Polynomial& operator+=(const Polynomial& other) noexcept {
                TMonomials monomials;
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

            Polynomial& operator*=(const Term& term) noexcept {
                for (size_t i = 0; i < monomials_.size(); i++) {
                    monomials_[i] *= term;
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

            friend Polynomial operator*(Polynomial left, const Term& right) noexcept {
                left *= right;
                return left;
            }

            friend Polynomial operator*(const Term& left, Polynomial right) noexcept {
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
                const TMonomials& monomials = polynomial.GetMonomials();
                Polynomial newPolynomial;
                for (size_t i = 0; i < monomials.size(); i++) {
                    Polynomial addPolynomial = *this;
                    addPolynomial *= monomials[i];
                    newPolynomial += addPolynomial;
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
                const TMonomials& monomials = polynomial.GetMonomials();
                for (size_t i = 0; i < monomials.size(); i++) {
                    if (i > 0 && monomials[i].GetCoef().IsPositive()) {
                        out << " + ";
                    }
                    out << monomials[i];
                }
                return out;
            }

        private:
            TMonomials monomials_;
        };

        template <typename TCoef, typename TComp>
        using TPolynomials = std::vector<Polynomial<TCoef, TComp>>;

        template <typename TCoef, typename TComp>
        std::ostream& operator<<(std::ostream& out, const TPolynomials<TCoef, TComp>& polynomials) noexcept {
            out << "{\n";
            for (const Polynomial<TCoef, TComp>& polynomial : polynomials) {
                out << polynomial << '\n';
            }
            return out << "}\n";
        }
    }
}
