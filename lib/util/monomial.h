#pragma once
#include "term.h"
#include "rational.h"

namespace NUtils {
    class Monomial {
    public:
        Monomial() = default;
        Monomial(TTerm term, Rational coef = (int64_t)1);

        const TTerm& GetTerm() const noexcept;
        const Rational& GetCoef() const noexcept;
        int64_t GetNumerator() const noexcept;
        int64_t GetDenominator() const noexcept;

        void AddCoef(const Monomial& other) noexcept;
        void SubCoef(const Monomial& other) noexcept;

        friend bool operator<(const Monomial&, const Monomial&) noexcept;
        friend bool operator>(const Monomial&, const Monomial&) noexcept;
        friend bool operator<=(const Monomial&, const Monomial&) noexcept;
        friend bool operator>=(const Monomial&, const Monomial&) noexcept;

        friend bool operator==(const Monomial&, const Monomial&) noexcept;
        friend bool operator!=(const Monomial&, const Monomial&) noexcept;

        Monomial operator+() const noexcept;
        Monomial operator-() const noexcept;

        Monomial& operator*=(const Monomial&) noexcept;
        friend Monomial operator*(Monomial, const Monomial&) noexcept;

        Monomial& operator*=(const Rational&) noexcept;
        friend Monomial operator*(Monomial, const Rational&) noexcept;

        Monomial& operator/=(const Monomial&);
        friend Monomial operator/(Monomial, const Monomial&);

        Monomial& operator/=(const Rational&);
        friend Monomial operator/(Monomial, const Rational&);

        friend std::ostream& operator<<(std::ostream&, const Monomial&) noexcept;

    private:
        TTerm term_;
        Rational coef_;
    };
}
