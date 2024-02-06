#pragma once
#include <vector>
#include "rational.h"

namespace NUtils {
    using TTerm = std::vector<int64_t>;
    class Monomial {
    public:
        Monomial(TTerm term, Rational coef = (int64_t)1);

        Monomial() = default;

        const TTerm& GetTerm() const noexcept;
        const Rational& GetCoef() const noexcept;

        void AddCoef(const Monomial& other) noexcept;
        void SubCoef(const Monomial& other) noexcept;

        bool IsDivisibleBy(const Monomial&) const noexcept __attribute__((pure));

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

        friend Monomial gcd(const Monomial&, const Monomial&) noexcept;

        friend Monomial lcm(const Monomial&, const Monomial&) noexcept;

        static void Normalize(TTerm&);

    private:
        TTerm term_;
        Rational coef_;
    };
    using TMonomials = std::vector<Monomial>;
}
