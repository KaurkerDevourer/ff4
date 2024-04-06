#pragma once
#include "monomial.h"

namespace NUtils {
    using TMonomials = std::vector<Monomial>;
    class Polynomial {
    public:
        Polynomial() = default;
        Polynomial(TMonomials&& monomials);

        const TMonomials& GetMonomials() const noexcept;
        const Monomial& GetHeadMonomial() const noexcept;
        bool IsZero() const noexcept;

        Polynomial operator+() const noexcept;
        Polynomial operator-() const noexcept;

        friend bool operator==(const Polynomial&, const Polynomial&) noexcept;
        friend bool operator!=(const Polynomial&, const Polynomial&) noexcept;

        Polynomial& operator+=(const Polynomial&) noexcept;
        friend Polynomial operator+(Polynomial, const Polynomial&) noexcept;

        Polynomial& operator-=(const Polynomial&) noexcept;
        friend Polynomial operator-(Polynomial, const Polynomial&) noexcept;

        Polynomial& operator*=(const Rational&) noexcept;
        friend Polynomial operator*(Polynomial, const Rational&) noexcept;
        friend Polynomial operator*(const Rational&, Polynomial) noexcept;

        Polynomial& operator/=(const Rational&) noexcept;
        friend Polynomial operator/(Polynomial, const Rational&) noexcept;

        Polynomial& operator*=(const Monomial&) noexcept;
        friend Polynomial operator*(Polynomial, const Monomial&) noexcept;
        friend Polynomial operator*(const Monomial&, Polynomial) noexcept;

        Polynomial& operator*=(const Polynomial&) noexcept;
        friend Polynomial operator*(Polynomial, const Polynomial&) noexcept;

        friend std::ostream& operator<<(std::ostream&, const Polynomial&) noexcept;

    private:
        TMonomials&& Normalize(TMonomials&& monomials);

    private:
        TMonomials monomials_;
    };
    using TPolynomials = std::vector<Polynomial>;
}
