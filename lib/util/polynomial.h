#pragma once
#include "monomial.h"

namespace NUtils {
    class Polynomial {
    public:
        explicit Polynomial(TMonomials&& monomials);

        const TMonomials& GetMonomials() const noexcept;
        const Monomial& GetHeadMonomial() const noexcept;
        bool IsZero() const noexcept;

        Polynomial operator+() const noexcept;
        Polynomial operator-() const noexcept;

        Polynomial& operator+=(const Polynomial&) noexcept;
        Polynomial& operator-=(const Polynomial&) noexcept;
        Polynomial& operator*=(const Monomial&) noexcept;

        friend Polynomial operator*(Polynomial, const Monomial&) noexcept;
        friend Polynomial operator-(Polynomial, const Polynomial&) noexcept;

        bool ReduceBy(const std::vector<Polynomial>& F) noexcept;

    private:
        TMonomials&& Normalize(TMonomials&& monomials);

    private:
        TMonomials monomials_;
    };
    using TPolynomials = std::vector<Polynomial>;
}
