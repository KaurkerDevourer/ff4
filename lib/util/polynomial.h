#pragma once
#include "monomial.h"

using TMonomials = std::vector<Monomial>;
class Polynomial {
public:
    explicit Polynomial(TMonomials&& monomials);

    TMonomials& GetMonomials() const noexcept;
    Monomial& GetHeadMonomial() const noexcept;

    Polynomial operator-() const noexcept;
    Polynomial operator+() const noexcept;

    Polynomial& operator+=(const Polynomial&) noexcept;
    Polynomial& operator-=(const Polynomial&) noexcept;
    Polynomial& operator*=(const Monomial&) noexcept;

private:
    void Normalize(TMonomials& monomials);

private:
    TMonomials monomials_;
}