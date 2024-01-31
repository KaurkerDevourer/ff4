#include "monomial.h"

Monomial::Monomial(std::vector<int64_t>&& term, Rational coef_)
: term_(std::move(term))
, coef_(coef)
{}


TTerm& Monomial::GetTerm() const noexcept {
    return term_;
}

Rational& Monomial::GetCoef() const noexcept {
    return coef_;
}

