#include "buchberger.cpp"
#include "monomial.cpp"
#include "polynomial.cpp"
#include "rational.cpp"
#include "term.cpp"

int main() {
    test_buchberger();
    test_monomial();
    test_polynomial();
    test_rational();
    test_term();
}
