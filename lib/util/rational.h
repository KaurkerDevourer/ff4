#pragma once

class Rational {
public:
    explicit Rational(int64_t numerator = 0, int64_t denominator = 1);

    int64_t GetNumerator() const noexcept;
    int64_t GetDenominator() const noexcept;

    friend bool operator<(const Rational&, const Rational&) noexcept;
    friend bool operator>(const Rational&, const Rational&) noexcept;
    friend bool operator<=(const Rational&, const Rational&) noexcept;
    friend bool operator>=(const Rational&, const Rational&) noexcept;

    friend bool operator==(const Rational&, const Rational&) noexcept;
    friend bool operator!=(const Rational&, const Rational&) noexcept;

    Rational operator+() const noexcept;
    Rational operator-() const noexcept;

    Rational& operator+=(const Rational&) noexcept;
    friend Rational operator+(Rational, const Rational&) noexcept;

    Rational& operator-=(const Rational&) noexcept;
    friend Rational operator-(Rational, const Rational&) noexcept;

    Rational& operator*=(const Rational&) noexcept;
    friend Rational operator*(Rational, const Rational&) noexcept;

    Rational& operator/=(const Rational&);
    friend Rational operator/(Rational, const Rational&);

    friend Rational pow(const Rational&, const gb::i64&) noexcept;

    friend std::ostream& operator<<(std::ostream&, const Rational&) noexcept;

private:
    void MakeIrreducible() noexcept;

private:
    int64_t numerator_;
    int64_t denominator_;
}