#pragma once
#include <vector>
#include <cstdint>
#include <iostream>

namespace FF4 {
    namespace NUtils {
        class Term {
            public:
            using Degree = uint64_t;
                inline Term() = default;

                Term(std::initializer_list<uint64_t>);

                uint64_t& operator[](size_t);
                const uint64_t& operator[](size_t) const;
                const std::vector<uint64_t>& GetData() const;

                void resize(size_t);
                size_t size() const noexcept;
                void reserve(size_t);
                void push_back(uint64_t);

                bool IsOne() const noexcept;
                bool IsDivisibleBy(const Term&) const noexcept;
                Degree TotalDegree() const noexcept;
                std::vector<Term> GetAllDivisors() const noexcept;
                Term& operator*=(const Term&) noexcept;
                friend Term operator*(Term, const Term&) noexcept;
                Term& operator/=(const Term&) noexcept;
                friend Term operator/(Term, const Term&) noexcept;

                friend bool operator<(const Term&, const Term&) noexcept;
                friend bool operator>(const Term&, const Term&) noexcept;
                friend bool operator<=(const Term&, const Term&) noexcept;
                friend bool operator>=(const Term&, const Term&) noexcept;
                friend bool operator==(const Term&, const Term&) noexcept;
                friend bool operator!=(const Term&, const Term&) noexcept;

                friend Term gcd(const Term&, const Term&) noexcept;
                friend Term lcm(const Term&, const Term&) noexcept;

                friend std::ostream& operator<<(std::ostream&, const Term&) noexcept;
            private:
                void Normalize() noexcept;
                friend void FillDivisors(Term&, std::vector<Term>&, size_t) noexcept;
                std::vector<uint64_t> data_;
                Degree sum_ = 0;
        };
    }
}

namespace std {
    template <>
    struct hash<FF4::NUtils::Term> {
        size_t operator()(const FF4::NUtils::Term& t) const {
            const std::vector<uint64_t>& data = t.GetData();
            size_t seed = data.size();
            for (auto x : data) {
                x += 0x9e3779b97f4a7c15; // https://github.com/skeeto/hash-prospector/blob/master/tests/splitmix64.c
                x ^= (x >> 30);
                x *= 0xbf58476d1ce4e5b9;
                x ^= (x >> 27);
                x *= 0x94d049bb133111eb;
                x ^= (x >> 31);
                seed ^= x + 0x9e3779b97f4a7492 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}
