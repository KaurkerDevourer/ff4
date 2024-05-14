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
                size_t size() const;
                void reserve(size_t);
                void push_back(uint64_t);

                bool IsOne() const noexcept;
                bool IsDivisibleBy(const Term&) const noexcept;
                Degree TotalDegree() const noexcept;
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
                void Normalize();
                std::vector<uint64_t> data_;
                Degree sum_ = 0;
        };
    }
}

namespace std {
    template <>
    struct hash<FF4::NUtils::Term> {
        size_t operator()(const FF4::NUtils::Term& t) const {
            const std::vector<uint64_t>& data = t.GetData(); // https://stackoverflow.com/a/27216842
            size_t seed = data.size();
            for (auto x : data) {
                x = ((x >> 32) ^ x) * 0x45d9f3b;
                x = ((x >> 32) ^ x) * 0x45d9f3b;
                x = (x >> 32) ^ x;
                seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}
