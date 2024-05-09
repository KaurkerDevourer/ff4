#pragma once
#include <vector>

namespace FF4 {
    namespace NUtils {
        template <typename TCoef>
        class Matrix {
        public:
            Matrix() = delete;

            Matrix(size_t n, size_t m)
            : N_(n)
            , M_(m)
            {
                Data_.resize(N_ * M_);
            }

            TCoef& operator()(size_t i, size_t j) {
                return Data_[i * M_ + j];
            }

            const TCoef& operator()(size_t i, size_t j) const {
                return Data_[i * M_ + j];
            }

            size_t N_;
            size_t M_;
        private:
            std::vector<TCoef> Data_;
        };
    }
}
