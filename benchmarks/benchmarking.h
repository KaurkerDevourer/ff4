#pragma once
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std::chrono;

using u32 = int32_t;

template <typename TFunction>
class TimerWrapper {
public:
    TimerWrapper(TFunction function, const std::string& text):
        call(function), text_(text), start_time_(high_resolution_clock::now()) {}

    operator TFunction() {
        return call;
    }

    ~TimerWrapper() {
        const time_point<system_clock> end_time_ = high_resolution_clock::now();
        const double diff = (end_time_ - start_time_).count() / 1'000'000;
        std::cout << text_;
        std::cout << diff << " milliseconds." << std::endl;
    }

    TFunction call;

private:
    const std::string text_;
    const time_point<system_clock> start_time_;
};

template <typename TFunction>
TimerWrapper<TFunction> test_time(TFunction function, const std::string& text) {
    return TimerWrapper<TFunction>(function, text);
}

// first argument - function to bench
// second argument - timeout in second
// third argument - args to function
template <typename TFunction, typename ...Args>
void FakeBenchmark(TFunction func, double timeout, std::string name, Args... args) noexcept {
    const auto start = high_resolution_clock::now();
    std::chrono::duration<double> diff;
    bool check = true;
    uint64_t expected = 0;
    uint64_t expected_prev = 0;
    bool is_good = false;
    for (uint64_t i = 1; ; ++i) {
        func(args...);
        if (check || (is_good && i >= expected)) {
            const auto now = high_resolution_clock::now();
            diff = now - start;
            if (diff.count() > timeout) {
                std::cout << name << " Executed " << i << " times in " << diff.count() << " seconds" << std::endl;
                double avg = diff.count() / i;
                avg *= 1'000'000;
                std::cout << "Avg time is " << std::setprecision(9) << avg << " microseconds" << std::endl;
                return;
            }
            double avg = diff.count() / i;
            expected_prev = expected;
            expected = timeout / avg;
            if (double(abs(expected - expected_prev)) / double(expected) < 0.01) {
                is_good = true;
            }
            check = false;
        }
        if (!is_good && i * 2 > expected) {
            check = true;
        }
    }
}
