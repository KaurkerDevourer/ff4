#pragma once

#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std::chrono;

template <typename TFunction>
class TimerWrapper {
public:
    TimerWrapper(TFunction function, std::string text):
        call(std::move(function)), text_(std::move(text)), start_time_(::clock()) {}

    ~TimerWrapper() {
        const clock_t end_time_ = ::clock();
        const clock_t diff = (end_time_ - start_time_);
        std::cout << "| " << text_ << " | ";
        std::cout <<  1000.0 * diff / CLOCKS_PER_SEC << " milliseconds. | " << std::endl;
    }

    TFunction call;

private:
    const std::string text_;
    const clock_t start_time_;
};

template <typename TFunction>
TimerWrapper<TFunction> test_time(TFunction function, const std::string& text) {
    return TimerWrapper<TFunction>(function, text);
}
