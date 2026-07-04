#include "timer.h"

Timer::Timer() {
    reset();
}

void Timer::reset() {
    start = std::chrono::steady_clock::now();
}

double Timer::time() {
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - start).count();
}
