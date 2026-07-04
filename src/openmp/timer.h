#pragma once

#include <chrono>

class Timer {
public:
    Timer();
    void reset();
    double time();
private:
    std::chrono::time_point<std::chrono::steady_clock> start;
};
