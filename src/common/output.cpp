#include "output.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;

std::string formatS(double v, int width, int precision) {
    std::ostringstream out;
    out << std::scientific << std::uppercase << std::setw(width) << std::setprecision(precision) << v;
    return out.str();
}

std::string formatF(double v, int width, int precision) {
    std::ostringstream out;
    out << std::fixed << std::uppercase << std::setw(width) << std::setprecision(precision) << v;
    return out.str();
}

std::string formatI(int v, int width) {
    std::ostringstream out;
    out << std::setw(width) << v;
    return out.str();
}

