#pragma once

#include "common.h"

#include <string>

std::string formatS(double v, int width = 12, int precision = 4);
std::string formatF(double v, int width = 10, int precision = 3);
std::string formatI(int v, int width = 10);

void outputDat(const std::string &filename, const DArray3 &data, int sx, int sy, int sz);
