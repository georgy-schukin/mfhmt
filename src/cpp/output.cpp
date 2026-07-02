#include "output.h"

#include <iostream>
#include <fstream>
#include <iomanip>

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

void outputDat(const std::string &filename, const DArray3 &data, int sx, int sy, int sz) {
    std::ofstream out(filename.c_str());
/*
    do k=1,km+2
        do l=1,lm+2
            do i=1,im+1
               write(17,101) i,l,k,jx(i,l,k)
    101          format(3i5,e12.4)
            enddo
        enddo
    enddo
*/
    for (int k = 0; k < sz; k++) {
        for (int l = 0; l < sy; l++) {
            for (int i = 0; i < sx; i++) {
                out << setw(5) << i + 1 << setw(5) << l + 1 << setw(5) << k + 1;
                out << setw(12) << std::scientific << std::uppercase << setprecision(4) << data(i, l, k);
                out << std::endl;
            }
        }
    }
}
