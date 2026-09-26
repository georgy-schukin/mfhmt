#pragma once

#include <string>
#include <fstream>
#include <iomanip>

std::string formatS(double v, int width = 12, int precision = 4);
std::string formatF(double v, int width = 10, int precision = 3);
std::string formatI(int v, int width = 10);

template <typename Array3D>
void outputDat(const std::string &filename, const Array3D &data, int sx, int sy, int sz) {
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
                out << std::setw(5) << i + 1 << std::setw(5) << l + 1 << std::setw(5) << k + 1;
                out << std::setw(12) << std::scientific << std::uppercase << std::setprecision(4) << data(i, l, k);
                out << std::endl;
            }
        }
    }
}

