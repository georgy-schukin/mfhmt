#include "output.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

/*void print(const DArray3 &data, int i_start, int i_end, int j_start, int j_end, ofstream &out) {
    out << setw(7) << "";
    for (int i = i_start; i < i_end; i++) {
        out << setw(3) << i + 1;
        if (i < i_end - 1) {
            out << setw(6) << "";
        }
    }
    out << std::endl;
    for (int j = j_end - 1; j >= j_start; j--) {
        out << setw(3) << j + 1 << setw(1) << "";
        for (int i = i_start; i < i_end; i++) {
            out << setw(10) << std::fixed << setprecision(3) << data(i, j);
        }
        out << std::endl;
    }
}

void output(const string &header, const DArray3 &data, int i_start, int i_end, int j_start, int j_end, ofstream &out) {
    out << "\n";
    out << " " << header << "\n";
    print(data, i_start, i_end, j_start, j_end, out);
}

void output(const string &header, const DArray3 &data, const std::array<int, 4> &range, ofstream &out) {
    output(header, data, range[0], range[1], range[2], range[3], out);
}

void outputFull(const string &header, const DArray3 &data, ofstream &out) {
    output(header, data, {0, data.size(0), 0, data.size(1)}, out);
}*/

std::string formatS(double v, int width, int precision) {
    std::ostringstream out;
    out << std::scientific << std::uppercase << std::setw(width) << std::setprecision(precision) << v;
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
                out << setw(12) << std::scientific << std::uppercase << setprecision(4) << data(i, l , k);
                out << std::endl;
            }
        }
    }
}
