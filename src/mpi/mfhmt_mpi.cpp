#include "defs.h"
#include "block_decomp.h"
#include "distributed_array3d.h"
#include "cart_topology.h"
#include "../common/output.h"
#include "../common/timer.h"

#include <mpi.h>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <tuple>

using namespace std;

using DDArray3 = DistributedArray3D;

void copySliceX(DDArray3 &arr, int dst, int src, int lend, int kend) {
    const auto dst_l = arr.range(0).toLocal(dst);
    const auto src_l = arr.range(0).toLocal(src);
    const auto lend_l = arr.range(1).localEnd(lend);
    const auto kend_l = arr.range(2).localEnd(kend);
    for (int l = 0; l < lend_l; l++) {
        for (int k = 0; k < kend_l; k++) {
            arr(dst_l, l, k) = arr(src_l, l, k);
        }
    }
}

void copySliceY(DDArray3 &arr, int dst, int src, int iend, int kend) {
    const auto dst_l = arr.range(1).toLocal(dst);
    const auto src_l = arr.range(1).toLocal(src);
    const auto iend_l = arr.range(0).localEnd(iend);
    const auto kend_l = arr.range(2).localEnd(kend);
    for (int i = 0; i < iend_l; i++) {
        for (int k = 0; k < kend_l; k++) {
            arr(i, dst_l, k) = arr(i, src_l, k);
        }
    }
}

void copySliceZ(DDArray3 &arr, int dst, int src, int iend, int lend) {
    const auto dst_l = arr.range(2).toLocal(dst);
    const auto src_l = arr.range(2).toLocal(src);
    const auto iend_l = arr.range(0).localEnd(iend);
    const auto lend_l = arr.range(1).localEnd(lend);
    for (int i = 0; i < iend_l; i++) {
        for (int l = 0; l < lend_l; l++) {
            arr(i, l, dst_l) = arr(i, l, src_l);
        }
    }
}

void addSliceZ(DDArray3 &arr, int dst, int src, int iend, int lend) {
    const auto dst_l = arr.range(2).toLocal(dst);
    const auto src_l = arr.range(2).toLocal(src);
    const auto iend_l = arr.range(0).localEnd(iend);
    const auto lend_l = arr.range(1).localEnd(lend);
    for (int i = 0; i < iend_l; i++) {
        for (int l = 0; l < lend_l; l++) {
            arr(i, l, dst_l) += arr(i, l, src_l);
        }
    }
}

void checkAndAdd(DDArray3 &arr, int i, int j, int k, double value) {
    if (arr.hasIndex(i, j, k)) {
        const auto l_ind = arr.toLocal(i, j, k);
        arr(l_ind[0], l_ind[1], l_ind[2]) += value;
    }
}

int main(int argc, char **argv) {

    const int IM_DEF = 40;
    const int LM_DEF = 40;
    const int KM_DEF = 20;
    const int NM_DEF = 1000;
    const int FOUT_DEF = 1;
    const int SOUT_DEF = 0;

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc > 1) {
        const auto s = string(argv[1]);
        if (s == "-h" || s == "--help") {
            std::cout << argv[0] <<
                " [im=" << IM_DEF << "]" <<
                " [lm=" << LM_DEF << "]" <<
                " [km=" << KM_DEF << "]" <<
                " [nm=" << NM_DEF << "]" <<
                " [file_output=" << FOUT_DEF << "]" <<
                " [screen_output=" << SOUT_DEF << "]" <<
                std::endl;
            MPI_Finalize();
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int lm = (argc > 2) ? stoi(argv[2]) : LM_DEF;
    const int km = (argc > 3) ? stoi(argv[3]) : KM_DEF;
    const int nm = (argc > 4) ? stoi(argv[4]) : NM_DEF;
    const bool file_output = (argc > 5) ? stoi(argv[5]) : FOUT_DEF;
    const bool screen_output = (argc > 6) ? stoi(argv[6]) : SOUT_DEF;

    // Can't distribute by Z because of data dependencies, so make 2D lattice only.
    CartTopology<3> cart_tp(rank, size, {0, 0, 1});

    const auto im_decomp = BlockDecomposition(im, cart_tp.dim(0)).grown(1); // decomp im + 2
    const auto lm_decomp = BlockDecomposition(lm, cart_tp.dim(1)).grown(1); // decomp lm + 2
    const auto km_decomp = BlockDecomposition(km, cart_tp.dim(2)).grown(1); // decomp km + 2
    BlockDecomposition3D decomp3d(im_decomp, lm_decomp, km_decomp);

    DDArray3 jx(decomp3d, cart_tp), jy(decomp3d, cart_tp), jz(decomp3d, cart_tp);
    DDArray3 ax(decomp3d, cart_tp), ay(decomp3d, cart_tp), az(decomp3d, cart_tp);
    DDArray3 bx(decomp3d, cart_tp), by(decomp3d, cart_tp), bz(decomp3d, cart_tp);

    const double qj = 10.0;
    const double pi = 3.14159265358979;
    const double c1 = 2.0 * pi / nm;
    const double xm = 4.0;
    const double ym = 4.0;
    const double zm = 2.0;
    const double x0 = 2.0;
    const double y0 = 2.0;
    const double hx = xm / im;
    const double hy = ym / lm;
    const double hz = zm / km;
    const double rhx = 1.0 / hx;
    const double rhy = 1.0 / hy;
    const double rhz = 1.0 / hz;
    const double r0 = 1.0;
    const double h0 = zm / (2.0 * pi);

    std::ofstream out_lst;

    if (file_output) {
        out_lst.open("output.lst");
        out_lst << "qj=" << formatF(qj) << std::endl;
        out_lst << "1. im,lm,km,nm=" << formatI(im) << formatI(lm) << formatI(km) << formatI(nm) << std::endl;
        out_lst << "1. hx,hy,hz=" << formatF(hx) << formatF(hy) << formatF(hz) << std::endl;
    }

    // Init functions

    auto gatherAndOutput = [rank, file_output](const std::string &filename, const DDArray3 &arr, int sx, int sy, int sz) {
        if (file_output) {
            auto gathered = arr.gather(0);
            if (rank == 0) {
                outputDat(filename, gathered, sx, sy, sz);
            }
        }
    };

    auto pqr = [hx, hy, hz, qj](DDArray3 &p, DDArray3 &q, DDArray3 &r, int i, int l, int k, double x, double y, double z, double x1, double y1, double z1) {
        double dx = 0.5 * (x + x1) - hx * (i - 1.5);
        double dy = 0.5 * (y + y1) - hy * (l - 1.5);
        double dz = 0.5 * (z + z1) - hz * (k - 1.5);
        double dx1 = hx - dx;
        double dy1 = hy - dy;
        double dz1 = hz - dz;
        double su = x1 - x;
        double sv = y1 - y;
        double sw = z1 - z;
        double s1 = sv * sw / 12.0;
        double s2 = su * sw / 12.0;
        double s3 = su * sv / 12.0;
        su = su * qj / (hy * hz);
        sv = sv * qj / (hx * hz);
        sw = sw * qj / (hx * hy);
        // Convert to C++ indices.
        i = i - 1;
        l = l - 1;
        k = k - 1;
        checkAndAdd(p, i, l, k, su * (dy1 * dz1 + s1));
        checkAndAdd(p, i, l, k+1, su * (dy1 * dz - s1));
        checkAndAdd(p, i, l+1, k, su * (dy * dz1 - s1));
        p(i,l+1,k+1) += su * (dy * dz + s1);
        q(i,l,k) += sv * (dx1 * dz1 + s2);
        q(i,l,k+1) += sv * (dx1 * dz - s2);
        q(i+1,l,k) += sv * (dx * dz1 - s2);
        q(i+1,l,k+1) += sv * (dx * dz + s2);
        r(i,l,k) += sw * (dx1 * dy1 + s3);
        r(i,l+1,k) += sw * (dx1 * dy - s3);
        r(i+1,l,k) += sw * (dx * dy1 - s3);
        r(i+1,l+1,k) += sw * (dx * dy + s3);
    };

    auto initHelicalCurrent = [x0, y0, r0, h0, nm, hx, hy, hz, rhx, rhy, rhz, c1, &pqr](DDArray3 &jx, DDArray3 &jy, DDArray3 &jz) {
        double x = x0 + r0;
        double y = y0;
        double z = 0.0;
        for (int n = 1; n <= nm; n++) {
            const double x1 = x0 + r0 * std::cos(c1 * n);
            const double y1 = y0 + r0 * std::sin(c1 * n);
            const double z1 = h0 * c1 * n;
            const int i1 = std::trunc(x * rhx + 1.5);
            const int l1 = std::trunc(y * rhy + 1.5);
            const int k1 = std::trunc(z * rhz + 1.5);
            const int i2 = std::trunc(x1 * rhx + 1.5);
            const int l2 = std::trunc(y1 * rhy + 1.5);
            const int k2 = std::trunc(z1 * rhz + 1.5);
            const int i = std::abs(i2 - i1);
            const int l = std::abs(l2 - l1);
            const int k = std::abs(k2 - k1);
            const int m = 4 * i + 2 * l + k;
            if ((m < 1) || (m > 7)) {
                pqr(jx, jy, jz, i1, l1, k1, x, y, z, x1, y1, z1);
            } else {
                double x2, y2, z2;
                double s;
                switch (m) {
                case 1:
                    z2 = hz * (0.5 * (k1 + k2) - 1.0);
                    s = (z2 - z) / (z1 - z);
                    x2 = x + (x1 - x) * s;
                    y2 = y + (y1 - y) * s;
                    break;
                case 2:
                    y2 = hy * (0.5 * (l1 + l2) - 1.0);
                    s = (y2 - y) / (y1 - y);
                    x2 = x + (x1 - x) * s;
                    z2 = z + (z1 - z) * s;
                    break;
                case 3:
                    y2 = hy * (0.5 * (l1 + l2) - 1.0);
                    z2 = hz * (0.5 * (k1 + k2) - 1.0);
                    s = ((z1 - z) * (z2 - z) + (y1 - y) * (y2 - y)) / ((z1 - z) * (z1 - z) + (y1 - y) * (y1 - y));
                    x2 = x + (x1 - x) * s;
                    break;
                case 4:
                    x2 = hx * (0.5 * (i1 + i2) - 1.0);
                    s = (x2 - x) / (x1 - x);
                    y2 = y + (y1 - y) * s;
                    z2 = z + (z1 - z) * s;
                    break;
                case 5:
                    x2 = hx * (0.5 * (i1 + i2) - 1.0);
                    z2 = hz * (0.5 * (k1 + k2) - 1.0);
                    s = ((z1 - z) * (z2 - z) + (x1 - x) * (x2 - x)) / ((z1 - z) * (z1 - z) + (x1 - x) * (x1 - x));
                    y2 = y + (y1 - y) * s;
                    break;
                case 6:
                    x2 = hx * (0.5 * (i1 + i2) - 1.0);
                    y2 = hy * (0.5 * (l1 + l2) - 1.0);
                    s=((y1 - y) * (y2 - y) + (x1 - x) * (x2 - x)) / ((y1 - y) * (y1 - y) + (x1 - x) * (x1 - x));
                    z2 = z + (z1 - z) * s;
                    break;
                case 7:
                    x2 = hx * (0.5 * (i1 + i2) - 1.0);
                    y2 = hy * (0.5 * (l1 + l2) - 1.0);
                    z2 = hz * (0.5 * (k1 + k2) - 1.0);
                    break;
                }
                pqr(jx, jy, jz, i1, l1, k1, x, y, z, x2, y2, z2);
                pqr(jx, jy, jz, i2, l2, k2, x2, y2, z2, x1, y1, z1);
            }
            x = x1;
            y = y1;
            z = z1;
        }
    };

    const double hx2 = hx * hx;
    const double hy2 = hy * hy;
    const double hz2 = hz * hz;
    const double rhx2 = 1.0 / hx2;
    const double rhy2 = 1.0 / hy2;
    const double rhz2 = 1.0 / hz2;
    const double eps = 1e-10;
    const double c2 = 2.0 * rhx2 + 2.0 * rhy2 + 2.0 * rhz2;
    const double rc2 = 1.0 / c2;
    const double c12 = hx / hy;
    const double c13 = hx / hz;
    const double c21 = hy / hx;
    const double c23 = hy / hz;

    auto updateCurrent = [im, lm, km](DDArray3 &jx, DDArray3 &jy, DDArray3 &jz) {
        addSliceZ(jx, 1, km + 1, im + 2, lm + 2); // z: Jx(1) <- Jx(km+1)
        addSliceZ(jx, km, 0, im + 2, lm + 2); // z: Jx(km) <- Jx(0)
        copySliceZ(jx, 0, km, im + 2, lm + 2); // z: Jx(0) <- Jx(km)
        copySliceZ(jx, km + 1, 1, im + 2, lm + 2); // z: Jx(km+1) <- Jx(1)

        addSliceZ(jy, 1, km + 1, im + 2, lm + 2); // z: Jy(1) <- Jy(km)
        addSliceZ(jy, km, 0, im + 2, lm + 2); // z: Jy(km) <- Jy(0)
        copySliceZ(jy, 0, km, im + 2, lm + 2); // z: Jy(0) <- Jy(km)
        copySliceZ(jy, km + 1, 1, im + 2, lm + 2); // z: Jy(km+1) <- Jy(1)

        addSliceZ(jz, km, 0, im + 2, lm + 2); // z: Jz(km) <- Jz(0)
        copySliceZ(jz, 0, km, im + 2, lm + 2); // z: Jz(0) <- Jz(km)
    };

    auto computeBoundaryX = [im, lm, km](DDArray3 &ax, DDArray3 &ay, DDArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        for (int l = 1; l < lm + 1; l++) {
            for (int k = 1; k < km + 1; k++) {
                ax(dst, l, k) = ax(src1, l, k) +
                                c1 * (ay(src2, l, k) - ay(src2, l - 1, k)) +
                                c2 * (az(src2, l, k) - az(src2, l, k - 1));
            }
        }
    };

    auto computeBoundaryY = [im, lm, km](DDArray3 &ax, DDArray3 &ay, DDArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        for (int i = 1; i < im + 1; i++) {
            for (int k = 1; k < km + 1; k++) {
                ay(i, dst, k) = ay(i, src1, k) +
                                c1 * (ax(i, src2, k) - ax(i - 1, src2, k)) +
                                c2 * (az(i, src2, k) - az(i, src2, k - 1));
            }
        }
    };

    auto computeStep = [rhx2, rhy2, rhz2, rc2](DDArray3 &arr, const DDArray3 &j, int iend, int lend, int kend) -> double {
        double maxdiff = 0.0;
        for (int k = 1; k < kend; k++) {
            for (int l = 1; l < lend; l++) {
                for (int i = 1; i < iend; i++) {
                    const double s = ((arr(i+1,l,k) + arr(i-1,l,k)) * rhx2 +
                                      (arr(i,l+1,k) + arr(i,l-1,k)) * rhy2 +
                                      (arr(i,l,k+1) + arr(i,l,k-1)) * rhz2 + j(i,l,k)) * rc2;
                    maxdiff = std::max(std::abs(arr(i,l,k) - s), maxdiff);
                    arr(i,l,k) = s;
                }
            }
        }
        return maxdiff;
    };

    auto computeStepZ = [rhx2, rhy2, rhz2, rc2](DDArray3 &arr, const DDArray3 &j, int iend, int lend, int kend) -> double {
        double maxdiff = 0.0;
        for (int l = 1; l < lend; l++) {
            for (int i = 1; i < iend; i++) {
                for (int k = 1; k < kend; k++) {
                    const double s = ((arr(i+1,l,k) + arr(i-1,l,k)) * rhx2 +
                                      (arr(i,l+1,k) + arr(i,l-1,k)) * rhy2 +
                                      (arr(i,l,k+1) + arr(i,l,k-1)) * rhz2 + j(i,l,k)) * rc2;
                    maxdiff = std::max(std::abs(arr(i,l,k) - s), maxdiff);
                    arr(i,l,k) = s;
                }
                const double s = ((arr(i+1,l,kend) + arr(i-1,l,kend)) * rhx2 +
                                  (arr(i,l+1,kend) + arr(i,l-1,kend)) * rhy2 +
                                  (arr(i,l,1) + arr(i,l,kend-1)) * rhz2 + j(i,l,kend)) * rc2;
                maxdiff = std::max(std::abs(arr(i,l,kend) - s), maxdiff);
                arr(i,l,kend) = s;
                arr(i,l,0) = s;
            }
        }
        return maxdiff;
    };

    auto computeB = [im, lm, km](int dim, const DDArray3 &a1, const DDArray3 &a2, DDArray3 &b, double rh1, double rh2) {
        const static Index3 i_end {2, 1, 1};
        const static Index3 l_end {1, 2, 1};
        const static Index3 k_end {1, 1, 2};
        const static std::array<Index3, 3> shift1 {Index3 {0, 1, 0}, Index3 {0, 0, 1}, Index3 {1, 0, 0}};
        const static std::array<Index3, 3> shift2 {Index3 {0, 0, 1}, Index3 {1, 0, 0}, Index3 {0, 1, 0}};

        const auto &s1 = shift1[dim];
        const auto &s2 = shift2[dim];
        for (int k = 0; k < km + k_end[dim]; k++) {
            for (int l = 0; l < lm + l_end[dim]; l++) {
                for (int i = 0; i < im + i_end[dim]; i++) {
                    b(i,l,k) = (a1(i+s1[0],l+s1[1],k+s1[2]) - a1(i,l,k)) * rh1 -
                                 (a2(i+s2[0],l+s2[1],k+s2[2]) - a2(i,l,k)) * rh2;
                }
            }
        }
    };

    auto outputMax = [&file_output, &out_lst, rank](const std::string &str, double m) {
        if (file_output && rank == 0) {
            out_lst << str << formatS(m) << std::endl;
        }
    };

    auto outputMaxI = [&file_output, &out_lst, rank](const std::string &str, double m, const std::array<int, 3> &ind) {
        if (file_output && rank == 0) {
            out_lst << str << formatS(m) << formatI(ind[0] + 1) << formatI(ind[1] + 1) << formatI(ind[2] + 1) << std::endl;
        }
    };

    auto computeDiv = [im, lm, km, rhx, rhy, rhz](const DDArray3 &x, const DDArray3 &y, const DDArray3 &z, int start, int shift1, int shift2)
        -> std::pair<double, Index3> {
        double maxval = 0.0;
        Index3 maxind {0, 0, 0};
        for (int k = start; k < km + 1; k++) {
            for (int l = start; l < lm + 1; l++) {
                for (int i = start; i < im + 1; i++) {
                    const double s = (x(i + shift1, l, k) - x(i + shift2, l, k)) * rhx +
                                     (y(i, l + shift1, k) - y(i, l + shift2, k)) * rhy +
                                     (z(i, l, k + shift1) - z(i, l, k + shift2)) * rhz;
                    if (std::abs(s) > maxval) {
                        maxval = std::abs(s);
                        maxind = {i, l, k};
                    }
                }
            }
        }
        return std::make_pair(maxval, maxind);
    };

    auto computeRotDiff = [im, lm, km](int dim, const DDArray3 &b1, const DDArray3 &b2, const DDArray3 &j, double rh1, double rh2)
        -> std::pair<double, Index3> {
        const static Index3 i_start {0, 1, 1};
        const static Index3 l_start {1, 0, 1};
        const static Index3 k_start {1, 1, 0};
        const static std::array<Index3, 3> shift1 {Index3 {0, -1, 0}, Index3 {0, 0, -1}, Index3 {-1, 0, 0}};
        const static std::array<Index3, 3> shift2 {Index3 {0, 0, -1}, Index3 {-1, 0, 0}, Index3 {0, -1, 0}};

        double maxval = 0.0;
        Index3 maxind {0, 0, 0};
        const auto &s1 = shift1[dim];
        const auto &s2 = shift2[dim];
        for (int k = k_start[dim]; k < km + 1; k++) {
            for (int l = l_start[dim]; l < lm + 1; l++) {
                for (int i = i_start[dim]; i < im + 1; i++) {
                    const double s = (b1(i,l,k) - b1(i+s1[0],l+s1[1],k+s1[2])) * rh1 -
                                     (b2(i,l,k) - b2(i+s2[0],l+s2[1],k+s2[2])) * rh2 - j(i,l,k);
                    if (std::abs(s) > maxval) {
                        maxval = std::abs(s);
                        maxind = Index3 {i, l, k};
                    }
                }
            }
        }
        return std::make_pair(maxval, maxind);
    };

    // The main program

    Timer work_timer;

    initHelicalCurrent(jx, jy, jz);
    updateCurrent(jx, jy, jz);

    int iter = 0;
    double sx = 0.0, sy = 0.0, sz = 0.0;

    do {
        iter++;

        sx = computeStep(ax, jx, im, lm + 1, km + 1); // Ax(1..im-1,1..lm,1..km) <- Ax(0..im,0..lm+1,0..km+1)

        computeBoundaryX(ax, ay, az, 0, 1, 1, c12, c13); // x: Ax(0) <- Ay(1), Az(1) (local)
        computeBoundaryX(ax, ay, az, im, im - 1, im, -c12, -c13); // x: Ax(im) <- Ay(im-1), Az(im) (local)

        copySliceZ(ax, 0, km, im + 1, lm + 2); // z: Ax(0) <- Ax(km) (remote)
        copySliceZ(ax, km + 1, 1, im + 1, lm + 2); // z: Ax(km+1) <- Ax(1) (remote)

        copySliceY(ax, 0, 1, im + 1, km + 2); // y: Ax(0) <- Ax(1) (local)
        copySliceY(ax, lm + 1, lm, im + 1, km + 2); // y: Ax(lm+1) <- Ax(lm) (local)

        sy = computeStep(ay, jy, im + 1, lm, km + 1); // Ay(1..im,1..lm-1,1..km) <- Ay(0..im+1,0..lm,0..km+1)

        computeBoundaryY(ax, ay, az, 0, 1, 1, c21, c23); // y: Ay(0) <- Ax(1), Ay(1) (local)
        computeBoundaryY(ax, ay, az, lm, lm - 1, lm, -c21, -c23); // y: Ay(lm) <- Ax(lm-1), Az(lm) (local)

        copySliceZ(ay, 0, km, im + 1, lm + 1); // z: Ay(0) <- Ay(km) (remote)
        copySliceZ(ay, km + 1, 1, im + 1, lm + 1); // z: Ay(km+1) <- Ay(1) (remote)

        copySliceX(ay, 0, 1, lm + 1, km + 2); // x: Ay(0) <- Ay(1) (local)
        copySliceX(ay, im + 1, im, lm + 1, km + 2); // x: Ay(im+1) <- Ay(im) (local)

        sz = computeStepZ(az, jz, im + 1, lm + 1, km); // Az(1..im,1..lm,0..km) <- Az(0..im+1,0..lm+1,0..km)

        copySliceZ(az, 0, km, im + 2, lm + 2); // z: Az(0) <- Az(km) (remote)
        copySliceZ(az, km+1, 1, im + 2, lm + 2); // z: Az(km+1) <- Az(1) (remote)

        if (screen_output) {
            std::cout << iter << " " << sx << " " << sy << " " << sz << std::endl;
        }
    } while (sx > eps || sy > eps || sz > eps);

    if (file_output && rank == 0) {
        out_lst << "n,sx,sy,sz=" << formatI(iter) << formatS(sx) << formatS(sy) << formatS(sz) << std::endl;
    }

    computeB(0, az, ay, bx, rhy, rhz);
    computeB(1, ax, az, by, rhz, rhx);
    computeB(2, ay, ax, bz, rhx, rhy);

    const auto mdj = computeDiv(jx, jy, jz, 1, 0, -1);
    const auto mdb = computeDiv(bx, by, bz, 0, 1, 0);
    const auto mda = computeDiv(ax, ay, az, 1, 0, -1);
    outputMax("max(divj)=", mdj.first);
    outputMax("max(divB)=", mdb.first);
    outputMaxI("max(divA)=", mda.first, mda.second);

    const auto mrx = computeRotDiff(0, bz, by, jx, rhy, rhz);
    const auto mry = computeRotDiff(1, bx, bz, jy, rhz, rhx);
    const auto mrz = computeRotDiff(2, by, bx, jz, rhx, rhy);
    outputMaxI("max(rotB_x-jx)=", mrx.first, mrx.second);
    outputMaxI("max(rotB_y-jy)=", mry.first, mry.second);
    outputMaxI("max(rotB_z-jz)=", mrz.first, mrz.second);

    const auto work_time = work_timer.time();

    if (file_output) {
        out_lst.close();
    }

    std::ostringstream out;

    if (rank == 0) {
        out << "Im: " << im << ", Lm: " << lm << ", Km: " << km << ", Nm: " << nm <<
            ", Nodes: " << size <<
            ", Grid: " << cart_tp.dim(0) << "x" << cart_tp.dim(1) << "x" << cart_tp.dim(2) <<
            std::endl;
        out << "TIME: " << work_time << std::endl;
        out << "Iters: " << iter << std::endl;
    }

    std::cout << out.str();

    gatherAndOutput("jx.dat", jx, im + 1, lm + 2, km + 2);
    gatherAndOutput("jy.dat", jy, im + 2, lm + 1, km + 2);
    gatherAndOutput("jz.dat", jz, im + 2, lm + 2, km + 1);
    gatherAndOutput("ax.dat", ax, im + 1, lm + 2, km + 2);
    gatherAndOutput("ay.dat", ay, im + 2, lm + 1, km + 2);
    gatherAndOutput("az.dat", az, im + 2, lm + 2, km + 1);
    gatherAndOutput("bx.dat", bx, im + 2, lm + 1, km + 1);
    gatherAndOutput("by.dat", by, im + 1, lm + 2, km + 1);
    gatherAndOutput("bz.dat", bz, im + 1, lm + 1, km + 2);

    MPI_Finalize();

    return 0;
}
