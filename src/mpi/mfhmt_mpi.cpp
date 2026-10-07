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
#include <map>

using namespace std;

using DDArray3 = DistributedArray3D;
using Double3 = std::array<double, 3>;

void copySliceX(DDArray3 &arr, int dst, int src, int lend, int kend) {
    if (!arr.range(0).hasIndex(dst)) {
        return;
    }
    const auto dst_l = arr.range(0).toLocal(dst);
    const auto src_l = arr.range(0).toLocal(src);
    const auto l_end_l = arr.range(1).localEnd(lend);
    const auto k_end_l = arr.range(2).localEnd(kend);
    for (int l = 0; l < l_end_l; l++) {
        for (int k = 0; k < k_end_l; k++) {
            arr(dst_l, l, k) = arr(src_l, l, k);
        }
    }
}

void copySliceY(DDArray3 &arr, int dst, int src, int iend, int kend) {
    if (!arr.range(1).hasIndex(dst)) {
        return;
    }
    const auto dst_l = arr.range(1).toLocal(dst);
    const auto src_l = arr.range(1).toLocal(src);
    const auto i_end_l = arr.range(0).localEnd(iend);
    const auto k_end_l = arr.range(2).localEnd(kend);
    for (int i = 0; i < i_end_l; i++) {
        for (int k = 0; k < k_end_l; k++) {
            arr(i, dst_l, k) = arr(i, src_l, k);
        }
    }
}

void copySliceZ(DDArray3 &arr, int dst, int src, int iend, int lend) {
    if (!arr.range(2).hasIndex(dst)) {
        return;
    }
    const auto dst_l = arr.range(2).toLocal(dst);
    const auto src_l = arr.range(2).toLocal(src);
    const auto i_end_l = arr.range(0).localEnd(iend);
    const auto l_end_l = arr.range(1).localEnd(lend);
    for (int i = 0; i < i_end_l; i++) {
        for (int l = 0; l < l_end_l; l++) {
            arr(i, l, dst_l) = arr(i, l, src_l);
        }
    }
}

void addSliceZ(DDArray3 &arr, int dst, int src, int iend, int lend) {
    if (!arr.range(2).hasIndex(dst)) {
        return;
    }
    const auto dst_l = arr.range(2).toLocal(dst);
    const auto src_l = arr.range(2).toLocal(src);
    const auto i_end_l = arr.range(0).localEnd(iend);
    const auto l_end_l = arr.range(1).localEnd(lend);
    for (int i = 0; i < i_end_l; i++) {
        for (int l = 0; l < l_end_l; l++) {
            arr(i, l, dst_l) += arr(i, l, src_l);
        }
    }
}

void checkAndAdd(DDArray3 &arr, int i, int j, int k, double value) {
    if (arr.hasIndex(i, j, k)) {
        const auto l_ind = arr.toLocal(i, j, k);
        arr(l_ind) += value;
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

    double comp_time = 0;
    double shadow_time = 0;
    double reduce_time = 0;

    enum ShadowSyncType {
        SHADOW_ALL = 0,
        SHADOW_PREV,
        SHADOW_NEXT
    };

    auto recvShadowsPrev = [&shadow_time](DDArray3 &arr) {
        Timer tm;
        auto ops = arr.recvShadowsPrev();
        ops.wait();
        shadow_time += tm.time();
    };

    auto recvShadowsNext = [&shadow_time](DDArray3 &arr) {
        Timer tm;
        auto ops = arr.recvShadowsNext();
        ops.wait();
        shadow_time += tm.time();
    };

    auto sendShadowsPrevAsync = [&shadow_time](DDArray3 &arr) -> AsyncOps {
        Timer tm;
        auto ops = arr.sendSnadowsPrev();
        shadow_time += tm.time();
        return ops;
    };

    auto sendShadowsNextAsync = [&shadow_time](DDArray3 &arr) -> AsyncOps {
        Timer tm;
        auto ops = arr.sendSnadowsNext();
        shadow_time += tm.time();
        return ops;
    };

    auto finishShadowSends = [&shadow_time](DDArray3 &arr) {
        Timer tm;
        arr.finishAllOps();
        shadow_time += tm.time();
    };

    auto syncShadows = [&shadow_time](DDArray3 &arr, ShadowSyncType type = SHADOW_ALL, int dim = -1) {
        Timer tm;
        AsyncOps ops;
        switch (type) {
            case SHADOW_ALL: ops = (dim >= 0 ? arr.syncShadows(dim) : arr.syncShadows()); break;
            case SHADOW_PREV: ops = (dim >= 0 ? arr.syncShadowsPrev(dim) : arr.syncShadowsPrev()); break;
            case SHADOW_NEXT: ops = (dim >= 0 ? arr.syncShadowsNext(dim) : arr.syncShadowsNext()); break;
        }
        ops.wait();
        shadow_time += tm.time();
    };

    auto syncShadowsI = [&syncShadows](DDArray3 &arr, const Index3 &shifts) {
        for (int dim = 0; dim < 3; dim++) {
            if (shifts[dim] > 0) {
                syncShadows(arr, SHADOW_NEXT, dim);
            } else if (shifts[dim] < 0) {
                syncShadows(arr, SHADOW_PREV, dim);
            }
        }
    };

    auto syncShadowsD = [&syncShadows](DDArray3 &arr, int dim, int shift) {
        if (shift > 0) {
            syncShadows(arr, SHADOW_NEXT, dim);
        } else if (shift < 0) {
            syncShadows(arr, SHADOW_PREV, dim);
        }
    };

    auto reduceMax = [&reduce_time](const Double3 &value) -> Double3 {
        Double3 global;
        Timer tm;
        MPI_Allreduce(value.data(), global.data(), 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        reduce_time += tm.time();
        return global;
    };

    auto reduceMaxInd = [rank, &reduce_time](double value, const Index3 &ind) -> std::pair<double, Index3> {
        struct {
            double value;
            int rank;
        } local, global;
        local.value = value;
        local.rank = rank;
        Index3 max_index = ind;
        Timer tm;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        MPI_Bcast(max_index.data(), 3, MPI_INT, global.rank, MPI_COMM_WORLD);
        reduce_time += tm.time();
        return std::make_pair(global.value, max_index);
    };

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
        checkAndAdd(p, i, l+1, k+1, su * (dy * dz + s1));
        checkAndAdd(q, i, l, k, sv * (dx1 * dz1 + s2));
        checkAndAdd(q, i, l, k+1, sv * (dx1 * dz - s2));
        checkAndAdd(q, i+1, l, k, sv * (dx * dz1 - s2));
        checkAndAdd(q, i+1, l, k+1, sv * (dx * dz + s2));
        checkAndAdd(r, i, l, k, sw * (dx1 * dy1 + s3));
        checkAndAdd(r, i, l+1, k, sw * (dx1 * dy - s3));
        checkAndAdd(r, i+1, l, k, sw * (dx * dy1 - s3));
        checkAndAdd(r, i+1, l+1, k, sw * (dx * dy + s3));
    };

    auto initHelicalCurrent = [x0, y0, r0, h0, nm, hx, hy, hz, rhx, rhy, rhz, c1, &pqr, &comp_time](DDArray3 &jx, DDArray3 &jy, DDArray3 &jz) {
        double x = x0 + r0;
        double y = y0;
        double z = 0.0;
        Timer tm;
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
        comp_time += tm.time();
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

    auto updateCurrent = [im, lm, km, &comp_time](DDArray3 &jx, DDArray3 &jy, DDArray3 &jz) {
        Timer tm;
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
        comp_time += tm.time();
    };

    auto computeBoundaryX = [im, lm, km, &comp_time](DDArray3 &ax, DDArray3 &ay, DDArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        if (!ax.range(0).hasIndex(dst)) {
            return;
        }
        const auto dst_l = ax.range(0).toLocal(dst);
        const auto src1_l = ax.range(0).toLocal(src1);
        const auto src2_l = ax.range(0).toLocal(src2);
        const auto l_start_l = ax.range(1).localStart(1);
        const auto l_end_l = ax.range(1).localEnd(lm + 1);
        Timer tm;
        for (int l = l_start_l; l < l_end_l; l++) {
            for (int k = 1; k < km + 1; k++) {
                ax(dst_l, l, k) = ax(src1_l, l, k) +
                                c1 * (ay(src2_l, l, k) - ay(src2_l, l - 1, k)) +
                                c2 * (az(src2_l, l, k) - az(src2_l, l, k - 1));
            }
        }
        comp_time += tm.time();
    };

    auto computeBoundaryY = [im, lm, km, &comp_time](DDArray3 &ax, DDArray3 &ay, DDArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        if (!ay.range(1).hasIndex(dst)) {
            return;
        }
        const auto dst_l = ay.range(1).toLocal(dst);
        const auto src1_l = ay.range(1).toLocal(src1);
        const auto src2_l = ay.range(1).toLocal(src2);
        const auto i_start_l = ax.range(0).localStart(1);
        const auto i_end_l = ax.range(0).localEnd(im + 1);
        Timer tm;
        for (int i = i_start_l; i < i_end_l; i++) {
            for (int k = 1; k < km + 1; k++) {
                ay(i, dst_l, k) = ay(i, src1_l, k) +
                                c1 * (ax(i, src2_l, k) - ax(i - 1, src2_l, k)) +
                                c2 * (az(i, src2_l, k) - az(i, src2_l, k - 1));
            }
        }
        comp_time += tm.time();
    };

    auto computeStep = [rhx2, rhy2, rhz2, rc2, &recvShadowsPrev, &recvShadowsNext, &comp_time](DDArray3 &arr, const DDArray3 &j, int iend, int lend, int kend) -> double {
        recvShadowsNext(arr);
        recvShadowsPrev(arr);
        arr.finishAllOps();
        double maxdiff = 0.0;
        const auto i_start_l = arr.range(0).localStart(1);
        const auto i_endl = arr.range(0).localEnd(iend);
        const auto l_start_l = arr.range(1).localStart(1);
        const auto l_end_l = arr.range(1).localEnd(lend);
        Timer tm;
        for (int i = i_start_l; i < i_endl; i++) {
            for (int l = l_start_l; l < l_end_l; l++) {
                for (int k = 1; k < kend; k++) {
                    const double s = ((arr(i+1,l,k) + arr(i-1,l,k)) * rhx2 +
                                      (arr(i,l+1,k) + arr(i,l-1,k)) * rhy2 +
                                      (arr(i,l,k+1) + arr(i,l,k-1)) * rhz2 + j(i,l,k)) * rc2;
                    maxdiff = std::max(std::abs(arr(i,l,k) - s), maxdiff);
                    arr(i,l,k) = s;
                }
            }
        }
        comp_time += tm.time();
        return maxdiff;
    };

    auto computeStepZ = [rhx2, rhy2, rhz2, rc2, &recvShadowsPrev, &recvShadowsNext, &comp_time](DDArray3 &arr, const DDArray3 &j, int iend, int lend, int kend) -> double {
        recvShadowsNext(arr);
        recvShadowsPrev(arr);
        arr.finishAllOps();
        double maxdiff = 0.0;
        const auto i_start_l = arr.range(0).localStart(1);
        const auto i_end_l = arr.range(0).localEnd(iend);
        const auto l_start_l = arr.range(1).localStart(1);
        const auto l_end_l = arr.range(1).localEnd(lend);
        Timer tm;
        for (int i = i_start_l; i < i_end_l; i++) {
            for (int l = l_start_l; l < l_end_l; l++) {
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
        comp_time += tm.time();
        return maxdiff;
    };

    auto computeIterationStep = [im, lm, km, c12, c13, c21, c23, &computeStep, &computeStepZ, &computeBoundaryX, &computeBoundaryY,
                                 &sendShadowsPrevAsync, sendShadowsNextAsync]
        (DDArray3 &ax, DDArray3 &ay, DDArray3 &az, const DDArray3 &jx, const DDArray3 &jy, const DDArray3 &jz)
        -> Double3 {
        const auto sx = computeStep(ax, jx, im, lm + 1, km + 1); // Ax(1..im-1,1..lm,1..km) <- Ax(0..im,0..lm+1,0..km+1)

        computeBoundaryX(ax, ay, az, 0, 1, 1, c12, c13); // x: Ax(0) <- Ay(1), Az(1) (local)
        computeBoundaryX(ax, ay, az, im, im - 1, im, -c12, -c13); // x: Ax(im) <- Ay(im-1), Az(im) (local)

        copySliceZ(ax, 0, km, im + 1, lm + 2); // z: Ax(0) <- Ax(km) (remote)
        copySliceZ(ax, km + 1, 1, im + 1, lm + 2); // z: Ax(km+1) <- Ax(1) (remote)

        copySliceY(ax, 0, 1, im + 1, km + 2); // y: Ax(0) <- Ax(1) (local)
        copySliceY(ax, lm + 1, lm, im + 1, km + 2); // y: Ax(lm+1) <- Ax(lm) (local)

        sendShadowsNextAsync(ax);
        sendShadowsPrevAsync(ax);

        const auto sy = computeStep(ay, jy, im + 1, lm, km + 1); // Ay(1..im,1..lm-1,1..km) <- Ay(0..im+1,0..lm,0..km+1)

        computeBoundaryY(ax, ay, az, 0, 1, 1, c21, c23); // y: Ay(0) <- Ax(1), Ay(1) (local)
        computeBoundaryY(ax, ay, az, lm, lm - 1, lm, -c21, -c23); // y: Ay(lm) <- Ax(lm-1), Az(lm) (local)

        copySliceZ(ay, 0, km, im + 1, lm + 1); // z: Ay(0) <- Ay(km) (remote)
        copySliceZ(ay, km + 1, 1, im + 1, lm + 1); // z: Ay(km+1) <- Ay(1) (remote)

        copySliceX(ay, 0, 1, lm + 1, km + 2); // x: Ay(0) <- Ay(1) (local)
        copySliceX(ay, im + 1, im, lm + 1, km + 2); // x: Ay(im+1) <- Ay(im) (local)

        sendShadowsNextAsync(ay);
        sendShadowsPrevAsync(ay);

        const auto sz = computeStepZ(az, jz, im + 1, lm + 1, km); // Az(1..im,1..lm,0..km) <- Az(0..im+1,0..lm+1,0..km)

        copySliceZ(az, 0, km, im + 2, lm + 2); // z: Az(0) <- Az(km) (remote)
        copySliceZ(az, km+1, 1, im + 2, lm + 2); // z: Az(km+1) <- Az(1) (remote)

        sendShadowsNextAsync(az);
        sendShadowsPrevAsync(az);

        return Double3 {sx, sy, sz};
    };

    auto computeIterationStepEmpty = [&recvShadowsPrev, &recvShadowsNext, &sendShadowsPrevAsync, sendShadowsNextAsync]
        (DDArray3 &ax, DDArray3 &ay, DDArray3 &az) {
        recvShadowsNext(ax);
        recvShadowsPrev(ax);
        recvShadowsNext(ay);
        recvShadowsPrev(ay);
        recvShadowsNext(az);
        recvShadowsPrev(az);
    };

    auto computeA = [screen_output, eps, &cart_tp, &computeIterationStep, &computeIterationStepEmpty, &sendShadowsNextAsync, &finishShadowSends, &reduceMax]
        (DDArray3 &ax, DDArray3 &ay, DDArray3 &az, const DDArray3 &jx, const DDArray3 &jy, const DDArray3 &jz) -> std::pair<int, Double3> {
        const auto max_iter_depth = cart_tp.dim(0) + cart_tp.dim(1) + cart_tp.dim(2) - 2;
        const auto my_coord = cart_tp.thisCoord();
        const auto my_rank = cart_tp.thisRank();
        const auto my_iter_depth = max_iter_depth - (my_coord[0] + my_coord[1] + my_coord[2]);

        const size_t buf_size = my_iter_depth;
        std::vector<DDArray3> ax_buf, ay_buf, az_buf;
        std::vector<Double3> diff_buf(buf_size);
        for (int i = 0; i < buf_size; i++) {
            ax_buf.push_back(ax);
            ay_buf.push_back(ay);
            az_buf.push_back(az);
        }

        int iter = 0, done_iter = 0;
        Double3 max_diff;
        size_t curr_buf_pos = 0;
        bool do_iter = true;
        int finishing = max_iter_depth - my_iter_depth;

        sendShadowsNextAsync(ax_buf[curr_buf_pos]);
        sendShadowsNextAsync(ay_buf[curr_buf_pos]);
        sendShadowsNextAsync(az_buf[curr_buf_pos]);

        do {
            auto &ax_curr = ax_buf[curr_buf_pos];
            auto &ay_curr = ay_buf[curr_buf_pos];
            auto &az_curr = az_buf[curr_buf_pos];

            done_iter++;
            const auto diff = computeIterationStep(ax_curr, ay_curr, az_curr, jx, jy, jz);
            diff_buf[curr_buf_pos] = diff;

            //std::cout << my_rank << ": done " << done_iter << " " << diff[0] << " " << diff[1] << " " << diff[2] << std::endl;

            if (done_iter >= my_iter_depth) {
                iter++;
                const size_t reduce_pos = (iter - 1) % buf_size;
                max_diff = reduceMax(diff_buf[reduce_pos]);
                do_iter = (max_diff[0] > eps || max_diff[1] > eps || max_diff[2] > eps);

                if (screen_output && my_rank == 0) {
                    std::cout << iter << " " << max_diff[0] << " " << max_diff[1] << " " << max_diff[2] << std::endl;
                }
            }

            if (do_iter) {
                const size_t next_buf_pos = (curr_buf_pos + 1) % buf_size;
                finishShadowSends(ax_buf[next_buf_pos]);
                finishShadowSends(ay_buf[next_buf_pos]);
                finishShadowSends(az_buf[next_buf_pos]);
                if (next_buf_pos != curr_buf_pos) {
                    ax_buf[next_buf_pos].copy(ax_curr);
                    ay_buf[next_buf_pos].copy(ay_curr);
                    az_buf[next_buf_pos].copy(az_curr);
                }
                curr_buf_pos = next_buf_pos;
            }
        } while (do_iter);

        while (finishing > 0) {
            computeIterationStepEmpty(ax, ay, az);
            finishing--;
        }

        const size_t solution_pos = (iter - 1) % buf_size;
        ax.copy(ax_buf[solution_pos]);
        ay.copy(ay_buf[solution_pos]);
        az.copy(az_buf[solution_pos]);

        return std::make_pair(iter, max_diff);
    };

    auto computeB = [im, lm, km, &syncShadowsI, &comp_time](int dim, DDArray3 &a1, DDArray3 &a2, DDArray3 &b, double rh1, double rh2) {
        const static Index3 i_end {2, 1, 1};
        const static Index3 l_end {1, 2, 1};
        const static Index3 k_end {1, 1, 2};
        const static std::array<Index3, 3> shift1 {Index3 {0, 1, 0}, Index3 {0, 0, 1}, Index3 {1, 0, 0}};
        const static std::array<Index3, 3> shift2 {Index3 {0, 0, 1}, Index3 {1, 0, 0}, Index3 {0, 1, 0}};

        const auto &s1 = shift1[dim];
        const auto &s2 = shift2[dim];

        syncShadowsI(a1, s1);
        syncShadowsI(a2, s2);

        const auto i_end_l = b.range(0).localEnd(im + i_end[dim]);
        const auto l_end_l = b.range(1).localEnd(lm + l_end[dim]);
        const auto k_end_l = b.range(2).localEnd(km + k_end[dim]);
        Timer tm;
        for (int i = 0; i < i_end_l; i++) {
            for (int l = 0; l < l_end_l; l++) {
                for (int k = 0; k < k_end_l; k++) {
                    b(i,l,k) = (a1(i+s1[0],l+s1[1],k+s1[2]) - a1(i,l,k)) * rh1 -
                                 (a2(i+s2[0],l+s2[1],k+s2[2]) - a2(i,l,k)) * rh2;
                }
            }
        }
        comp_time += tm.time();
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

    auto computeDiv = [im, lm, km, rhx, rhy, rhz, &syncShadowsD, &reduceMaxInd, &comp_time](DDArray3 &x, DDArray3 &y, DDArray3 &z, int start, int shift1, int shift2)
        -> std::pair<double, Index3> {
        syncShadowsD(x, 0, shift1);
        syncShadowsD(x, 0, shift2);
        syncShadowsD(y, 1, shift1);
        syncShadowsD(y, 1, shift2);
        syncShadowsD(z, 2, shift1);
        syncShadowsD(z, 2, shift2);

        double maxval = 0.0;
        Index3 maxind {0, 0, 0};

        const auto i_start_l = x.range(0).localStart(start);
        const auto l_start_l = x.range(1).localStart(start);
        const auto k_start_l = x.range(2).localStart(start);
        const auto i_end_l = x.range(0).localEnd(im + 1);
        const auto l_end_l = x.range(1).localEnd(lm + 1);
        const auto k_end_l = x.range(2).localEnd(km + 1);
        Timer tm;
        for (int i = i_start_l; i < i_end_l; i++) {
            for (int l = l_start_l; l < l_end_l; l++) {
                for (int k = k_start_l; k < k_end_l; k++) {
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
        comp_time += tm.time();
        return reduceMaxInd(maxval, maxind);
    };

    auto computeRotDiff = [im, lm, km, &syncShadowsI, &reduceMaxInd, &comp_time](int dim, DDArray3 &b1, DDArray3 &b2, const DDArray3 &j, double rh1, double rh2)
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

        syncShadowsI(b1, s1);
        syncShadowsI(b2, s2);

        const auto i_start_l = b1.range(0).localStart(i_start[dim]);
        const auto l_start_l = b1.range(1).localStart(l_start[dim]);
        const auto k_start_l = b1.range(2).localStart(k_start[dim]);
        const auto i_end_l = b1.range(0).localEnd(im + 1);
        const auto l_end_l = b1.range(1).localEnd(lm + 1);
        const auto k_end_l = b1.range(2).localEnd(km + 1);
        Timer tm;
        for (int i = i_start_l; i < i_end_l; i++) {
            for (int l = l_start_l; l < l_end_l; l++) {
                for (int k = k_start_l; k < k_end_l; k++) {
                    const double s = (b1(i,l,k) - b1(i+s1[0],l+s1[1],k+s1[2])) * rh1 -
                                     (b2(i,l,k) - b2(i+s2[0],l+s2[1],k+s2[2])) * rh2 - j(i,l,k);
                    if (std::abs(s) > maxval) {
                        maxval = std::abs(s);
                        maxind = Index3 {i, l, k};
                    }
                }
            }
        }
        comp_time += tm.time();
        return reduceMaxInd(maxval, maxind);
    };

    // The main program

    Timer work_timer;

    initHelicalCurrent(jx, jy, jz);
    updateCurrent(jx, jy, jz);

    int iter = 0;
    Double3 max_diff;
    std::tie(iter, max_diff) = computeA(ax, ay, az, jx, jy, jz);

    if (file_output && rank == 0) {
        out_lst << "n,sx,sy,sz=" << formatI(iter) << formatS(max_diff[0]) << formatS(max_diff[1]) << formatS(max_diff[2]) << std::endl;
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

    out << "Node " << rank <<
        ": Comp time: " << comp_time <<
        ", Shadow time: " << shadow_time <<
        ", Reduce time: " << reduce_time << std::endl;

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

    /*BlockDecomposition3D dd(8, cart_tp.dim(0), 8, cart_tp.dim(1), 1, 1);
    DDArray3 tmp(dd, cart_tp), tmp2(dd, cart_tp);
    for (int i = 0; i < tmp.range(0).localEnd(8);i++)
        for (int j = 0; j < tmp.range(1).localEnd(8);j++)
        {
            const auto ii = tmp.range(0).toGlobal(i);
            const auto jj = tmp.range(1).toGlobal(j);
            tmp(i,j,0) = ii + jj;
        }

    tmp.sendSnadowsPrev();
    tmp.recvShadowsPrev();
    tmp.finishAllOps();

    for (int i = tmp.range(0).localStart(1); i < tmp.range(0).localEnd(8);i++)
        for (int j = tmp.range(1).localStart(1); j < tmp.range(1).localEnd(8);j++)
        {
            tmp2(i,j,0) = tmp(i,j,0) + tmp(i-1,j,0) + tmp(i,j-1,0);
        }

    gatherAndOutput("tmp.dat", tmp, 8, 8, 1);
    gatherAndOutput("tmp2.dat", tmp2, 8, 8, 1);*/

    MPI_Finalize();

    return 0;
}
