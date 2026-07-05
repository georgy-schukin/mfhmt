#include "common.h"
#include "output.h"
#include "timer.h"

#include <omp.h>

#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <algorithm>
#include <tuple>
#include <functional>

using namespace std;

int main(int argc, char **argv) {
    const int NTHREADS_DEF = 1;
    const int IM_DEF = 40;
    const int LM_DEF = 40;
    const int KM_DEF = 20;
    const int NM_DEF = 1000;
    const int FOUT_DEF = 1;
    const int SOUT_DEF = 0;
    const int TASKS_PER_DIM_DEF = 10;
    const int MIN_TASK_SIZE_DEF = 10;

    if (argc > 1) {
        const auto s = string(argv[1]);
        if (s == "-h" || s == "--help") {
            std::cout << argv[0] <<
                " [nthreads=" << NTHREADS_DEF << "]" <<
                " [im=" << IM_DEF << "]" <<
                " [lm=" << LM_DEF << "]" <<
                " [km=" << KM_DEF << "]" <<
                " [nm=" << NM_DEF << "]" <<
                " [file_output=" << FOUT_DEF << "]" <<
                " [screen_output=" << SOUT_DEF << "]" <<
                " [tasks_per_dim=" << TASKS_PER_DIM_DEF << "]" <<
                " [min_task_size=" << MIN_TASK_SIZE_DEF << "]" <<
                std::endl;
            return 0;
        }
    }

    auto intArg = [&argc,&argv](int index, int dft) {
        return (argc > index) ? stoi(argv[index]) : dft;
    };

    const int num_of_threads = intArg(1, NTHREADS_DEF);
    const int im = intArg(2, IM_DEF);
    const int lm = intArg(3, LM_DEF);
    const int km = intArg(4, KM_DEF);
    const int nm = intArg(5, NM_DEF);
    const bool file_output = intArg(6, FOUT_DEF);
    const bool screen_output = intArg(7, SOUT_DEF);
    const int tasks_per_dim = intArg(8, TASKS_PER_DIM_DEF);
    const int min_task_size = intArg(9, MIN_TASK_SIZE_DEF);

    omp_set_num_threads(num_of_threads);

    const size_t ims = im + 2;
    const size_t lms = lm + 2;
    const size_t kms = km + 2;

    DArray3 jx(ims, lms, kms), jy(ims, lms, kms), jz(ims, lms, kms);
    DArray3 ax(ims, lms, kms), ay(ims, lms, kms), az(ims, lms, kms);
    DArray3 bx(ims, lms, kms), by(ims, lms, kms), bz(ims, lms, kms);

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
        out_lst.open("smt45.lst");
        out_lst << "qj=" << formatF(qj) << std::endl;
        out_lst << "1. im,lm,km,nm=" << formatI(im) << formatI(lm) << formatI(km) << formatI(nm) << std::endl;
        out_lst << "1. hx,hy,hz=" << formatF(hx) << formatF(hy) << formatF(hz) << std::endl;
    }

    auto pqr = [hx, hy, hz, qj](DArray3 &p, DArray3 &q, DArray3 &r, int i, int l, int k, double x, double y, double z, double x1, double y1, double z1) {
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
        p(i,l,k) += su * (dy1 * dz1 + s1);
        p(i,l,k+1) += su * (dy1 * dz - s1);
        p(i,l+1,k) += su * (dy * dz1 - s1);
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

    Timer timer;

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

    auto copySliceX = [](DArray3 &arr, int dst, int src, int lend, int kend) {
        #pragma omp parallel for
        for (int l = 0; l < lend; l++) {
            for (int k = 0; k < kend; k++) {
                arr(dst, l, k) = arr(src, l, k);
            }
        }
    };

    auto copySliceY = [](DArray3 &arr, int dst, int src, int iend, int kend) {
        #pragma omp parallel for
        for (int i = 0; i < iend; i++) {
            for (int k = 0; k < kend; k++) {
                arr(i, dst, k) = arr(i, src, k);
            }
        }
    };

    auto copySliceZ = [](DArray3 &arr, int dst, int src, int iend, int lend) {
        #pragma omp parallel for
        for (int i = 0; i < iend; i++) {
            for (int l = 0; l < lend; l++) {
                arr(i, l, dst) = arr(i, l, src);
            }
        }
    };

    auto addSliceZ = [](DArray3 &arr, int dst, int src, int iend, int lend) {
        #pragma omp parallel for
        for (int i = 0; i < iend; i++) {
            for (int l = 0; l < lend; l++) {
                arr(i, l, dst) += arr(i, l, src);
            }
        }
    };

    auto computeBoundaryX = [im, lm, km](DArray3 &ax, const DArray3 &ay, const DArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        #pragma omp parallel for
        for (int l = 1; l < lm + 1; l++) {
            for (int k = 1; k < km + 1; k++) {
                ax(dst, l, k) = ax(src1, l, k) +
                                c1 * (ay(src2, l, k) - ay(src2, l - 1, k)) +
                                c2 * (az(src2, l, k) - az(src2, l, k - 1));
            }
        }
    };

    auto computeBoundaryY = [im, lm, km](const DArray3 &ax, DArray3 &ay, const DArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        #pragma omp parallel for
        for (int i = 1; i < im + 1; i++) {
            for (int k = 1; k < km + 1; k++) {
                ay(i, dst, k) = ay(i, src1, k) +
                                c1 * (ax(i, src2, k) - ax(i - 1, src2, k)) +
                                c2 * (az(i, src2, k) - az(i, src2, k - 1));
            }
        }
    };

    auto computeStepBlock = [rhx2, rhy2, rhz2, rc2](DArray3 &arr, const DArray3 &j, int is, int ie, int ls, int le, int ks, int ke) -> double {
        double maxdiff = 0.0;
        for (int i = is; i < ie; i++) {
            for (int l = ls; l < le; l++) {
                for (int k = ks; k < ke; k++) {
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

    auto computeStepZBlock = [rhx2, rhy2, rhz2, rc2, &computeStepBlock](DArray3 &arr, const DArray3 &j, int is, int ie, int ls, int le, int ks, int ke, int kend) -> double {
        double maxdiff = computeStepBlock(arr, j, is, ie, ls, le, ks, ke);
        if (ke == kend) {
            for (int i = is; i < ie; i++) {
                for (int l = ls; l < le; l++) {
                    const double s = ((arr(i+1,l,kend) + arr(i-1,l,kend)) * rhx2 +
                                      (arr(i,l+1,kend) + arr(i,l-1,kend)) * rhy2 +
                                      (arr(i,l,1) + arr(i,l,kend-1)) * rhz2 + j(i,l,kend)) * rc2;
                    maxdiff = std::max(std::abs(arr(i,l,kend) - s), maxdiff);
                    arr(i,l,kend) = s;
                    arr(i,l,0) = s;
                }
            }
        }
        return maxdiff;
    };

    using TaskFuncType = std::function<double(DArray3&,const DArray3&,int,int,int,int,int,int)>;

    std::vector<int> task_counter(num_of_threads, 0);
    std::vector<double> task_work_time(num_of_threads, 0);
    std::vector<double> task_crit_time(num_of_threads, 0);

    auto computeWaveTasks = [&](DArray3 &arr, const DArray3 &j, int iend, int lend, int kend, TaskFuncType &task_func) -> double {
        static const int PLACED = 1;
        static const int DONE = 2;
        const int istart = 1, lstart = 1, kstart = 1;
        const int task_size_x = std::max((iend - istart) / tasks_per_dim, min_task_size);
        const int task_size_y = std::max((lend - lstart) / tasks_per_dim, min_task_size);
        const int task_size_z = std::max((kend - kstart) / tasks_per_dim, min_task_size);
        const int num_tasks_x = std::ceil(float(iend - istart) / task_size_x);
        const int num_tasks_y = std::ceil(float(lend - lstart) / task_size_y);
        const int num_tasks_z = std::ceil(float(kend - kstart) / task_size_z);
        double maxdiff = 0.0;
        Array3D<int> state(num_tasks_x, num_tasks_y, num_tasks_z, 0);

        std::function<void(int,int,int)> doTask;
        std::function<void(int,int,int)> placeTask;

        doTask = [&](int ii, int ll, int kk) {
            Timer tm;
            const int is = istart + ii * task_size_x;
            const int ie = std::min(is + task_size_x, iend);
            const int ls = lstart + ll * task_size_y;
            const int le = std::min(ls + task_size_y, lend);
            const int ks = kstart + kk * task_size_z;
            const int ke = std::min(ks + task_size_z, kend);
            const auto mx = task_func(arr, j, is, ie, ls, le, ks, ke);
            task_work_time[omp_get_thread_num()] += tm.time();
            tm.reset();
            #pragma omp critical (doTask)
            {
                maxdiff = std::max(maxdiff, mx);
                state(ii, ll, kk) = DONE;
                // When complete, can spawn further tasks.
                placeTask(ii + 1, ll, kk);
                placeTask(ii, ll + 1, kk);
                placeTask(ii, ll, kk + 1);
            }
            task_counter[omp_get_thread_num()]++;
            task_crit_time[omp_get_thread_num()] += tm.time();
        };

        placeTask = [&](int ii, int ll, int kk) {
            if (ii >= num_tasks_x || ll >= num_tasks_y || kk >= num_tasks_z) {
                return;
            }
            // Check that task wasn't already placed.
            const auto &st = state(ii, ll, kk);
            if (st == PLACED || st == DONE) {
                return;
            }
            // Check that other required tasks completed already.
            if ((ii > 0 && state(ii - 1, ll, kk) != DONE) ||
                (ll > 0 && state(ii, ll - 1, kk) != DONE) ||
                (kk > 0 && state(ii, ll, kk - 1) != DONE)) {
                return;
            }
            state(ii, ll, kk) = PLACED;
            #pragma omp task
            doTask(ii, ll, kk);
        };

        #pragma omp parallel
        {
            #pragma omp single
            {
                // Spawn the first task.
                state(0, 0, 0) = PLACED;
                #pragma omp task
                doTask(0, 0, 0);
            }
            #pragma omp taskwait
        }
        return maxdiff;
    };

    auto computeStep = [&](DArray3 &arr, const DArray3 &j, int iend, int lend, int kend) -> double {
        TaskFuncType func = [&](DArray3 &arr, const DArray3 &j, int is, int ie, int ls, int le, int ks, int ke) {
            return computeStepBlock(arr, j, is,ie, ls, le, ks, ke);
        };
        return computeWaveTasks(arr, j, iend, lend, kend, func);
    };

    auto computeStepZ = [&](DArray3 &arr, const DArray3 &j, int iend, int lend, int kend) -> double {
        TaskFuncType func = [&](DArray3 &arr, const DArray3 &j, int is, int ie, int ls, int le, int ks, int ke) {
            return computeStepZBlock(arr, j, is,ie, ls, le, ks, ke, kend);
        };
        return computeWaveTasks(arr, j, iend, lend, kend, func);
    };

    addSliceZ(jx, 1, km + 1, im + 2, lm + 2);
    addSliceZ(jx, km, 0, im + 2, lm + 2);
    copySliceZ(jx, 0, km, im + 2, lm + 2);
    copySliceZ(jx, km + 1, 1, im + 2, lm + 2);

    addSliceZ(jy, 1, km + 1, im + 2, lm + 2);
    addSliceZ(jy, km, 0, im + 2, lm + 2);
    copySliceZ(jy, 0, km, im + 2, lm + 2);
    copySliceZ(jy, km + 1, 1, im + 2, lm + 2);

    addSliceZ(jz, km, 0, im + 2, lm + 2);
    copySliceZ(jz, 0, km, im + 2, lm + 2);

    const auto init_time = timer.time();
    timer.reset();

    int n = 0;
    double sx = 0.0, sy = 0.0, sz = 0.0;

    do {
        n++;

        sx = computeStep(ax, jx, im, lm + 1, km + 1);

        computeBoundaryX(ax, ay, az, 0, 1, 1, c12, c13);
        computeBoundaryX(ax, ay, az, im, im - 1, im, -c12, -c13);

        copySliceZ(ax, 0, km, im + 1, lm + 2);
        copySliceZ(ax, km + 1, 1, im + 1, lm + 2);

        copySliceY(ax, 0, 1, im + 1, km + 2);
        copySliceY(ax, lm + 1, lm, im + 1, km + 2);

        sy = computeStep(ay, jy, im + 1, lm, km + 1);

        computeBoundaryY(ax, ay, az, 0, 1, 1, c21, c23);
        computeBoundaryY(ax, ay, az, lm, lm - 1, lm, -c21, -c23);

        copySliceZ(ay, 0, km, im + 1, lm + 1);
        copySliceZ(ay, km + 1, 1, im + 1, lm + 1);

        copySliceX(ay, 0, 1, lm + 1, km + 2);
        copySliceX(ay, im + 1, im, lm + 1, km + 2);

        sz = computeStepZ(az, jz, im + 1, lm + 1, km);

        copySliceZ(az, 0, km, im + 2, lm + 2);
        copySliceZ(az, km+1, 1, im + 2, lm + 2);

        if (screen_output) {
            std::cout << n << " " << sx << " " << sy << " " << sz << std::endl;
        }
    } while (sx > eps || sy > eps || sz > eps);

    if (file_output) {
        out_lst << "n,sx,sy,sz=" << formatI(n) << formatS(sx) << formatS(sy) << formatS(sz) << std::endl;
    }

    const auto a_time = timer.time();
    timer.reset();

    auto computeB = [im, lm, km](int dim, const DArray3 &a1, const DArray3 &a2, DArray3 &b, double rh1, double rh2) {
        const static Index3 i_end {2, 1, 1};
        const static Index3 l_end {1, 2, 1};
        const static Index3 k_end {1, 1, 2};
        const static std::array<Index3, 3> shift1 {Index3 {0, 1, 0}, Index3 {0, 0, 1}, Index3 {1, 0, 0}};
        const static std::array<Index3, 3> shift2 {Index3 {0, 0, 1}, Index3 {1, 0, 0}, Index3 {0, 1, 0}};

        const auto &s1 = shift1[dim];
        const auto &s2 = shift2[dim];
        #pragma omp parallel for
        for (int k = 0; k < km + k_end[dim]; k++) {
            for (int l = 0; l < lm + l_end[dim]; l++) {
                for (int i = 0; i < im + i_end[dim]; i++) {
                    b(i,l,k) = (a1(i+s1[0],l+s1[1],k+s1[2]) - a1(i,l,k)) * rh1 -
                               (a2(i+s2[0],l+s2[1],k+s2[2]) - a2(i,l,k)) * rh2;
                }
            }
        }
    };

    computeB(0, az, ay, bx, rhy, rhz);
    computeB(1, ax, az, by, rhz, rhx);
    computeB(2, ay, ax, bz, rhx, rhy);

    const auto b_time = timer.time();
    timer.reset();

    auto outputMax = [&file_output, &out_lst](const std::string &str, double m) {
        if (file_output) {
            out_lst << str << formatS(m) << std::endl;
        }
    };

    auto outputMaxI = [&file_output, &out_lst](const std::string &str, double m, const std::array<int, 3> &ind) {
        if (file_output) {
            out_lst << str << formatS(m) << formatI(ind[0] + 1) << formatI(ind[1] + 1) << formatI(ind[2] + 1) << std::endl;
        }
    };

    auto computeDiv = [im, lm, km, rhx, rhy, rhz](const DArray3 &x, const DArray3 &y, const DArray3 &z, int start, int shift1, int shift2)
        -> std::pair<double, Index3> {
        double maxval = 0.0;
        Index3 maxind {0, 0, 0};
        #pragma omp parallel
        {
            double maxval_l = 0.0;
            Index3 maxind_l {0, 0, 0};
            #pragma omp for
            for (int k = start; k < km + 1; k++) {
                for (int l = start; l < lm + 1; l++) {
                    for (int i = start; i < im + 1; i++) {
                        const double s = (x(i + shift1, l, k) - x(i + shift2, l, k)) * rhx +
                                         (y(i, l + shift1, k) - y(i, l + shift2, k)) * rhy +
                                         (z(i, l, k + shift1) - z(i, l, k + shift2)) * rhz;
                        if (std::abs(s) > maxval_l) {
                            maxval_l = std::abs(s);
                            maxind_l = {i, l, k};
                        }
                    }
                }
            }
            #pragma omp critical
            {
                if (maxval_l > maxval) {
                    maxval = maxval_l;
                    maxind = maxind_l;
                }
            }
        }
        return std::make_pair(maxval, maxind);
    };

    const auto mdj = computeDiv(jx, jy, jz, 1, 0, -1);
    outputMax("max(divj)=", mdj.first);

    const auto mdb = computeDiv(bx, by, bz, 0, 1, 0);
    outputMax("max(divB)=", mdb.first);

    const auto mda = computeDiv(ax, ay, az, 1, 0, -1);
    outputMaxI("max(divA)=", mda.first, mda.second);

    const auto div_time = timer.time();
    timer.reset();

    auto computeRotDiff = [im, lm, km](int dim, const DArray3 &b1, const DArray3 &b2, const DArray3 &j, double rh1, double rh2) ->
    std::pair<double, Index3> {
        const static Index3 i_start {0, 1, 1};
        const static Index3 l_start {1, 0, 1};
        const static Index3 k_start {1, 1, 0};
        const static std::array<Index3, 3> shift1 {Index3 {0, -1, 0}, Index3 {0, 0, -1}, Index3 {-1, 0, 0}};
        const static std::array<Index3, 3> shift2 {Index3 {0, 0, -1}, Index3 {-1, 0, 0}, Index3 {0, -1, 0}};

        double maxval = 0.0;
        Index3 maxind {0, 0, 0};
        const auto &s1 = shift1[dim];
        const auto &s2 = shift2[dim];
        #pragma omp parallel
        {
            double maxval_l = 0.0;
            Index3 maxind_l {0, 0, 0};
            #pragma omp for
            for (int k = k_start[dim]; k < km + 1; k++) {
                for (int l = l_start[dim]; l < lm + 1; l++) {
                    for (int i = i_start[dim]; i < im + 1; i++) {
                        const double s = (b1(i,l,k) - b1(i+s1[0],l+s1[1],k+s1[2])) * rh1 -
                                         (b2(i,l,k) - b2(i+s2[0],l+s2[1],k+s2[2])) * rh2 - j(i,l,k);                        
                        if (std::abs(s) > maxval_l) {
                            maxval_l = std::abs(s);
                            maxind_l = Index3 {i, l, k};
                        }
                    }
                }
            }
            #pragma omp critical (maxRotDiff)
            {                
                if (maxval_l > maxval) {
                    maxval = maxval_l;
                    maxind = maxind_l;
                }
            }
        }
        return std::make_pair(maxval, maxind);
    };

    const auto mrx = computeRotDiff(0, bz, by, jx, rhy, rhz);
    outputMaxI("max(rotB_x-jx)=", mrx.first, mrx.second);

    const auto mry = computeRotDiff(1, bx, bz, jy, rhz, rhx);
    outputMaxI("max(rotB_y-jy)=", mry.first, mry.second);

    const auto mrz = computeRotDiff(2, by, bx, jz, rhx, rhy);
    outputMaxI("max(rotB_z-jz)=", mrz.first, mrz.second);

    const auto rot_time = timer.time();

    const auto work_time = init_time + a_time + b_time + div_time + rot_time;
    std::cout << "THREADS: " << num_of_threads << ", TIME: " << work_time << endl;
    std::cout << "Init: " << init_time <<
        ", A: " << a_time <<
        ", B: " << b_time <<
        ", Div: " << div_time <<
        ", Rot: " << rot_time << std::endl;
    for (int i = 0; i < num_of_threads; i++) {
        std::cout << "Thread " << i << ": tasks: " << task_counter[i] <<
            ", work time: " << task_work_time[i] <<
            ", crit time: " << task_crit_time[i] <<
            std::endl;
    }

    if (file_output) {
        out_lst.close();
    }

    if (file_output) {
        outputDat("jx.dat", jx, im + 1, lm + 2, km + 2);
        outputDat("jy.dat", jy, im + 2, lm + 1, km + 2);
        outputDat("jz.dat", jz, im + 2, lm + 2, km + 1);
        outputDat("ax.dat", ax, im + 1, lm + 2, km + 2);
        outputDat("ay.dat", ay, im + 2, lm + 1, km + 2);
        outputDat("az.dat", az, im + 2, lm + 2, km + 1);
        outputDat("bx.dat", bx, im + 2, lm + 1, km + 1);
        outputDat("by.dat", by, im + 1, lm + 2, km + 1);
        outputDat("bz.dat", bz, im + 1, lm + 1, km + 2);
    }

    return 0;
}
