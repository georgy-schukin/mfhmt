#include "defs.h"
#include "../common/output.h"
#include "../common/timer.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>

using namespace std;

int main(int argc, char **argv) {

    const int IM_DEF = 40;
    const int LM_DEF = 40;
    const int KM_DEF = 20;
    const int NM_DEF = 1000;
    const int FOUT_DEF = 1;
    const int SOUT_DEF = 0;

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
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int lm = (argc > 2) ? stoi(argv[2]) : LM_DEF;
    const int km = (argc > 3) ? stoi(argv[3]) : KM_DEF;
    const int nm = (argc > 4) ? stoi(argv[4]) : NM_DEF;
    const bool file_output = (argc > 5) ? stoi(argv[5]) : FOUT_DEF;
    const bool screen_output = (argc > 6) ? stoi(argv[6]) : SOUT_DEF;

    const size_t ims = im + 2;
    const size_t lms = lm + 2;
    const size_t kms = km + 2;

/*
    real*8 jx(im+2,lm+2,km+2),jy(im+2,lm+2,km+2),jz(im+2,lm+2,km+2),
    ax(im+2,lm+2,km+2),ay(im+2,lm+2,km+2),az(im+2,lm+2,km+2),
    bx(im+2,lm+2,km+2),by(im+2,lm+2,km+2),bz(im+2,lm+2,km+2)
*/

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
        out_lst.open("output.lst");
        out_lst << "qj=" << formatF(qj) << std::endl;
        out_lst << "1. im,lm,km,nm=" << formatI(im) << formatI(lm) << formatI(km) << formatI(nm) << std::endl;
        out_lst << "1. hx,hy,hz=" << formatF(hx) << formatF(hy) << formatF(hz) << std::endl;
    }

    // Init functions

/*
    subroutine pqr(i,l,k,h1,h2,h3,x,y,z,x1,y1,z1,a)
      parameter(im=40,lm=40,km=20)
      real*8 p(im+2,lm+2,km+2),q(im+2,lm+2,km+2),r(im+2,lm+2,km+2)
      integer i,l,k
      real*8 h1,h2,h3,x,y,z,x1,y1,z1,a,
     *dx,dy,dz,dx1,dy1,dz1,su,sv,sw,s1,s2,s3
      common/j/p,q,r
      dx=0.5*(x+x1)-h1*(i-1.5)
      dy=0.5*(y+y1)-h2*(l-1.5)
      dz=0.5*(z+z1)-h3*(k-1.5)
      dx1=h1-dx
      dy1=h2-dy
      dz1=h3-dz
      su=x1-x
      sv=y1-y
      sw=z1-z
      s1=sv*sw/12.
      s2=su*sw/12.
      s3=su*sv/12.
      su=su*a/(h2*h3)   !
      sv=sv*a/(h1*h3)   !
      sw=sw*a/(h1*h2)   !
      p(i,l,k)=p(i,l,k)+su*(dy1*dz1+s1)
      p(i,l,k+1)=p(i,l,k+1)+su*(dy1*dz-s1)
      p(i,l+1,k)=p(i,l+1,k)+su*(dy*dz1-s1)
      p(i,l+1,k+1)=p(i,l+1,k+1)+su*(dy*dz+s1)
      q(i,l,k)=q(i,l,k)+sv*(dx1*dz1+s2)
      q(i,l,k+1)=q(i,l,k+1)+sv*(dx1*dz-s2)
      q(i+1,l,k)=q(i+1,l,k)+sv*(dx*dz1-s2)
      q(i+1,l,k+1)=q(i+1,l,k+1)+sv*(dx*dz+s2)
      r(i,l,k)=r(i,l,k)+sw*(dx1*dy1+s3)
      r(i,l+1,k)=r(i,l+1,k)+sw*(dx1*dy-s3)
      r(i+1,l,k)=r(i+1,l,k)+sw*(dx*dy1-s3)
      r(i+1,l+1,k)=r(i+1,l+1,k)+sw*(dx*dy+s3)
      return
      end
*/
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

/*
    zadanie vintovogo toka

      do k=1,km+2
         do l=1,lm+2
            do i=1,im+2
               jx(i,l,k)=0.d0
               jy(i,l,k)=0.d0
               jz(i,l,k)=0.d0
            enddo
         enddo
      enddo

      x=x0+r0
      y=y0
      z=0.d0
      do 111 n=1,nm
      x1=x0+r0*dcos(c1*n)
      y1=y0+r0*dsin(c1*n)
      z1=h0*c1*n
      s2=x/hx
      i1=idint(s2+1.5)
      s4=y/hy
      l1=idint(s4+1.5)
      s6=z/hz
      k1=idint(s6+1.5)
      s2=x1/hx
      i2=idint(s2+1.5)
      s4=y1/hy
      l2=idint(s4+1.5)
      s6=z1/hz
      k2=idint(s6+1.5)
      i=abs(i2-i1)
      l=abs(l2-l1)
      k=abs(k2-k1)
      m=4*i+2*l+k
      goto(1,2,3,4,5,6,7),m
      call pqr(i1,l1,k1,hx,hy,hz,x,y,z,x1,y1,z1,qj)
      goto 18
    1 z2=hz*(0.5*(k1+k2)-1.)
      s=(z2-z)/(z1-z)
      x2=x+(x1-x)*s
      y2=y+(y1-y)*s
      goto 11
    2 y2=hy*(0.5*(l1+l2)-1.)
      s=(y2-y)/(y1-y)
      x2=x+(x1-x)*s
      z2=z+(z1-z)*s
      goto 11
    3 y2=hy*(0.5*(l1+l2)-1.)
      z2=hz*(0.5*(k1+k2)-1.)
      s=((z1-z)*(z2-z)+(y1-y)*(y2-y))/((z1-z)**2+(y1-y)**2)
      x2=x+(x1-x)*s
      goto 11
    4 x2=hx*(0.5*(i1+i2)-1.)
      s=(x2-x)/(x1-x)
      y2=y+(y1-y)*s
      z2=z+(z1-z)*s
      goto 11
    5 x2=hx*(0.5*(i1+i2)-1.)
      z2=hz*(0.5*(k1+k2)-1.)
      s=((z1-z)*(z2-z)+(x1-x)*(x2-x))/((z1-z)**2+(x1-x)**2)
      y2=y+(y1-y)*s
      goto 11
    6 x2=hx*(0.5*(i1+i2)-1.)
      y2=hy*(0.5*(l1+l2)-1.)
      s=((y1-y)*(y2-y)+(x1-x)*(x2-x))/((y1-y)**2+(x1-x)**2)
      z2=z+(z1-z)*s
      goto 11
    7 x2=hx*(0.5*(i1+i2)-1.)
      y2=hy*(0.5*(l1+l2)-1.)
      z2=hz*(0.5*(k1+k2)-1.)
   11 call pqr(i1,l1,k1,hx,hy,hz,x,y,z,x2,y2,z2,qj)
      call pqr(i2,l2,k2,hx,hy,hz,x2,y2,z2,x1,y1,z1,qj)
   18 continue

      x=x1
      y=y1
      z=z1

  111 continue
*/
    auto initHelicalCurrent = [x0, y0, r0, h0, nm, hx, hy, hz, rhx, rhy, rhz, c1, &pqr](DArray3 &jx, DArray3 &jy, DArray3 &jz) {
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

/*
    vychislenie vektornogo potentsiala

    eps=1.d-10
    c2=2.d0/hx**2+2.d0/hy**2+2.d0/hz**2
    c12=hx/hy
    c13=hx/hz
    c21=hy/hx
    c23=hy/hz
*/

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
        for (int l = 0; l < lend; l++) {
            for (int k = 0; k < kend; k++) {
                arr(dst, l, k) = arr(src, l, k);
            }
        }
    };

    auto copySliceY = [](DArray3 &arr, int dst, int src, int iend, int kend) {
        for (int i = 0; i < iend; i++) {
            for (int k = 0; k < kend; k++) {
                arr(i, dst, k) = arr(i, src, k);
            }
        }
    };

    auto copySliceZ = [](DArray3 &arr, int dst, int src, int iend, int lend) {
        for (int i = 0; i < iend; i++) {
            for (int l = 0; l < lend; l++) {
                arr(i, l, dst) = arr(i, l, src);
            }
        }
    };

    auto addSliceZ = [](DArray3 &arr, int dst, int src, int iend, int lend) {
        for (int i = 0; i < iend; i++) {
            for (int l = 0; l < lend; l++) {
                arr(i, l, dst) += arr(i, l, src);
            }
        }
    };

/*
    do l=1,lm+2
        do i=1,im+2
            jx(i,l,2)=jx(i,l,2)+jx(i,l,km+2)
            jx(i,l,km+1)=jx(i,l,km+1)+jx(i,l,1)
            jx(i,l,1)=jx(i,l,km+1)
            jx(i,l,km+2)=jx(i,l,2)

            jy(i,l,2)=jy(i,l,2)+jy(i,l,km+2)
            jy(i,l,km+1)=jy(i,l,km+1)+jy(i,l,1)
            jy(i,l,1)=jy(i,l,km+1)
            jy(i,l,km+2)=jy(i,l,2)

            jz(i,l,km+1)=jz(i,l,km+1)+jz(i,l,1)
            jz(i,l,1)=jz(i,l,km+1)
        enddo
    enddo
*/

    auto updateCurrent = [im, lm, km, &addSliceZ, &copySliceZ](DArray3 &jx, DArray3 &jy, DArray3 &jz) {
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
    };

/*
      do k=2,km+1
         do l=2,lm+1
            ax(1,l,k)=ax(2,l,k)+c12*(ay(2,l,k)-ay(2,l-1,k))+
     =               c13*(az(2,l,k)-az(2,l,k-1))
            ax(im+1,l,k)=ax(im,l,k)-c12*(ay(im+1,l,k)-ay(im+1,l-1,k))-
     =               c13*(az(im+1,l,k)-az(im+1,l,k-1))
         enddo
      enddo
*/

    auto computeBoundaryX = [im, lm, km](DArray3 &ax, DArray3 &ay, DArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        for (int l = 1; l < lm + 1; l++) {
            for (int k = 1; k < km + 1; k++) {
                ax(dst, l, k) = ax(src1, l, k) +
                                c1 * (ay(src2, l, k) - ay(src2, l - 1, k)) +
                                c2 * (az(src2, l, k) - az(src2, l, k - 1));
            }
        }
    };

/*
      do k=2,km+1
         do i=2,im+1
            ay(i,1,k)=ay(i,2,k)+c21*(ax(i,2,k)-ax(i-1,2,k))+
     =          c23*(az(i,2,k)-az(i,2,k-1))
            ay(i,lm+1,k)=ay(i,lm,k)-c21*(ax(i,lm+1,k)-ax(i-1,lm+1,k))-
     =          c23*(az(i,lm+1,k)-az(i,lm+1,k-1))
         enddo
      enddo
*/

    auto computeBoundaryY = [im, lm, km](DArray3 &ax, DArray3 &ay, DArray3 &az, int dst, int src1, int src2, double c1, double c2) {
        for (int i = 1; i < im + 1; i++) {
            for (int k = 1; k < km + 1; k++) {
                ay(i, dst, k) = ay(i, src1, k) +
                                c1 * (ax(i, src2, k) - ax(i - 1, src2, k)) +
                                c2 * (az(i, src2, k) - az(i, src2, k - 1));
            }
        }
    };

/*
      sx=0.d0
      do k=2,km+1
         do l=2,lm+1
            do i=2,im
               s=((ax(i+1,l,k)+ax(i-1,l,k))/hx**2+
     =            (ax(i,l+1,k)+ax(i,l-1,k))/hy**2+
     =            (ax(i,l,k+1)+ax(i,l,k-1))/hz**2+jx(i,l,k))/c2
               s2=dabs(ax(i,l,k)-s)
               if(s2.gt.sx) sx=s2
               ax(i,l,k)=s
            enddo
         enddo
      enddo
*/

    auto computeStep = [rhx2, rhy2, rhz2, rc2](DArray3 &arr, const DArray3 &j, int iend, int lend, int kend) -> double {
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

/*
      sz=0.d0
      do l=2,lm+1
         do i=2,im+1
            do k=2,km
               s=((az(i+1,l,k)+az(i-1,l,k))/hx**2+
     =            (az(i,l+1,k)+az(i,l-1,k))/hy**2+
     =            (az(i,l,k+1)+az(i,l,k-1))/hz**2+jz(i,l,k))/c2
               s2=dabs(az(i,l,k)-s)
               if(s2.gt.sz) sz=s2
               az(i,l,k)=s
            enddo
            s=((az(i+1,l,km+1)+az(i-1,l,km+1))/hx**2+
     =         (az(i,l+1,km+1)+az(i,l-1,km+1))/hy**2+
     =         (az(i,l,2)+az(i,l,km))/hz**2+jz(i,l,km+1))/c2
            s2=dabs(az(i,l,km+1)-s)
            if(s2.gt.sz) sz=s2
            az(i,l,km+1)=s
            az(i,l,1)=s
         enddo
      enddo
*/

    auto computeStepZ = [rhx2, rhy2, rhz2, rc2](DArray3 &arr, const DArray3 &j, int iend, int lend, int kend) -> double {
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

/*
      do k=1,km+1
         do l=1,lm+1
            do i=1,im+2
               bx(i,l,k)=(az(i,l+1,k)-az(i,l,k))/hy-
     =                   (ay(i,l,k+1)-ay(i,l,k))/hz
            enddo
         enddo
      enddo
*/

    auto computeB = [im, lm, km](int dim, const DArray3 &a1, const DArray3 &a2, DArray3 &b, double rh1, double rh2) {
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

/*
      s=0.d0
      do k=2,km+1
         do l=2,lm+1
            do i=2,im+1
               s1=(jx(i,l,k)-jx(i-1,l,k))/hx+
     =            (jy(i,l,k)-jy(i,l-1,k))/hy+
     =            (jz(i,l,k)-jz(i,l,k-1))/hz
               if(dabs(s1).gt.s) s=s1
            enddo
         enddo
      enddo
*/

    auto computeDiv = [im, lm, km, rhx, rhy, rhz](const DArray3 &x, const DArray3 &y, const DArray3 &z, int start, int shift1, int shift2)
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

/*
    s1=0.d0
    s2=0.d0
    s3=0.d0

    do k=2,km+1
     do l=2,lm+1
        do i=1,im+1
           s4=(bz(i,l,k)-bz(i,l-1,k))/hy-
    =            (by(i,l,k)-by(i,l,k-1))/hz-jx(i,l,k)
           s=dabs(s4)
           if(s.gt.s1) then
              s1=s4
              i1=i
              l1=l
              k1=k
           endif
        enddo
     enddo
    enddo
*/

    auto computeRotDiff = [im, lm, km](int dim, const DArray3 &b1, const DArray3 &b2, const DArray3 &j, double rh1, double rh2)
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
            std::cout << iter << " " << sx << " " << sy << " " << sz << std::endl;
        }
    } while (sx > eps || sy > eps || sz > eps);

    if (file_output) {
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

    std::cout << "Im: " << im << ", Lm: " << lm << ", Km: " << km << ", Nm: " << nm << std::endl;
    std::cout << "TIME: " << work_time << std::endl;
    std::cout << "Iters: " << iter << std::endl;

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
