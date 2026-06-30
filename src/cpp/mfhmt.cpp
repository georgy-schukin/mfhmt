#include "common.h"
#include "output.h"

#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <algorithm>

using namespace std;

int main(int argc, char **argv) {
/*
    program simit44
*/

    const int IM_DEF = 40;
    const int LM_DEF = 40;
    const int KM_DEF = 20;
    const int FOUT_DEF = 1;
    const int FULL_OUTPUT_DEF = 0;

    if (argc > 1) {
        const auto s = string(argv[1]);
        if (s == "-h" || s == "--help") {
            cout << argv[0] <<
                " [im=" << IM_DEF << "]" <<
                " [lm=" << LM_DEF << "]" <<
                " [km=" << KM_DEF << "]" <<
                " [file_output=" << FOUT_DEF << "]" <<
                " [full_output=" << FULL_OUTPUT_DEF << "]" <<
                endl;
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int lm = (argc > 2) ? stoi(argv[2]) : LM_DEF;
    const int km = (argc > 3) ? stoi(argv[3]) : KM_DEF;
    const bool file_output = (argc > 4) ? stoi(argv[4]) : FOUT_DEF;
    const bool full_output = (argc > 5) ? stoi(argv[5]) : FULL_OUTPUT_DEF;

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
    const int nm = 1000;
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

/*
      open(25,file='smt45.lst',form='formatted')

      write(25,*)'qj=',qj

      print*,'1. im,lm,km,nm=',im,lm,km,nm
      write(25,100) im,lm,km,nm
  100 format('1. im,lm,km,nm=',4i6)
      write(25,102) hx,hy,hz
  102 format('1. hx,hy,hz=',3f10.3)
*/

/*
c==============================================================
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
        p(i,l,k) = p(i,l,k) + su * (dy1 * dz1 + s1);
        p(i,l,k+1) = p(i,l,k+1) + su * (dy1 * dz - s1);
        p(i,l+1,k) = p(i,l+1,k) + su * (dy * dz1 - s1);
        p(i,l+1,k+1) = p(i,l+1,k+1) + su * (dy * dz + s1);
        q(i,l,k) = q(i,l,k) + sv * (dx1 * dz1 + s2);
        q(i,l,k+1) = q(i,l,k+1) + sv * (dx1 * dz - s2);
        q(i+1,l,k) = q(i+1,l,k) + sv * (dx * dz1 - s2);
        q(i+1,l,k+1) = q(i+1,l,k+1) + sv * (dx * dz + s2);
        r(i,l,k) = r(i,l,k) + sw * (dx1 * dy1 + s3);
        r(i,l+1,k) = r(i,l+1,k) + sw * (dx1 * dy - s3);
        r(i+1,l,k) = r(i+1,l,k) + sw * (dx * dy1 - s3);
        r(i+1,l+1,k) = r(i+1,l+1,k) + sw * (dx * dy + s3);
    };


/*
c zadanie vintovogo toka

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

    for (int l = 0; l < lm + 2; l++) {
        for (int i = 0; i < im + 2; i++) {
            jx(i,l,1) = jx(i,l,1) + jx(i,l,km+1);
            jx(i,l,km) = jx(i,l,km) + jx(i,l,0);
            jx(i,l,0) = jx(i,l,km);
            jx(i,l,km + 1) = jx(i,l,1);

            jy(i,l,1) = jy(i,l,1) + jy(i,l,km+1);
            jy(i,l,km) = jy(i,l,km) + jy(i,l,0);
            jy(i,l,0) = jy(i,l,km);
            jy(i,l,km+1) = jy(i,l,1);

            jz(i,l,km) = jz(i,l,km) + jz(i,l,0);
            jz(i,l,0) = jz(i,l,km);
        }
    }

/*
c   vychislenie vektornogo potentsiala

      eps=1.d-10
      c2=2.d0/hx**2+2.d0/hy**2+2.d0/hz**2
      c12=hx/hy
      c13=hx/hz
      c21=hy/hx
      c23=hy/hz

      do k=1,km+2
         do l=1,lm+2
            do i=1,im+2
               ax(i,l,k)=0.d0
               ay(i,l,k)=0.d0
               az(i,l,k)=0.d0
            enddo
         enddo
      enddo

      n=0

    8 n=n+1
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

/*
c  ax
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
    double sx = 0.0;
    for (int k = 1; k < km + 1; k++) {
        for (int l = 1; l < lm + 1; l++) {
            for (int i = 1; i < im; i++) {
                const double s = ((ax(i+1,l,k) + ax(i-1,l,k)) * rhx2 +
                            (ax(i,l+1,k) + ax(i,l-1,k)) * rhy2 +
                            (ax(i,l,k+1) + ax(i,l,k-1)) * rhz2 + jx(i,l,k)) * rc2;
                sx = std::max(std::abs(ax(i,l,k) - s), sx);
                ax(i,l,k) = s;
            }
        }
    }


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
    for (int k = 1; k < km + 1; k++) {
        for (int l = 1; l < lm + 1; l++) {
            ax(0,l,k) = ax(1,l,k) +
                        c12 * (ay(1,l,k) - ay(1,l-1,k)) +
                        c13 * (az(1,l,k) - az(1,l,k-1));
            ax(im,l,k) = ax(im,l,k) -
                        c12 * (ay(im,l,k) - ay(im,l-1,k)) -
                        c13 * (az(im,l,k) - az(im,l,k-1));
        }
    }

/*
      do l=1,lm+2
         do i=1,im+1
            ax(i,l,1)=ax(i,l,km+1)
            ax(i,l,km+2)=ax(i,l,2)
         enddo
      enddo
*/
    for (int l = 0; l < lm + 2; l++) {
        for (int i = 0; i < im + 1; i++) {
            ax(i,l,0) = ax(i,l,km);
            ax(i,l,km+1) = ax(i,l,1);
        }
    }

/*
      do k=1,km+2
         do i=1,im+1
            ax(i,1,k)=ax(i,2,k)             ! ?
            ax(i,lm+2,k)=ax(i,lm+1,k)       ! ?
         enddo
      enddo
*/
    for (int k = 0; k < km + 2; k++) {
        for (int i = 0; i < im + 1; i++) {
            ax(i,0,k) = ax(i,1,k);
            ax(i,lm+1,k) = ax(i,lm,k);
        }
    }

/*
c  ay
      sy=0.d0
      do k=2,km+1
         do l=2,lm
            do i=2,im+1
               s=((ay(i+1,l,k)+ay(i-1,l,k))/hx**2+
     =            (ay(i,l+1,k)+ay(i,l-1,k))/hy**2+
     =            (ay(i,l,k+1)+ay(i,l,k-1))/hz**2+jy(i,l,k))/c2
               s2=dabs(ay(i,l,k)-s)
               if(s2.gt.sy) sy=s2
               ay(i,l,k)=s
            enddo
         enddo
      enddo
*/
    double sy = 0.0;
    for (int k = 1; k < km + 1; k++) {
        for (int l = 1; l < lm; l++) {
            for (int i = 1; i < im + 1; i++) {
                const double s = ((ay(i+1,l,k) + ay(i-1,l,k)) * rhx2 +
                                  (ay(i,l+1,k) + ay(i,l-1,k)) * rhy2 +
                                  (ay(i,l,k+1) + ay(i,l,k-1)) * rhz2 + jy(i,l,k)) * rc2;
                sy = std::max(std::abs(ay(i,l,k) - s), sy);
                ay(i,l,k) = s;
            }
        }
    }

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
    for (int k = 1; k < km + 1; k++) {
        for (int i = 1; i < im + 1; i++) {
            ay(i,0,k) = ay(i,1,k) +
                        c21 * (ax(i,1,k) - ax(i-1,1,k)) +
                        c23 * (az(i,1,k) - az(i,1,k-1));
            ay(i,lm,k) = ay(i,lm,k) -
                        c21 * (ax(i,lm,k) - ax(i-1,lm,k)) -
                        c23 * (az(i,lm,k) - az(i,lm,k-1));
        }
    }

/*
      do l=1,lm+1
         do i=1,im+2
            ay(i,l,1)=ay(i,l,km+1)
            ay(i,l,km+2)=ay(i,l,2)
         enddo
      enddo
*/
    for (int l = 0; l < lm + 1; l++) {
        for (int i = 0; i < im + 1; i++) {
            ay(i,l,0) = ay(i,l,km);
            ay(i,l,km+1) = ay(i,l,1);
        }
    }

/*
      do k=1,km+2
         do l=1,lm+1
            ay(1,l,k)=ay(2,l,k)        ! ?
            ay(im+2,l,k)=ay(im+1,l,k)  ! ?
         enddo
      enddo
*/
    for (int k = 0; k < km + 2; k++) {
        for (int l = 0; l < lm + 1; l++) {
            ay(0,l,k) = ay(1,l,k);
            ay(im+1,l,k) = ay(im,l,k);
        }
    }

/*
c  az
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
    double sz = 0.0;
    for (int l = 1; l < lm + 1; l++) {
        for (int i = 1; i < im + 1; i++) {
            for (int k = 1; k < km; k++) {
                const double s = ((az(i+1,l,k) + az(i-1,l,k)) * rhx2 +
                                  (az(i,l+1,k) + az(i,l-1,k)) * rhy2 +
                                  (az(i,l,k+1) + az(i,l,k-1)) * rhz2 + jz(i,l,k)) * rc2;
                sz = std::max(std::abs(az(i,l,k) - s), sz);
                az(i,l,k) = s;
            }
            const double s = ((az(i+1,l,km) + az(i-1,l,km)) * rhx2 +
                              (az(i,l+1,km) + az(i,l-1,km)) * rhy2 +
                              (az(i,l,1) + az(i,l,km)) * rhz2 + jz(i,l,km)) * rc2;
            sz = std::max(std::abs(az(i,l,km) - s), sz);
            az(i,l,km) = s;
            az(i,l,0) = s;
        }
    }

/*
      do l=1,lm+2
         do i=1,im+2
            az(i,l,1)=az(i,l,km+1)
            az(i,l,km+2)=az(i,l,2)
         enddo
      enddo
*/
    for (int l = 0; l < lm + 2; l++) {
        for (int i = 0; i < im + 2; i++) {
            az(i,l,0) = az(i,l,km);
            az(i,l,km+1) = az(i,l,1);
        }
    }

/*
c      write(25,*) 'ax. n,sx=',n,sx
      print 104,n,sx,sy,sz
  104 format('n,sx,sy,sz=',i6,3e12.4)

      if((sx.gt.eps).or.(sy.gt.eps).or.(sz.gt.eps)) goto 8
c      if(n.lt.200) goto 8

      write(25,104) n,sx,sy,sz
*/

/*
c========================= bx, by, bz ============

      do k=1,km+1
         do l=1,lm+1
            do i=1,im+2
               bx(i,l,k)=(az(i,l+1,k)-az(i,l,k))/hy-
     =                   (ay(i,l,k+1)-ay(i,l,k))/hz
            enddo
         enddo
      enddo
*/
    for (int k = 0; k < km + 1; k++) {
        for (int l = 0; l < lm + 1; l++) {
            for (int i = 0; i < im + 2; i++) {
                bx(i,l,k) = (az(i,l+1,k) - az(i,l,k)) * rhy -
                            (ay(i,l,k+1) - ay(i,l,k)) * rhz;
            }
        }
    }

/*
      do k=1,km+1
         do l=1,lm+2
            do i=1,im+1
               by(i,l,k)=(ax(i,l,k+1)-ax(i,l,k))/hz-
     =                   (az(i+1,l,k)-az(i,l,k))/hx
            enddo
         enddo
      enddo
*/
    for (int k = 0; k < km + 1; k++) {
        for (int l = 0; l < lm + 2; l++) {
            for (int i = 0; i < im + 1; i++) {
                by(i,l,k) = (ax(i,l,k+1) - ax(i,l,k)) * rhz -
                            (az(i+1,l,k) - az(i,l,k)) * rhx;
            }
        }
    }

/*
      do k=1,km+2
         do l=1,lm+1
            do i=1,im+1
               bz(i,l,k)=(ay(i+1,l,k)-ay(i,l,k))/hx-
     =                   (ax(i,l+1,k)-ax(i,l,k))/hy
            enddo
         enddo
      enddo
*/
    for (int k = 0; k < km + 2; k++) {
        for (int l = 0; l < lm + 1; l++) {
            for (int i = 0; i < im + 1; i++) {
                bz(i,l,k) = (ay(i+1,l,k) - ay(i,l,k)) * rhx -
                            (ax(i,l+1,k) - ax(i,l,k)) * rhy;
            }
        }
    }

/*
c===================================== divj
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
      write(25,*) 'max(divj)=',s
*/
    double s = 0.0;
    for (int k = 1; k < km + 1; k++) {
        for (int l = 1; l < lm + 1; l++) {
            for (int i = 1; i < im + 1; i++) {
                const double s1 = (jx(i,l,k) - jx(i-1,l,k)) * rhx +
                                (jy(i,l,k) - jy(i,l-1,k)) * rhy +
                                (jz(i,l,k) - jz(i,l,k-1)) * rhz;
                if (std::abs(s1) > s) {
                    s = s1;
                }
            }
        }
    }
/*
c===================================== divB
      s=0.d0
      do k=1,km+1
         do l=1,lm+1
            do i=1,im+1
               s1=(bx(i+1,l,k)-bx(i,l,k))/hx+
     =            (by(i,l+1,k)-by(i,l,k))/hy+
     =            (bz(i,l,k+1)-bz(i,l,k))/hz
               if(dabs(s1).gt.s) s=s1
            enddo
         enddo
      enddo
      write(25,*) 'max(divB)=',s
*/
    s = 0.0;
    for (int k = 0; k < km + 1; k++) {
        for (int l = 0; l < lm + 1; l++) {
            for (int i = 0; i < im + 1; i++) {
                const double s1 = (bx(i+1,l,k) - bx(i,l,k)) * rhx +
                                  (by(i,l+1,k) - by(i,l,k)) * rhy +
                                  (bz(i,l,k+1) - bz(i,l,k)) * rhz;
                if (std::abs(s1) > s) {
                    s = s1;
                }
            }
        }
    }
/*
c===================================== divA
      s1=0.d0
      do k=2,km+1
         do l=2,lm+1
            do i=2,im+1
               s2=(ax(i,l,k)-ax(i-1,l,k))/hx+
     =            (ay(i,l,k)-ay(i,l-1,k))/hy+
     =            (az(i,l,k)-az(i,l,k-1))/hz
               s=dabs(s2)
               if(s.gt.s1) then
                  s1=s
                  i1=i
                  l1=l
                  k1=k
               endif
            enddo
         enddo
      enddo
      write(25,*) 'max(divA)=',s1,i1,l1,k1
*/
    s = 0.0;
    std::array<int, 3> maxind {0, 0, 0};
    for (int k = 1; k < km + 1; k++) {
        for (int l = 1; l < lm + 1; l++) {
            for (int i = 1; i < im + 1; i++) {
                const double s1 = (ax(i,l,k) - ax(i-1,l,k)) * rhx +
                                  (ay(i,l,k) - ay(i,l-1,k)) * rhy +
                                  (az(i,l,k) - az(i,l,k-1)) * rhz;
                if (std::abs(s1) > s) {
                    s = s1;
                    maxind = {i, l, k};
                }
            }
        }
    }

/*
c===================================== rotB-j
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
      write(25,105) s1,i1,l1,k1
  105 format(' max(rotB_x-jx)=',e12.4,3i4)
*/
    s = 0.0;
    maxind = {0, 0, 0};
    for (int k = 1; k < km + 1; k++) {
        for (int l = 1; l < lm + 1; l++) {
            for (int i = 0; i < im + 1; i++) {
                const double s1 = (bz(i,l,k) - bz(i,l-1,k)) * rhy +
                                  (by(i,l,k) - by(i,l,k-1)) * rhz - jx(i,l,k);
                if (std::abs(s1) > s) {
                    s = s1;
                    maxind = {i, l, k};
                }
            }
        }
    }

/*
      do k=2,km+1
         do l=1,lm+1
            do i=2,im+1
               s4=(bx(i,l,k)-bx(i,l,k-1))/hz-
     =            (bz(i,l,k)-bz(i-1,l,k))/hx-jy(i,l,k)
               if(dabs(s4).gt.s2) then
                  s2=s4
                  i1=i
                  l1=l
                  k1=k
               endif
            enddo
         enddo
      enddo
      write(25,106) s2,i1,l1,k1
  106 format(' max(rotB_y-jy)=',e12.4,3i4)
*/
    s = 0.0;
    maxind = {0, 0, 0};
    for (int k = 1; k < km + 1; k++) {
        for (int l = 0; l < lm + 1; l++) {
            for (int i = 1; i < im + 1; i++) {
                const double s1 = (bx(i,l,k) - bx(i,l,k-1)) * rhz +
                                  (bz(i,l,k) - bz(i-1,l,k)) * rhx - jy(i,l,k);
                if (std::abs(s1) > s) {
                    s = s1;
                    maxind = {i, l, k};
                }
            }
        }
    }

/*
      do k=1,km+1
         do l=2,lm+1
            do i=2,im+1
               s4=(by(i,l,k)-by(i-1,l,k))/hx-
     =            (bx(i,l,k)-bx(i,l-1,k))/hy-jz(i,l,k)
               if(dabs(s4).gt.s3) then
                  s3=s4
                  i1=i
                  l1=l
                  k1=k
               endif
            enddo
         enddo
      enddo
      write(25,107) s3,i1,l1,k1
  107 format(' max(rotB_z-jz)=',e12.4,3i4)
*/
    s = 0.0;
    maxind = {0, 0, 0};
    for (int k = 0; k < km + 1; k++) {
        for (int l = 1; l < lm + 1; l++) {
            for (int i = 1; i < im + 1; i++) {
                const double s1 = (by(i,l,k) - by(i-1,l,k)) * rhx +
                                  (bx(i,l,k) - by(i,l-1,k)) * rhy - jz(i,l,k);
                if (std::abs(s1) > s) {
                    s = s1;
                    maxind = {i, l, k};
                }
            }
        }
    }

    outputDat("jx.dat", jx, im + 1, lm + 2, km + 2);
    outputDat("jy.dat", jy, im + 2, lm + 1, km + 2);
    outputDat("jz.dat", jz, im + 2, lm + 2, km + 1);
    outputDat("ax.dat", ax, im + 1, lm + 2, km + 2);
    outputDat("ay.dat", ay, im + 2, lm + 1, km + 2);
    outputDat("az.dat", az, im + 2, lm + 2, km + 1);
    outputDat("bx.dat", bx, im + 2, lm + 1, km + 1);
    outputDat("by.dat", by, im + 1, lm + 2, km + 1);
    outputDat("bz.dat", bz, im + 1, lm + 1, km + 2);

    return 0;
}
