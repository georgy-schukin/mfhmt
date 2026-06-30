      program simit44
      implicit none
      integer im,lm,km,nm,i,l,k,n,m,i1,l1,k1,i2,l2,k2
      parameter(im=40,lm=40,km=20)

      real*8 jx(im+2,lm+2,km+2),jy(im+2,lm+2,km+2),jz(im+2,lm+2,km+2),
     =ax(im+2,lm+2,km+2),ay(im+2,lm+2,km+2),az(im+2,lm+2,km+2),
     =bx(im+2,lm+2,km+2),by(im+2,lm+2,km+2),bz(im+2,lm+2,km+2)

      real*8 qj,pi,c1,xm,ym,zm,x0,y0,hx,hy,hz,r0,h0,eps,
     =x,y,z,x1,y1,z1,x2,y2,z2,s2,s4,s6,s,c2,s1,s3,r,c12,c13,c21,c23,
     =sx,sy,sz


      common/j/jx,jy,jz

      open(25,file='smt45.lst',form='formatted')

      qj=10.d0                    ! ?
      pi=3.14159265358979d0
      nm=1000
      c1=2.d0*pi/nm
      xm=4.d0
      ym=4.d0
      zm=2.d0
      x0=2.d0
      y0=2.d0
c      im=40
c      lm=40
c      km=20
      hx=xm/im
      hy=ym/lm
      hz=zm/km
      
      r0=1.d0
      h0=zm/(2.d0*pi)
      write(25,*)'qj=',qj

      print*,'1. im,lm,km,nm=',im,lm,km,nm
      write(25,100) im,lm,km,nm
  100 format('1. im,lm,km,nm=',4i6)
      write(25,102) hx,hy,hz
  102 format('1. hx,hy,hz=',3f10.3)


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

      eps=1.d-10
      c2=2.d0/hx**2+2.d0/hy**2+2.d0/hz**2

c=========================================================
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

      do k=2,km+1
         do l=2,lm+1
            ax(1,l,k)=ax(2,l,k)+c12*(ay(2,l,k)-ay(2,l-1,k))+
     =               c13*(az(2,l,k)-az(2,l,k-1))
            ax(im+1,l,k)=ax(im,l,k)-c12*(ay(im+1,l,k)-ay(im+1,l-1,k))-
     =               c13*(az(im+1,l,k)-az(im+1,l,k-1))
         enddo
      enddo

      do l=1,lm+2
         do i=1,im+1
            ax(i,l,1)=ax(i,l,km+1)
            ax(i,l,km+2)=ax(i,l,2)
         enddo
      enddo

      do k=1,km+2
         do i=1,im+1
            ax(i,1,k)=ax(i,2,k)             ! ?
            ax(i,lm+2,k)=ax(i,lm+1,k)       ! ?
         enddo
      enddo

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

      do k=2,km+1
         do i=2,im+1
            ay(i,1,k)=ay(i,2,k)+c21*(ax(i,2,k)-ax(i-1,2,k))+
     =          c23*(az(i,2,k)-az(i,2,k-1))
            ay(i,lm+1,k)=ay(i,lm,k)-c21*(ax(i,lm+1,k)-ax(i-1,lm+1,k))-
     =          c23*(az(i,lm+1,k)-az(i,lm+1,k-1))
         enddo
      enddo

      do l=1,lm+1
         do i=1,im+2
            ay(i,l,1)=ay(i,l,km+1)
            ay(i,l,km+2)=ay(i,l,2)
         enddo
      enddo

      do k=1,km+2
         do l=1,lm+1
            ay(1,l,k)=ay(2,l,k)        ! ?
            ay(im+2,l,k)=ay(im+1,l,k)  ! ?
         enddo
      enddo

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

      do l=1,lm+2
         do i=1,im+2
            az(i,l,1)=az(i,l,km+1)
            az(i,l,km+2)=az(i,l,2)
         enddo
      enddo

c      write(25,*) 'ax. n,sx=',n,sx
      print 104,n,sx,sy,sz
  104 format('n,sx,sy,sz=',i6,3e12.4)

      if((sx.gt.eps).or.(sy.gt.eps).or.(sz.gt.eps)) goto 8
c      if(n.lt.200) goto 8

      write(25,104) n,sx,sy,sz
c========================= bx, by, bz ============

      do k=1,km+1
         do l=1,lm+1
            do i=1,im+2
               bx(i,l,k)=(az(i,l+1,k)-az(i,l,k))/hy-
     =                   (ay(i,l,k+1)-ay(i,l,k))/hz
            enddo
         enddo
      enddo

      do k=1,km+1
         do l=1,lm+2
            do i=1,im+1
               by(i,l,k)=(ax(i,l,k+1)-ax(i,l,k))/hz-
     =                   (az(i+1,l,k)-az(i,l,k))/hx
            enddo
         enddo
      enddo

      do k=1,km+2
         do l=1,lm+1
            do i=1,im+1
               bz(i,l,k)=(ay(i+1,l,k)-ay(i,l,k))/hx-
     =                   (ax(i,l+1,k)-ax(i,l,k))/hy
            enddo
         enddo
      enddo

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

c===================================== write

      open(17,file='jx.dat',form='formatted')
      do k=1,km+2
         do l=1,lm+2
            do i=1,im+1
               write(17,101) i,l,k,jx(i,l,k)
  101          format(3i5,es12.4)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='jy.dat',form='formatted')
      do k=1,km+2
         do l=1,lm+1
            do i=1,im+2
               write(17,101) i,l,k,jy(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='jz.dat',form='formatted')
      do k=1,km+1
         do l=1,lm+2
            do i=1,im+2
               write(17,101) i,l,k,jz(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='ax.dat',form='formatted')
      do k=1,km+2
         do l=1,lm+2
            do i=1,im+1
               write(17,101) i,l,k,ax(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='ay.dat',form='formatted')
      do k=1,km+2
         do l=1,lm+1
            do i=1,im+2
               write(17,101) i,l,k,ay(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='az.dat',form='formatted')
      do k=1,km+1
         do l=1,lm+2
            do i=1,im+2
               write(17,101) i,l,k,az(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='bx.dat',form='formatted')
      do k=1,km+1
         do l=1,lm+1
            do i=1,im+2
               write(17,101) i,l,k,bx(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='by.dat',form='formatted')
      do k=1,km+1
         do l=1,lm+2
            do i=1,im+1
               write(17,101) i,l,k,by(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      open(17,file='bz.dat',form='formatted')
      do k=1,km+2
         do l=1,lm+1
            do i=1,im+1
               write(17,101) i,l,k,bz(i,l,k)
            enddo
         enddo
      enddo
      close(17)

      end
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