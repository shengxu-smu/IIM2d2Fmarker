c-----------------------------------------------------------------------
c
      subroutine difference_winded(l,m,n,x,f,df)
      integer l,m,n,i,j,kp
      real*8 x(l),f(l),r(l),c(l),h(l-1),g(l-1),s(l-1,l-1)
      real*8 dx,df,x0,g0,factorial

      x0=x(m)
      g0=f(m)
      kp=1
      do i=1,l
       if(i.ne.m) then
         dx=x(i)-x0
         if(dx.ne.0.0d0) then
           h(kp)=dx
           g(kp)=f(i)
           kp=kp+1
         elseif(n.eq.0) then
           df=f(i)
           return
         endif
       endif
      enddo
      kp=kp-1

      do j=1,kp-1
        do i=j+1,kp
          s(i,j)=h(kp-i+j+1)-h(kp-i+1)
        enddo
      enddo

      do j=1,kp
        s(1,j)=h(j)
      enddo 
      do i=2,kp
        do j=i,kp
          s(i,j)=s(kp-i+2,j-i+1)*s(i-1,j)
        enddo
      enddo

      do i=1,kp
        r(i)=0.0d0
        if(i.eq.n) r(i)=factorial(n)
      enddo
      r(kp+1)=0.0d0
      if(n.eq.0) r(kp+1)=1.0d0

      do j=1,kp-1
        do i=kp,j+1,-1
          r(i)=r(i)-h(j)*r(i-1)
        enddo
      enddo

      c(kp)=r(kp)/s(kp,kp)
      do i=kp-1,1,-1
        dx=0.0d0
        do j=i+1,kp
          dx=dx+s(i,j)*c(j)
        enddo
        c(i)=(r(i)-dx)/s(i,i)
      enddo
      dx=0.0d0
      do i=1,kp
        dx=dx-c(i)
      enddo
      c(kp+1)=dx

      df=c(kp+1)*g0
      do i=1,kp
        df=df+c(i)*g(i)
      enddo


      return
      end


c-----------------------------------------------------------------------
c

      function factorial(n)
      integer n,nf,i
      real*8 factorial

      nf=1
      if(n.gt.0) then
        do i=1,n
          nf=nf*i
        enddo
      endif
      factorial=dble(nf)

      return
      end

c-----------------------------------------------------------------------
