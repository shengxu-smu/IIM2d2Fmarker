c-----------------------------------------------------------------------
c
      subroutine difference(m,n,x,f,df)
      integer m,n,i,ip
      real*8 x(5),f(5),h(4),g(4)
      real*8 df,x0,g0,r1,r2,r3,r4,d21,d31,d32,d41,d42,d43
      real*8 hh2,hh3,hh4,hhh3,hhh4,hhhh4,c0,c1,c2,c3,c4

      r1=0.0d0
      r2=0.0d0
      if(n.eq.1) r1=1.0d0
      if(n.eq.2) r2=2.0d0  

      x0=x(m)
      g0=f(m)
      ip=1
      do i=1,5
       if(i.ne.m) then
         h(ip)=x(i)-x0
         g(ip)=f(i)
         ip=ip+1
       endif
      enddo
      d21=h(2)-h(1)
      d31=h(3)-h(1)
      d32=h(3)-h(2)
      d41=h(4)-h(1)
      d42=h(4)-h(2)
      d43=h(4)-h(3)

      hh2=h(2)*d21
      hh3=h(3)*d31
      hh4=h(4)*d41
      hhh3=hh3*d32
      hhh4=hh4*d42
      hhhh4=hhh4*d43

      r4=0.0d0
      r3=-h(1)*r2
      r2=r2-h(1)*r1

      r4=r4-h(2)*r3
      r3=r3-h(2)*r2

      r4=r4-h(3)*r3

      c4=r4/hhhh4
      c3=(r3-hhh4*c4)/hhh3
      c2=(r2-hh4*c4-hh3*c3)/hh2
      c1=(r1-h(4)*c4-h(3)*c3-h(2)*c2)/h(1)
      c0=0.0d0-c4-c3-c2-c1

      df=c0*g0+c1*g(1)+c2*g(2)+c3*g(3)+c4*g(4)


      return
      end

c-----------------------------------------------------------------------
