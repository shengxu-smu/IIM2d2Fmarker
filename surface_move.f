c-----------------------------------------------------------------------  
c
      subroutine surface_move(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      integer krk
      real*8 fac

      if(krk.eq.1) t0=t0+0.5d0*dt
      if(krk.eq.3) t0=t0+0.5d0*dt

c      theta(1)=(pi-0.25d0*pi*(1.0d0-dsin(0.8d0*t0)))*0.0d0
c      xsc(1)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dcos(pi/3.0d0)*0.0d0
c      ysc(1)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dsin(pi/3.0d0)*0.0d0

c      thetat(1)=0.25d0*pi*0.8d0*dcos(0.8d0*t0)*0.0d0
c      xsct(1)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dcos(pi/3.0d0)*0.0d0
c      ysct(1)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dsin(pi/3.0d0)*0.0d0

c      thetatt(1)=-0.25d0*pi*0.64d0*dsin(0.8d0*t0)*0.0d0
c      xsctt(1)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dcos(pi/3.0d0)*0.0d0
c      ysctt(1)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dsin(pi/3.0d0)*0.0d0

      theta(1)=0
      thetat(1)=0.0d0

c      theta(2)=(pi-0.25d0*pi*(1.0d0-dsin(0.8d0*t0)))

c      thetat(2)=0.25d0*pi*0.8d0*dcos(0.8d0*t0)
c      xsct(2)=0.0d0
c      ysct(2)=0.0d0

c      thetatt(2)=-0.25d0*pi*0.64d0*dsin(0.8d0*t0)
c      xsctt(2)=0.0d0
c      ysctt(2)=0.0d0

c      theta(3)=0.0d0

c      thetat(3)=0.0d0
c      xsct(3)=0.0d0
c      ysct(3)=0.0d0

c      thetatt(3)=0.0d0
c      xsctt(3)=0.0d0
c      ysctt(3)=0.0d0

      DO l=1,ms

      do m=0,ns-1
        xs(m,l)=xsc(l)+xs0(m,l)*dcos(theta(l))-ys0(m,l)*dsin(theta(l))
        ys(m,l)=ysc(l)+xs0(m,l)*dsin(theta(l))+ys0(m,l)*dcos(theta(l))
        us(m,l)=xsct(l)-thetat(l)*(ys(m,l)-ysc(l))
        vs(m,l)=ysct(l)+thetat(l)*(xs(m,l)-xsc(l))
      enddo
      xs(ns,l)=xs(0,l)
      ys(ns,l)=ys(0,l)
      us(ns,l)=us(0,l)
      vs(ns,l)=vs(0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
