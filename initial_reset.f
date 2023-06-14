c-----------------------------------------------------------------------
c
      subroutine initial_reset
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer nn
      real*8 xx,yy,ang,r,rc,ds

      ds=0.0d0*dsqrt(dx*dx+dy*dy)

      DO l=1,ms

      do j=1,ny
        do i=1,nx-1
          xx=xe(i)-xsc(l)
          yy=yc(j)-ysc(l)
          ang=datan2(yy,xx)+pi-theta(l)
          if(ang.lt.0.0d0) ang=ang+2.0d0*pi
          r=dsqrt(xx*xx+yy*yy)
          nn=int(ang/dar)
          rc=ra(nn,l)
          if(r.le.rc+ds) then
            u(i,j)=xsct(l)-thetat(l)*yy
          endif
        enddo
      enddo

      do j=1,ny-1
        do i=1,nx
          xx=xc(i)-xsc(l)
          yy=ye(j)-ysc(l)
          ang=datan2(yy,xx)+pi-theta(l)
          if(ang.lt.0.0d0) ang=ang+2.0d0*pi
          r=dsqrt(xx*xx+yy*yy)
          nn=int(ang/dar)
          rc=ra(nn,l)
          if(r.le.rc+ds) then
            v(i,j)=ysct(l)+thetat(l)*xx
          endif
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------

