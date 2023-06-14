c-----------------------------------------------------------------------
c
      subroutine prhs_reset
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      integer nn,ic,ic1,jc,jc1
      real*8 xx,yy,ang,r,rc,signx,signy,gama

      DO l=1,ms

      do j=1,ny
        do i=1,nx
          xx=xc(i)-xsc(l)
          yy=yc(j)-ysc(l)
          ang=datan2(yy,xx)+pi-theta(l)
          if(ang.lt.0.0d0) ang=ang+2.0d0*pi
          r=dsqrt(xx*xx+yy*yy)
          nn=int(ang/dar)
          rc=ra(nn,l)
          if(r.lt.rc) then
            prhs(i,j)= rho_neg*(
     .                -2.0d0*(dvx(i,j)*duy(i,j)-dux(i,j)*dvy(i,j)))
          endif
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------

