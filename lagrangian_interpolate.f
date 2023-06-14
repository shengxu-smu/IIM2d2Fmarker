c-----------------------------------------------------------------------
c
      subroutine lagrangian_interpolate(nc,nn,abxy,fabxy,ccxy,tmp)
      include 'parameter.inc'
      include 'surface.inc'
      integer nc,nn,ip
      real*8 bl,bu
      real*8 abxy(0:nn),fabxy(2,0:nn),ccxy(2,0:2*nn)
      real*8 ffm(2),ffmp(2),ccm(2),ccmp(2),foo(2),dfoo(2)
      real*8 tmp(2,0:ns)

      if(abxy(0).gt.0.0d0) then
        ip=-1
      else
        ip=0
      endif
      do m=0,ns-1
5       continue
        if(ip.eq.-1) then
          bl=abxy(nc-1)-2.0d0*pi
          if(abs(bl).lt.1.0d-5) bl=0.0d0
          bu=abxy(0)
          do n=1,2
            ffm(n)=fabxy(n,nc-1)
            ffmp(n)=fabxy(n,nc)
            ccm(n)=ccxy(n,nc-1)
            ccmp(n)=ccxy(n,nc)
          enddo
        else
          bl=abxy(ip)
          bu=abxy(ip+1)
          do n=1,2
            ffm(n)=fabxy(n,ip)
            ffmp(n)=fabxy(n,ip+1)
            ccm(n)=ccxy(n,ip)
            ccmp(n)=ccxy(n,ip+1)
          enddo
        endif
        if(alfa(m).ge.bl.and.alfa(m).lt.bu) then
          call fdf(bl,bu,ffm,ffmp,
     .             ccm,ccmp,alfa(m),foo,dfoo,2)
          tmp(1,m)=foo(1)
          tmp(2,m)=foo(2)
        else
          ip=ip+1
          goto 5
        endif
      enddo
      tmp(1,ns)=tmp(1,0)
      tmp(2,ns)=tmp(2,0)


      return
      end


c-----------------------------------------------------------------------
