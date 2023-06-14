c-----------------------------------------------------------------------
c
      subroutine pressure(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'old.inc'
      integer krk
      real*8 fac,pw(ny),pe(ny),ps(nx),pn(nx)
      real*8 funrho
      real*8 hola(1:nx,1:ny)
      real*8 Jdpdn(0:ns-1,ms)
      real*8 xt, yt, gacobi
      real*8 Rep
      real*8 nada

      do j=1,ny
        do i=1,nx
          d(i,j)=dx1*(u(i,j)-u(i-1,j))+dy1*(v(i,j)-v(i,j-1))
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            d(ixc(i,l),jexc(i,l)+1)=d(ixc(i,l),jexc(i,l)+1)+pcvdy(i,l)
          enddo
          do j=0,ncyc(l)
            d(ieyc(j,l)+1,jyc(j,l))=d(ieyc(j,l)+1,jyc(j,l))+pcudx(j,l)
          enddo
        enddo
      endif
      
      !------------------------------!
      !    ---  rhs pressure ---     !
      !------------------------------!

      do j=1,ny
        do i=1,nx
          Rep = Re_neg*funHp(i,j)+Re_pos*(1.0d0-funHp(i,j))
          prhs(i,j)= (
     .             dn(i,j)/(fac*dt)
     .            +(dx2/Rep)*(d(i+1,j)-2.0d0*d(i,j)+d(i-1,j))
     .            +(dy2/Rep)*(d(i,j+1)-2.0d0*d(i,j)+d(i,j-1))
     .            -2.0d0*(d(i,j)*d(i,j)+
     .             (ucc(i,j)/dx)*0.5d0*(d(i+1,j)-d(i-1,j))+
     .             (vcc(i,j)/dy)*0.5d0*(d(i,j+1)-d(i,j-1)))
     .              -2.0d0*(dvx(i,j)*duy(i,j)-dux(i,j)*dvy(i,j)))
        enddo
      enddo
c      call prhs_reset
 
      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            prhs(ixc(i,l),jcxc(i,l))=prhs(ixc(i,l),jcxc(i,l))-
     .                               pcdyy(1,i,l)
            prhs(ixc(i,l),jcxc(i,l)+1)=prhs(ixc(i,l),jcxc(i,l)+1)-
     .                                 pcdyy(2,i,l)
          enddo
          do j=0,ncyc(l)
            prhs(icyc(j,l),jyc(j,l))=prhs(icyc(j,l),jyc(j,l))-
     .                               pcdxx(1,j,l)
            prhs(icyc(j,l)+1,jyc(j,l))=prhs(icyc(j,l)+1,jyc(j,l))-
     .                                 pcdxx(2,j,l)
          enddo
        enddo
      endif

      if(ipfft.eq.1) then
        call pbc_fft(pw,pe,ps,pn)
        call poisson_fft(pw,pe,ps,pn,
     .                   lw_pbc,le_pbc,ls_pbc,ln_pbc,1)
      endif

      do j=1,ny
        do i=1,nx
          o(i,j)=p(i,j)
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
