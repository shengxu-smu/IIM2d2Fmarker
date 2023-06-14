c-----------------------------------------------------------------------
c
      subroutine chorin_iterate(krk,fac,iflag)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,iflag
      real*8 fac,alamda,divmax,divlimiter

      do j=1,ny
        do i=1,nx
          d(i,j)=(1.0d0/dx)*(u(i,j)-u(i-1,j))+
     .           (1.0d0/dy)*(v(i,j)-v(i,j-1))
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

      divmax=0.0
      divlimiter=1.0d-4
      do j=1,ny
        do i=1,nx
          divmax=max(divmax,abs(d(i,j)))
        enddo
      enddo
      
      iflag=0
      if(divmax.ge.divlimiter) then
        iflag=1
        alamda=0.9d0/((fac*dt)*(dx2+dy2))
        do j=1,ny
          do i=1,nx
            p(i,j)=p(i,j)-alamda*d(i,j)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
