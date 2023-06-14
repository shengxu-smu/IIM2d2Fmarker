c-----------------------------------------------------------------------
c
      subroutine pgrad(krk)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer krk
      real*8 grad
      real*8 funrho, rr

      do j=1,ny
        do i=1,nx-1
          grad=-dx1*(p(i+1,j)-p(i,j))
          urk(i,j,krk)=urk(i,j,krk)+grad
        enddo
      enddo

      !---------------------------!
      !            v              !
      !---------------------------!

      do j=1,ny-1
        do i=1,nx
          grad=-dy1*(p(i,j+1)-p(i,j))
          vrk(i,j,krk)=vrk(i,j,krk)+grad
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
