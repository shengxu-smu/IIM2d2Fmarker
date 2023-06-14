c-----------------------------------------------------------------------
c
      subroutine IterpressureSolver
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer iter
      real*8 fac,err,tol,gama

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !               ---  Initial Guess  ---                 !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
 
      do j=1,ny
        do i=1,nx
          p(i,j)=0.0d0
        enddo
      enddo

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                ---  Iterations  ---                   !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      iter=0
      gama=1.0d-2
      tol=dx*dx*dx
      err=1.0d6
 9    continue
      if(err.gt.tol.and.iter.lt.200) then
        call jc_dpdn
        call diffusion_pressure(gama)
        err=0.0d0
        do j=1,ny
          do i=1,nx
            err=max(err,abs(p(i,j)-o(i,j)))
          enddo
        enddo
        iter=iter+1
c        write(*,*)' iter = ', iter,' err = ', err
        goto 9
      endif
      write(*,*),'New solver iterations ',iter, ' err = ',err

      return
      end


c-----------------------------------------------------------------------
