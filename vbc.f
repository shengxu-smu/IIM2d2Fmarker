c-----------------------------------------------------------------------
c
      subroutine vbc(krk,fac)
      include 'parameter.inc'
      include 'field.inc'
      integer krk
      real*8 fac
      real*8 a,b,foo,rr
      integer ProblemBubble

      a = 1.0d0 
      b = 0.0d0 

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                 ---  PERIODIC ---                     !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      if(lw_vbc.eq.0.and.le_vbc.eq.0) then
        do j=0,ny
          v(0,j)=v(nx-1,j)
          v(nx+1,j)=v(2,j)
        enddo
      endif

      if(ls_vbc.eq.0.and.ln_vbc.eq.0) then
        do i=0,nx+1
          v(i,0)=v(i,ny-1)
          v(i,ny)=v(i,1)
        enddo
      endif

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                 ---  DIRICHLET ---                    !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      ProblemBubble=1
      IF (ProblemBubble==1) THEN
        do j=0,ny
          foo=0.0d0 
          v(0,j)   =2.0d0*foo-v(2,j)
          v(nx+1,j)=2.0d0*foo-v(nx-1,j)
        enddo       
        do i=0,nx+1
          foo=0.0d0
          v(i,0) =2.0d0*foo-v(i,1)
          v(i,ny)=2.0d0*foo-v(i,ny-1)
        enddo

      ELSE

      if(lw_vbc.eq.1) then
        do j=0,ny
          rr=xc(1)*xc(1)+ye(j)*ye(j)
          foo=(a+b/rr)*xc(1)
          v(0,j)=a*xc(0) !added
c          v(1,j)=a*xc(1) !added
c          v(0,j)=2.0d0*foo-v(2,j)
        enddo
      endif

      if(le_vbc.eq.1) then
        do j=0,ny
          rr=xc(nx)*xc(nx)+ye(j)*ye(j)
          foo=(a+b/rr)*xc(nx)
          v(nx+1,j)=a*xc(nx+1) !added
c          v(nx,j)  =a*xc(nx)   !added
c          v(nx+1,j)=2.0d0*foo-v(nx-1,j)
        enddo
      endif

      if(ls_vbc.eq.1) then
        do i=0,nx+1
          rr=xc(i)*xc(i)+yc(1)*yc(1)
          foo=(a+b/rr)*xc(i)
          v(i,0)=a*xc(i) !added
c          v(i,1)=a*xc(i) !added
c          v(i,0)=2.0d0*foo-v(i,1)
        enddo
      endif

      if(ln_vbc.eq.1) then
        do i=0,nx+1
          rr=xc(i)*xc(i)+yc(ny)*yc(ny)
          foo=(a+b/rr)*xc(i)
c          v(i,ny-1)=a*xc(i) !added
          v(i,ny)  =a*xc(i) !added
c          v(i,ny)=2.0d0*foo-v(i,ny-1)
        enddo
      endif
     
      ENDIF

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                  ---  NEWMANN ---                     !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      if(lw_vbc.eq.2) then
        do j=0,ny
          v(0,j)=v(2,j)
        enddo
      endif

      if(le_vbc.eq.2) then
        do j=0,ny
          v(nx+1,j)=v(nx-1,j)
        enddo
      endif

      if(ls_vbc.eq.2) then
        do i=0,nx+1
          rr=xc(i)*xc(i)+yc(1)*yc(1)
          foo=-2.0d0*b*xc(i)*yc(1)/(rr*rr)
          v(i,0)=v(i,1)-dy*foo
        enddo
      endif

      if(ln_vbc.eq.2) then
        do i=0,nx+1
          rr=xc(i)*xc(i)+yc(ny)*yc(ny)
          foo=-2.0d0*b*xc(i)*yc(ny)/(rr*rr)
          v(i,ny)=v(i,ny-1)+dy*foo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
