c-----------------------------------------------------------------------
c
      subroutine ubc(krk,fac)
      include 'parameter.inc'
      include 'field.inc'
      include 'old.inc'
      integer krk
      real*8 fac,conv,visc,uip,uim,vip,vim
      real*8 o1,o2,r1,r2,a,b,foo,rr
      real*8 Om2
      integer ProblemBubble

      Om2 = 2.0d0
      o1=1.0d0
      o2=-1.0d0
      r1=0.5d0
      r2=2.0d0
      a = 1.0d0  
      b = 0.0d0  

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                 ---  PERIODIC ---                     !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      if(ls_ubc.eq.0.and.ln_ubc.eq.0) then
        do i=0,nx
          u(i,0)=u(i,ny-1)
          u(i,ny+1)=u(i,2)
        enddo
      endif

      if(lw_ubc.eq.0.and.le_ubc.eq.0) then
        do j=0,ny+1
          u(0,j)=u(nx-1,j)
          u(nx,j)=u(1,j)
        enddo
      endif

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                 ---  DIRICHLET  ---                   !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      ProblemBubble=1
      IF (ProblemBubble==1) THEN
        do j=0,ny+1
          foo=0.0d0 
          u(0,j) =2.0d0*foo-u(1,j)
          u(nx,j)=2.0d0*foo-u(nx-1,j)
        enddo       
        do i=0,nx
          foo=0.0d0
          u(i,0)   =2.0d0*foo-u(i,2)
          u(i,ny+1)=2.0d0*foo-u(i,ny-1)
        enddo

      ELSE

      if(lw_ubc.eq.1) then
        do j=0,ny+1
          rr=xc(1)*xc(1)+yc(j)*yc(j)
          foo=-(a+b/rr)*yc(j)
          u(0,j)=-a*yc(j)
c          u(1,j)=-a*yc(j)
c          u(0,j)=2.0d0*foo-u(1,j)
        enddo
      endif

      if(le_ubc.eq.1) then
        do j=0,ny+1
          rr=xc(nx)*xc(nx)+yc(j)*yc(j)
          foo=-(a+b/rr)*yc(j)
c          u(nx-1,j)=-a*yc(j)
          u(nx,j)  =-a*yc(j)
c          u(nx,j)=2.0d0*foo-u(nx-1,j)
        enddo
      endif

      if(ls_ubc.eq.1) then
        do i=0,nx
          rr=xe(i)*xe(i)+yc(1)*yc(1)
          foo=-(a+b/rr)*yc(1)
c          u(i,1)=-a*yc(1) !added
          u(i,0)=-a*yc(0) !added
c          u(i,0)=2.0d0*foo-u(i,2)
        enddo
      endif

      if(ln_ubc.eq.1) then
        do i=0,nx
          rr=xe(i)*xe(i)+yc(ny)*yc(ny)
          foo=-(a+b/rr)*yc(ny)
c          u(i,ny)=-a*yc(ny)    !added
          u(i,ny+1)=-a*yc(ny+1)!added
c          u(i,ny+1)=2.0d0*foo-u(i,ny-1)
        enddo
      endif
     
      ENDIF
      
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                  ---  NEWMANN  ---                    !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      if(lw_ubc.eq.2) then
        do j=0,ny+1
          u(0,j)=u(1,j)-dx*foo
        enddo
      endif

      if(le_ubc.eq.2) then
        do j=0,ny+1
          u(nx,j)=u(nx-1,j)+dx*foo
        enddo
      endif

      if(ls_ubc.eq.2) then
        do i=0,nx
          u(i,0)=u(i,2)
        enddo
      endif

      if(ln_ubc.eq.2) then
        do i=0,nx
          u(i,ny+1)=u(i,ny-1)
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
