c-----------------------------------------------------------------------
c
      subroutine initial
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer ni,nj
      real*8 o1,o2,r1,r2,aa,bb,rr,xcr,ycr
      real*8 unone,vnone,pnone

      call vrffti(ns,wsave)
      call surface_parametrization

      do l=1,ms
        do m=0,ns
          us(m,l)=xsct(l)-thetat(l)*(ys(m,l)-ysc(l))
          vs(m,l)=ysct(l)+thetat(l)*(xs(m,l)-xsc(l))
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx
            xcr=xe(i)
            ycr=yc(j)
            u(i,j)=0.0d0
            !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
            !       ---   Exact solution u  ---          !
            !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
c            call ExactSol(xcr,ycr,u(i,j),vnone,pnone)
            !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
        enddo
      enddo

      do j=0,ny
        do i=0,nx+1
            xcr=xc(i)
            ycr=ye(j)
            v(i,j)=0.0d0
            !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
            !       ---  Exact Solution v ---             !
            !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
c            call ExactSol(xcr,ycr,unone,v(i,j),pnone)
            !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          p(i,j)=0.0d0
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          d(i,j)=0.0d0
          dn(i,j)=0.0d0
        enddo
      enddo

      if(lw_pbc.eq.0.and.le_pbc.eq.0.and.
     .   ls_pbc.eq.0.and.ln_pbc.eq.0) then
        ni=nx-1
        nj=ny-1
        call vrffti(ni,wsavei)
        call vrffti(nj,wsavej)
      endif

      if(lw_pbc.eq.1.and.le_pbc.eq.1.and.
     .   ls_pbc.eq.1.and.ln_pbc.eq.1) then
        ni=nx-2
        nj=ny-2
        call vsinti(ni,wsavei)
        call vsinti(nj,wsavej)
      endif

      if(lw_pbc.eq.2.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.2.and.ln_pbc.eq.2) then
        ni=nx
        nj=ny
        call vcosti(ni,wsavei)
        call vcosti(nj,wsavej)
      endif

      if(lw_pbc.eq.1.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.2.and.ln_pbc.eq.2) then
        ni=nx-1
        nj=ny
        call vsinqi(ni,wsavei)
        call vcosti(nj,wsavej)
      endif      

      if(lw_pbc.eq.1.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.0.and.ln_pbc.eq.0) then
        ni=nx-1
        nj=ny-1
        call vsinqi(ni,wsavei)
        call vrffti(nj,wsavej)
      endif

      return
      end


c-----------------------------------------------------------------------
