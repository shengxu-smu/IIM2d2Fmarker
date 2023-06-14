c-----------------------------------------------------------------------
c
      subroutine data_plot
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'
      integer ih,jh,iref,igap,iend,jgap,jend,ii,jj,mark1,mark2,iexact
      real*8 foo,unorm,vnorm,pnorm1,pnorm2,cp1,cp2,xcr,ycr,d0,pmax,pmin
      real*8 au(0:nx,0:ny+1),av(0:nx+1,0:ny)
      real*8 o1,o2,r1,r2,aa,bb,rr
      real*8 tao_x, tao_y, nor_x, nor_y, norm_tao, Delta_n
      real*8 foo1
      real*8 xn,yn,gacobi,ExactJdpdn(0:ns,ms),ApproxJdpdn(0:ns,ms)
      real*8 coE,coI,rho
      real*8 unone,vnone,pnone,Ep1,Ep2,Ep3,Ep4,Epmid
      integer mid

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      !     ---  Data files: xe, ye, xc, yc  ---         !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      call field_interpolate
      call streamfunction

200   format(1x,10e16.6)

      open(unit=24,file='DAT/xe.dat',status='unknown')
      open(unit=34,file='DAT/ye.dat',status='unknown')
      do i=1,nx
        write(24,200)xe(i)
      enddo
      do j=1,ny
        write(34,200)ye(j)
      enddo
      close(24)
      close(34)
      open(unit=44,file='DAT/xc.dat',status='unknown')
      open(unit=54,file='DAT/yc.dat',status='unknown')
      do i=1,nx
        write(44,200)xc(i)
      enddo
      do j=1,ny
        write(54,200)yc(j)
      enddo
      close(44)
      close(54)


      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      !  ---  Printing the approx together: gx, gy ---   !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      ih=(0.0d0-x0)/dx
      jh=(0.0d0-y0)/dy
      open(unit=14,file='DAT/gy.dat',status='unknown')
      open(unit=24,file='DAT/gx.dat',status='unknown')
      do j=1,ny 
        write(14,200)yc(j),ucc(ih,j),vcc(ih,j),p(ih,j),o(ih,j)
      enddo
      do i=1,nx 
        write(24,200)xc(i),ucc(i,jh),vcc(i,jh),p(i,jh),o(i,jh)
      enddo
      close(14)
      close(24)


      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      !  ---  Printing the approximations separtely ---  !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

400   format(1x,2000e16.6)

      open(unit=38,file='DAT/p.dat',status='unknown')
      open(unit=48,file='DAT/d.dat',status='unknown')
      open(unit=58,file='DAT/uc.dat',status='unknown')
      open(unit=68,file='DAT/vc.dat',status='unknown')
      open(unit=78,file='DAT/wo.dat',status='unknown')
      open(unit=88,file='DAT/ph.dat',status='unknown')
      do j=1,ny
        do i=1,nx
          rho = rho_neg*funHp(i,j) + rho_pos*(1.0d0-funHp(i,j))
          p(i,j)=rho*p(i,j)
        enddo
      enddo
      do j=1,ny
        write(38,400)(p(i,j),i=1,nx)
        write(48,400)(d(i,j),i=1,nx)
        write(58,400)(ucc(i,j),i=1,nx)
        write(68,400)(vcc(i,j),i=1,nx)
        write(78,400)(o(i,j),i=1,nx)
        write(88,400)(-ph(i,j),i=1,nx)
      enddo
      close(38)
      close(48)
      close(58)
      close(68)
      close(78)
      close(88)

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      ! - Calculating and printing the exact solution -  !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      iexact=1
      if(iexact.eq.1) then
        
        ! Looking the constants of the pressure

        mid = int(nx/2.0d0)
        call ExactSol(xc(1),yc(1),unone,vnone,Ep1)
        call ExactSol(xc(mid),yc(mid),unone,vnone,Epmid)

        coE = Ep1-p(1,1)
        coI = Epmid-p(mid,mid)

        !  Calculating the exact solutions

        do j=1,ny
          do i=1,nx
            xcr=xc(i)
            ycr=yc(j)
            rr = dsqrt(xcr**2+ycr**2)
            ap(i,j)= (0.5d0*(rr**2)*rho_neg-coI)*funHp(i,j)+
     .               (0.5d0*(rr**2)*rho_pos-coE)*(1-funHp(i,j))
c            if (rr>0.5d0-0.05d0) then 
c               ap(i,j)= (0.5d0*(rr**2)*rho_pos-coE)
c            else
c               ap(i,j)= (0.5d0*(rr**2)*rho_neg-coI)
c            endif

            call ExactSol(xcr,ycr,au(i,j),av(i,j),pnone)
          enddo
        enddo
      
      ! Printing the analytical solution

        open(unit=44,file='DAT/ap.dat',status='unknown')
        open(unit=54,file='DAT/au.dat',status='unknown')
        open(unit=64,file='DAT/av.dat',status='unknown')

        do j=1,ny
          write(54,400)(au(i,j),i=1,nx)
          write(64,400)(av(i,j),i=1,nx)  
          write(44,400)(ap(i,j),i=1,nx)     
        enddo

        close(44)
        close(54) 
        close(64)

      endif
                       
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      !               ---  Errors ---                   !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      unorm = 0.0d0
      vnorm = 0.0d0
      pnorm1= 0.0d0
      pnorm2= 0.0d0
      pmax  = 0.0d0
      pmin  = 0.0d0
      do j=1,ny
        do i=1,nx
          if(iexact.eq.1) then
            ucc(i,j)= ucc(i,j)- au(i,j)
            vcc(i,j)= vcc(i,j)- av(i,j)
            p(i,j)  = p(i,j)  - ap(i,j)
          endif
          unorm = max(unorm,abs(ucc(i,j)))
          vnorm = max(vnorm,abs(vcc(i,j)))
          pnorm1= max(pnorm1,abs(p(i,j)))
          pmax  = max(pmax,p(i,j))
          pmin  = min(pmin,p(i,j))
        enddo
      enddo
      pnorm2=0.5d0*(pmax-pmin)

      write(*,*)
      write(*,*)'unorm  = ',unorm
      write(*,*)'vnorm  = ',vnorm
      write(*,*)'pnorm1 = ',pnorm1
      write(*,*)'pnorm2 = ',pnorm2

      ! Printing Errors
      open(unit=14,file='DAT/eu.dat',status='unknown')
      open(unit=24,file='DAT/ev.dat',status='unknown')
      open(unit=34,file='DAT/ep.dat',status='unknown')

      do j=1,ny
        write(14,400)(abs(ucc(i,j)),i=1,nx)
        write(24,400)(abs(vcc(i,j)),i=1,nx)  
        write(34,400)(abs(p(i,j)),i=1,nx)     
      enddo      

      close(14)
      close(24)
      close(34)

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      !               ---   xs,ys ---                   !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      open(unit=54,file='DAT/xs.dat',status='unknown')
      open(unit=64,file='DAT/ys.dat',status='unknown')

      DO l=1,ms
        do m=0,ns
           write(54,400),xs(m,l)
           write(64,400),ys(m,l)
        enddo
      ENDDO

      close(54)
      close(64)

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      !       ---   Printing fy0=[dp/dn] ---            !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      open(unit=54,file='DAT/Jumps/Ejc_dpdn.dat',status='unknown')
      open(unit=64,file='DAT/Jumps/Njc_dpdn.dat',status='unknown')

      DO l=1,ms
        do m=0,ns-1
           write(54,400),0.5d0*(rho_pos-rho_neg)
           write(64,400),fy0(m,l)
        enddo
      ENDDO
      close(54)
      close(64)

      !=====================================!

500   format(1x,e36.18)

      return
      end


c-----------------------------------------------------------------------
