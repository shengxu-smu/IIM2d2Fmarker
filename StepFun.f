c-----------------------------------------------------------------------
c
      subroutine StepFun
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'old.inc'
      integer ie,ic,ic1,je,jc,jc1,krk
      real*8 tmp(1,0:ns),ccc(1,0:ns),foo,fac,grad,xt,yt
      real*8 gacobi,gacobi2,dgacobi,dnx,dny,r1,r2,r3,ts,tm
      real*8 fn,tx,ty,dft,dfn,dtx,dty,ddft,ddfn
      real*8 dpxjc,dpyjc,ddpxjc,ddpyjc
      real*8 alfa0,alfa1,alpha,signx,signy
      real*8 aa(6),bb(6),f(6),sx,sy
      real*8 pw(ny),pe(ny),ps(nx),pn(nx),pero(nx,ny)
      real*8 DisplayError,rr,maxx,dd

c jc_firstsecond

      DO l=1,ms

      do m=0,ns-1
        tx=taox(m,l)
        ty=taoy(m,l)
        dtx=dtaox(m,l)
        dty=dtaoy(m,l)
c  -- for 2-fluid pressure
        fn=-1.0d0      ![p]
        dft=0.0d0      ![dp/dn]
        dfn=0.0d0      ![Lp]
c        ddft=0.0d0
c        ddfn=0.0d0
c  -- end for --

        gacobi=sqrt(tx*tx+ty*ty)
        gacobi2=tx*tx+ty*ty

        dpxjc=0.0d0
        dpyjc=0.0d0
        ddpxjc=0.0d0
        ddpyjc=0.0d0

        pjc(1,m,l)=dpxjc
        pjc(3,m,l)=ddpxjc
        pjc(5,m,l)=fn

        pjc(2,m,l)=dpyjc
        pjc(4,m,l)=ddpyjc
        pjc(6,m,l)=fn
      enddo

      pjc(1,ns,l)=pjc(1,0,l)
      pjc(3,ns,l)=pjc(3,0,l)
      pjc(5,ns,l)=pjc(5,0,l)

      pjc(2,ns,l)=pjc(2,0,l)
      pjc(4,ns,l)=pjc(4,0,l)
      pjc(6,ns,l)=pjc(6,0,l)

      ENDDO

c mac_distribute
      
      DO l=1,ms

      do i=0,ncxc(l) 
        alpha=falfaxc(0,i,l)
        signx=falfaxc(2,i,l)
        signy=falfaxc(3,i,l)

        m=int(alpha/dalfa)
        if(m.eq.-1) m=0
        if(m.eq.ns) m=ns-1

        alfa0=alfa(m)
        alfa1=alfa(m+1)
        do n=1,6
          aa(n)=(pjc(n,m+1,l)-pjc(n,m,l))/dalfa
          bb(n)=(pjc(n,m,l)*alfa1-pjc(n,m+1,l)*alfa0)/dalfa
        enddo

        do n=1,6
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,6
          pjcxc(n,i,l)=f(n)
        enddo

        pjcxc(1,i,l)=signx*pjcxc(1,i,l)
        pjcxc(3,i,l)=signx*pjcxc(3,i,l)
        pjcxc(5,i,l)=signx*pjcxc(5,i,l)

        pjcxc(2,i,l)=signy*pjcxc(2,i,l)
        pjcxc(4,i,l)=signy*pjcxc(4,i,l)
        pjcxc(6,i,l)=signy*pjcxc(6,i,l)

      enddo

      do j=0,ncyc(l)
        alpha=falfayc(0,j,l)
        signx=falfayc(2,j,l)
        signy=falfayc(3,j,l)

        m=int(alpha/dalfa)
        if(m.eq.-1) m=0
        if(m.eq.ns) m=ns-1

        alfa0=alfa(m)
        alfa1=alfa(m+1)
        do n=1,6
          aa(n)=(pjc(n,m+1,l)-pjc(n,m,l))/dalfa
          bb(n)=(pjc(n,m,l)*alfa1-pjc(n,m+1,l)*alfa0)/dalfa
        enddo

        do n=1,6
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,6
          pjcyc(n,j,l)=f(n)
        enddo

        pjcyc(1,j,l)=signx*pjcyc(1,j,l)
        pjcyc(3,j,l)=signx*pjcyc(3,j,l)
        pjcyc(5,j,l)=signx*pjcyc(5,j,l)

        pjcyc(2,j,l)=signy*pjcyc(2,j,l)
        pjcyc(4,j,l)=signy*pjcyc(4,j,l)
        pjcyc(6,j,l)=signy*pjcyc(6,j,l)

      enddo

      ENDDO

c correction_difference

      DO l=1,ms

      do i=0,ncxc(l)
        sy=falfaxc(1,i,l)
        je=jexc(i,l)
        jc=jcxc(i,l)
        jc1=jc+1
        pcdyy(1,i,l)=-dy2*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc1)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
        pcdyy(2,i,l)=dy2*(pjcxc(6,i,l)+
     .                    pjcxc(2,i,l)*(yc(jc)-sy)+
     .              0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
      enddo

      do j=0,ncyc(l)
        sx=falfayc(1,j,l)
        ie=ieyc(j,l)
        ic=icyc(j,l)
        ic1=ic+1
        pcdxx(1,j,l)=-dx2*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic1)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
        pcdxx(2,j,l)=dx2*(pjcyc(5,j,l)+
     .                    pjcyc(1,j,l)*(xc(ic)-sx)+
     .              0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
      enddo

      ENDDO

c  pressure

      !----------------------------------------------------------------!
      !     --- Solving the poisson equation: rhs, bc & fft ---        !
      !----------------------------------------------------------------!

      do j=1,ny
        do i=1,nx
          prhs(i,j)=0.0d0 
        enddo
      enddo

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

      do j=1,ny
        pw(j)=0.0d0
        pe(j)=0.0d0
      enddo

      do i=1,nx
        ps(i)=0.0d0
        pn(i)=0.0d0
      enddo

      call poisson_fft(pw,pe,ps,pn,2,2,2,2,1)  

      do j=1,ny
        do i=1,nx
          funHp(i,j) = p(i,j)
        enddo
      enddo
      funHp = funHp - p(1,1)
      !funHp = funHp - p(nx/2,ny/2)

      !----------------------------------------------------------------!
      !      --- Step functions for the velocities position ---        !
      !----------------------------------------------------------------!
     
      funHv = funHp 
      funHu = funHp 

      DO l=1,ms

        ! irregular points at lines xc (velocity v)
        do m=0,ncxc(l)
           if (jcxc(m,l)==jexc(m,l)+1) then
             funHv(ixc(m,l),jexc(m,l)+1) = funHp(ixc(m,l),jcxc(m,l)+1)  
           endif
        enddo

        ! Irregular points at lines yc (velocity u)
        do m=0,ncyc(l)
           if (icyc(m,l)==ieyc(m,l)+1) then
             funHu(ieyc(m,l)+1,jyc(m,l)) = funHp(icyc(m,l)+1,jyc(m,l))  
           endif
        enddo

      ENDDO

      !----------------------------------------------------------------!
      !                  --- Calculating Errors ---                    !
      !----------------------------------------------------------------!
 
      DisplayError=0
      IF (DisplayError==1) THEN
      maxx=0.0d0
      do i=1,nx
        do j=1,ny
          rr=dsqrt(xc(i)**2+yc(j)**2)
          if (rr>=0.5d0) then
             dd = abs(funHp(i,j)-0.0d0)
          else
             dd = abs(funHp(i,j)-1.0d0)
          endif
          maxx = max(maxx,dd)
        enddo
      enddo
      write(*,*),'ErrorFunHp=',maxx
      maxx=0.0d0
      do i=1,nx
        do j=1,ny
          rr=dsqrt(xe(i)**2+yc(j)**2)
          if (rr>=0.5d0) then
             dd = abs(funHu(i,j)-0.0d0)
          else
             dd = abs(funHu(i,j)-1.0d0)
          endif
          maxx = max(maxx,dd)
        enddo
      enddo
      write(*,*),'ErrorFunHu=',maxx
      maxx=0.0d0
      do i=1,nx
        do j=1,ny
          rr=dsqrt(xc(i)**2+ye(j)**2)
          if (rr>=0.5d0) then
             dd = abs(funHv(i,j)-0.0d0)
          else
             dd = abs(funHv(i,j)-1.0d0)
          endif
          maxx = max(maxx,dd)
        enddo
      enddo
      write(*,*),'ErrorFunHv=',maxx
      ENDIF

c    !----------------------------------------------------------------!
c    !      ---  Writting the step functions in a data file ---       !
c    !----------------------------------------------------------------!

400   format(1x,2000e16.6)
      open(unit=19,file='DAT/StepFun/Heaviside_p.dat',status='unknown')
      open(unit=29,file='DAT/StepFun/Heaviside_v.dat',status='unknown')
      open(unit=39,file='DAT/StepFun/Heaviside_u.dat',status='unknown')

      do j=1,ny
        write(19,400)(funHp(i,j),i=1,nx)
        write(29,400)(funHv(i,j),i=1,nx)
        write(39,400)(funHu(i,j),i=1,nx)
      enddo
      close(19)
      close(29)
      close(39)

c    !-----------------------------------------------------------------!

      return
      end
      
