c-----------------------------------------------------------------------
c
      subroutine jc_firstsecond
      include 'parameter.inc'
      include 'surface.inc'
      real*8 gacobi,gacobi2,dgacobi,dnx,dny,r1,r2,r3,ts,tm
      real*8 ft,fn,tx,ty,dft,dfn,dtx,dty,ddft,ddfn,sd
      real*8 duxjc,duyjc,dduxjc,dduyjc,dduxyjc,ddpxjc0,ddpxjc1
      real*8 dvxjc,dvyjc,ddvxjc,ddvyjc,ddvxyjc,ddpyjc0,ddpyjc1
      real*8 dpxjc,dpyjc
      real*8 foo,tmp(2,0:ns),ccc(2,0:ns)
      real*8 djdudn(0:ns,ms), djdvdn(0:ns,ms)
      
      real*8 xi,yj,rr
      real*8 aa,bb,dd, Omega1,Omega2, xxc, yyc
      real*8 Q, Apos, Bpos, Aneg, Bneg 
      real*8 dudxp, dudxm, dudyp, dudym 
      real*8 dvdxp, dvdxm, dvdyp, dvdym
      real*8 ddudxp, ddudxm, ddudyp, ddudym
      real*8 ddvdxp, ddvdxm, ddvdyp, ddvdym      

      DO l=1,ms
     
      !------------------------------------------!
      !  d([du/dn])/dalpha  &  d([dv/dn])/dalpha !
      !------------------------------------------!

      do m=0,ns
        tmp(1,m)=Prin_jc(m,l,2)  ! [du/dn]
        tmp(2,m)=Prin_jc(m,l,5)  ! [dv/dn]
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,2)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,djdudn(m,l),1)
        call fdf(alfa(m),alfa(m+1),tmp(2,m),tmp(2,m+1),
     .           ccc(2,m),ccc(2,m+1),alfa(m),foo,djdvdn(m,l),1)
      enddo
      djdudn(ns,l)=djdudn(0,l)
      djdvdn(ns,l)=djdvdn(0,l)

400   format(1x,2000e16.6)

      open(unit=69,file='DAT/dduN.dat',status='unknown')
      open(unit=79,file='DAT/dduE.dat',status='unknown')

      !------------------------------------------!
      !    Formulas of the Cartesian jumps       !
      !------------------------------------------!

      do m=0,ns-1
        tx=taox(m,l)           ! x tan vector (not unitary)
        ty=taoy(m,l)           ! y tan vector (not unitary)
        dtx=dtaox(m,l)         ! d(tx)/d(alfa)
        dty=dtaoy(m,l)         ! d(ty)/d(alfa)
        ft=fx(m,l)             !
        fn=fy(m,l)             ! This is [p] 
        dft=dfx(m,l)           ! J[dp/dn]
        dfn=dfy(m,l)           ! d[p]/d(alfa)
        ddft=ddfx(m,l)
        ddfn=ddfy(m,l)

        gacobi=sqrt(tx*tx+ty*ty)
        gacobi2=tx*tx+ty*ty
        dgacobi=(tx*dtx+ty*dty)/gacobi
        dnx=dty/gacobi-ty*dgacobi/gacobi2
        dny=-dtx/gacobi+tx*dgacobi/gacobi2
        ts=(tx*tx-ty*ty)/gacobi2
        tm=2.0d0*tx*ty/gacobi

        !________________________________________ 
        !      Velocity first derivative         !
        duxjc = ty/gacobi*Prin_jc(m,l,2)
        duyjc =-tx/gacobi*Prin_jc(m,l,2)
        dvxjc = ty/gacobi*Prin_jc(m,l,5)
        dvyjc =-tx/gacobi*Prin_jc(m,l,5)

        !________________________________________
        !     Velocity  second derivative        !

        r1 = -dtx*duxjc-dty*duyjc
        r2 = djdudn(m,l)-dnx*duxjc-dny*duyjc
        r3 = Prin_jc(m,l,3)
        dduxjc  = ( r1*ts+r2*tm+r3*ty*ty)/gacobi2
        dduyjc  = (-r1*ts-r2*tm+r3*tx*tx)/gacobi2
        dduxyjc = ( r1*tm-r2*ts-r3*tx*ty)/gacobi2

        r1= -dtx*dvxjc-dty*dvyjc
        r2= djdvdn(m,l)-dnx*dvxjc-dny*dvyjc
        r3= Prin_jc(m,l,6)
        ddvxjc =(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddvyjc =(-r1*ts-r2*tm+r3*tx*tx)/gacobi2
        ddvxyjc=( r1*tm-r2*ts-r3*tx*ty)/gacobi2

        !________________________________________
        !     Pressure first derivative          !
       
        dpxjc=(tx*dfn+ty*dft)/gacobi2
        dpyjc=(ty*dfn-tx*dft)/gacobi2

        !________________________________________
        !     Pressure second derivative         !
        !         (Aumented approach)            ! 
        !         r3=[Sp]=Prin_jc(m,l,9)         !      

                   
        r1 = ddfn-dtx*dpxjc-dty*dpyjc
        r2=ddft/gacobi-dgacobi*dft/gacobi2-dnx*dpxjc-dny*dpyjc
        r3=Prin_jc(m,l,9)  
        ddpxjc1=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddpyjc1=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2

        !========================================!

        ujc(1,m,l)=duxjc
        vjc(1,m,l)=dvxjc
        pjc(1,m,l)=dpxjc
        ujc(3,m,l)=dduxjc
        vjc(3,m,l)=ddvxjc
        pjc(3,m,l)=ddpxjc1
        ujc(5,m,l)=dduxyjc
        vjc(5,m,l)=ddvxyjc
        pjc(5,m,l)=fn

        ujc(2,m,l)=duyjc
        vjc(2,m,l)=dvyjc
        pjc(2,m,l)=dpyjc
        ujc(4,m,l)=dduyjc
        vjc(4,m,l)=ddvyjc
        pjc(4,m,l)=ddpyjc1
        ujc(6,m,l)=dduxyjc
        vjc(6,m,l)=ddvxyjc
        pjc(6,m,l)=fn
      enddo

      close(69)
      close(79)

      ujc(1,ns,l)=ujc(1,0,l)
      vjc(1,ns,l)=vjc(1,0,l)
      pjc(1,ns,l)=pjc(1,0,l)
      ujc(3,ns,l)=ujc(3,0,l)
      vjc(3,ns,l)=vjc(3,0,l)
      pjc(3,ns,l)=pjc(3,0,l)
      ujc(5,ns,l)=ujc(5,0,l)
      vjc(5,ns,l)=vjc(5,0,l)
      pjc(5,ns,l)=pjc(5,0,l)

      ujc(2,ns,l)=ujc(2,0,l)
      vjc(2,ns,l)=vjc(2,0,l)
      pjc(2,ns,l)=pjc(2,0,l)
      ujc(4,ns,l)=ujc(4,0,l)
      vjc(4,ns,l)=vjc(4,0,l)
      pjc(4,ns,l)=pjc(4,0,l)
      ujc(6,ns,l)=ujc(6,0,l)
      vjc(6,ns,l)=vjc(6,0,l)
      pjc(6,ns,l)=pjc(6,0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
