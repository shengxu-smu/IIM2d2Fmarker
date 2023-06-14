c-----------------------------------------------------------------------
c
      subroutine pressure_refresh
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ic,ic1,jc,jc1
      real*8 tmp(1,0:ns),ccc(1,0:ns),foo
      real*8 gacobi,gacobi2,dgacobi,dnx,dny,r1,r2,r3,ts,tm
      real*8 fn,tx,ty,dft,dfn,dtx,dty,ddft,ddfn
      real*8 dpxjc,dpyjc,ddpxjc0,ddpyjc0,ddpxjc1,ddpyjc1
      real*8 alfa0,alfa1,alpha,signx,signy,sx,sy,hp,hm
      real*8 aa(8),bb(8),f(8)
      real*8 duyp,duym,dvyp,dvym,duxp,duxm,dvxp,dvxm

c jc_rhs:

      DO l=1,ms

      do m=0,ns
        tmp(1,m)=fy(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,1)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,dfy(m,l),1)
      enddo
      dfy(ns,l)=dfy(0,l)

      do m=0,ns
        tmp(1,m)=dfy(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,1)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,ddfy(m,l),1)
      enddo
      ddfy(ns,l)=ddfy(0,l)

      ENDDO

c jc_firstsecond

      DO l=1,ms

      do m=0,ns-1
        tx=taox(m,l)
        ty=taoy(m,l)
        dtx=dtaox(m,l)
        dty=dtaoy(m,l)
        fn=fy(m,l)
        dft=dfx(m,l)
        dfn=dfy(m,l)
        ddft=ddfx(m,l)
        ddfn=ddfy(m,l)

        gacobi=sqrt(tx*tx+ty*ty)
        gacobi2=tx*tx+ty*ty
        dgacobi=(tx*dtx+ty*dty)/gacobi
        dnx=dty/gacobi-ty*dgacobi/gacobi2
        dny=-dtx/gacobi+tx*dgacobi/gacobi2
        ts=(tx*tx-ty*ty)/gacobi2
        tm=2.0d0*tx*ty/gacobi

        dpxjc=(tx*dfn+ty*dft)/gacobi2
        dpyjc=(ty*dfn-tx*dft)/gacobi2

        r1=ddfn-dtx*dpxjc-dty*dpyjc
        r2=ddft/gacobi-dgacobi*dft/gacobi2-
     .     dnx*dpxjc-dny*dpyjc
        r3=0.0d0
        ddpxjc0=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddpyjc0=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2
        r3=1.0d0
        ddpxjc1=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddpyjc1=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2


        pjc(1,m,l)=dpxjc
        pjc(3,m,l)=ddpxjc1
        pjc(5,m,l)=fn
        pjc(7,m,l)=ddpxjc0

        pjc(2,m,l)=dpyjc
        pjc(4,m,l)=ddpyjc1
        pjc(6,m,l)=fn
        pjc(8,m,l)=ddpyjc0
      enddo

      pjc(1,ns,l)=pjc(1,0,l)
      pjc(3,ns,l)=pjc(3,0,l)
      pjc(5,ns,l)=pjc(5,0,l)
      pjc(7,ns,l)=pjc(7,0,l)

      pjc(2,ns,l)=pjc(2,0,l)
      pjc(4,ns,l)=pjc(4,0,l)
      pjc(6,ns,l)=pjc(6,0,l)
      pjc(8,ns,l)=pjc(8,0,l)

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
        do n=1,8
          aa(n)=(pjc(n,m+1,l)-pjc(n,m,l))/dalfa
          bb(n)=(pjc(n,m,l)*alfa1-pjc(n,m+1,l)*alfa0)/dalfa
        enddo

        do n=1,8
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,8
          pjcxc(n,i,l)=f(n)
        enddo

c        sy=falfaxc(1,i,l)
c        ic=ixc(i,l)
c        jc=jcxc(i,l)
c        jc1=jc+1
c        hm=yc(jc)-sy
c        hp=yc(jc1)-sy
c        if(hm*signy.gt.0.0d0) then
c          r3=prhs(ic,jc)
c        else
c          r3=prhs(ic,jc1)
c        endif

        duxp=falfaxc(4,i,l)
        duxm=falfaxc(5,i,l)
        dvxp=falfaxc(6,i,l)
        dvxm=falfaxc(7,i,l)
        duyp=falfaxc(8,i,l)
        duym=falfaxc(9,i,l)
        dvyp=falfaxc(10,i,l)
        dvym=falfaxc(11,i,l)
        r3=2.0d0*signy*((duxp*dvyp-duxm*dvym)-(dvxp*duyp-dvxm*duym))

        pjcxc(3,i,l)=pjcxc(7,i,l)+(pjcxc(3,i,l)-pjcxc(7,i,l))*r3
        pjcxc(4,i,l)=pjcxc(8,i,l)+(pjcxc(4,i,l)-pjcxc(8,i,l))*r3

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
        do n=1,8
          aa(n)=(pjc(n,m+1,l)-pjc(n,m,l))/dalfa
          bb(n)=(pjc(n,m,l)*alfa1-pjc(n,m+1,l)*alfa0)/dalfa
        enddo

        do n=1,8
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,8
          pjcyc(n,j,l)=f(n)
        enddo

c        sx=falfayc(1,j,l)
c        jc=jyc(j,l)
c        ic=icyc(j,l)
c        ic1=ic+1
c        hm=xc(ic)-sx
c        hp=xc(ic1)-sx
c        if(hm*signx.gt.0.0d0) then
c          r3=prhs(ic,jc)
c        else
c          r3=prhs(ic1,jc)
c        endif

        duxp=falfayc(4,j,l)
        duxm=falfayc(5,j,l)
        dvxp=falfayc(6,j,l)
        dvxm=falfayc(7,j,l)
        duyp=falfayc(8,j,l)
        duym=falfayc(9,j,l)
        dvyp=falfayc(10,j,l)
        dvym=falfayc(11,j,l)
        r3=2.0d0*signx*((duxp*dvyp-duxm*dvym)-(dvxp*duyp-dvxm*duym))

        pjcyc(3,j,l)=pjcyc(7,j,l)+(pjcyc(3,j,l)-pjcyc(7,j,l))*r3
        pjcyc(4,j,l)=pjcyc(8,j,l)+(pjcyc(4,j,l)-pjcyc(8,j,l))*r3

        pjcyc(1,j,l)=signx*pjcyc(1,j,l)
        pjcyc(3,j,l)=signx*pjcyc(3,j,l)
        pjcyc(5,j,l)=signx*pjcyc(5,j,l)

        pjcyc(2,j,l)=signy*pjcyc(2,j,l)
        pjcyc(4,j,l)=signy*pjcyc(4,j,l)
        pjcyc(6,j,l)=signy*pjcyc(6,j,l)
      enddo

      ENDDO

c correction_pressure: 
      call correction_pressure

      return
      end


c-----------------------------------------------------------------------
