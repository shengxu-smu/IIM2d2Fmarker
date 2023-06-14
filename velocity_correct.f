c-----------------------------------------------------------------------
c
      subroutine velocity_correct(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'old.inc'
      integer ie,ic,ic1,je,jc,jc1,krk
      real*8 tmp(1,0:ns),ccc(1,0:ns),foo,fac,grad
      real*8 gacobi,gacobi2,dgacobi,dnx,dny,r1,r2,r3,ts,tm
      real*8 fn,tx,ty,dft,dfn,dtx,dty,ddft,ddfn
      real*8 dpxjc,dpyjc,ddpxjc,ddpyjc
      real*8 alfa0,alfa1,alpha,signx,signy
      real*8 aa(6),bb(6),f(6),sx,sy
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)

c jc_rhs:

      DO l=1,ms

      do m=0,ns
        tmp(1,m)=fy0(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,1)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,dfy0(m,l),1)
      enddo
      dfy0(ns,l)=dfy0(0,l)

      do m=0,ns
        tmp(1,m)=dfy0(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,1)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,ddfy0(m,l),1)
      enddo
      ddfy0(ns,l)=ddfy0(0,l)

      ENDDO

c jc_firstsecond

      DO l=1,ms

      do m=0,ns-1
        tx=taox(m,l)
        ty=taoy(m,l)
        dtx=dtaox(m,l)
        dty=dtaoy(m,l)
        fn=fy0(m,l)
        dft=0.0d0
        dfn=dfy0(m,l)
        ddft=0.0d0
        ddfn=ddfy0(m,l)

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
        ddpxjc=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddpyjc=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2

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
        if(jc.eq.je) then
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc1)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
        else 
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
        endif
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
        if(ic.eq.ie) then
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic1)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
        else
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
        endif
      enddo

      ENDDO

c  pressure

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

c correction:

      do j=1,ny
        do i=1,nx-1
          grad=-dx1*(p(i+1,j)-p(i,j))
          urk(i,j,krk)=urk(i,j,krk)+grad
        enddo
      enddo

      do j=1,ny-1
        do i=1,nx
          grad=-dy1*(p(i,j+1)-p(i,j))
          vrk(i,j,krk)=vrk(i,j,krk)+grad
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do j=0,ncyc(l)
            urk(icyc(j,l),jyc(j,l),krk)=urk(icyc(j,l),jyc(j,l),krk)
     .                                  -ucpdx(j,l)
          enddo
          do i=0,ncxc(l)
            vrk(ixc(i,l),jcxc(i,l),krk)=vrk(ixc(i,l),jcxc(i,l),krk)
     .                                  -vcpdy(i,l)
          enddo
        enddo
      endif

      do j=1,ny
        do i=1,nx
          p(i,j)=o(i,j)+p(i,j)
        enddo
      enddo

      do l=1,ms
        do m=0,ns
          fy(m,l)=fy(m,l)+fy0(m,l)
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
