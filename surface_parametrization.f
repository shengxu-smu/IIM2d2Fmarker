c-----------------------------------------------------------------------
c
      subroutine surface_parametrization
      include 'parameter.inc'
      include 'surface.inc'
      integer nd,nu,nn
      real*8 angle,ssl,xsi,ysi,tmp,tangle,tsl,xf,yf,ang1,ang2,ra1,ra2
      real*8 r(2,ns),w(2,ns)

      open(unit=13,file='DAT/shape0.dat',status='unknown')
      open(unit=14,file='DAT/rangle.dat',status='unknown')

c 1: rounded plate; 2: (x/a)^4+(y/b)^2=1; 3: (x/a)^2+(y/b)^2=1

      DO l=1,ms

        xf=a(l)/2
        yf=b(l)/2
        ssl=(2*a(l)+pi*b(l))/ns
        angle=-ssl/yf
        xsi=-xf
        ysi=yf

        do m=0,ns
          IF(jectob(l).eq.1) THEN
          if( xsi+ssl .LE. xf .AND. xsi .GE. -xf .AND.
     +        ysi .GT. 0.0d0 ) then
               xsi=xsi+ssl
          elseif( xsi-ssl .GE. -xf .AND. xsi .LE. xf
     +            .AND. ysi .LT. 0.0d0 ) then
               xsi=xsi-ssl
          elseif( xsi .GE. xf .AND.
     +    abs(datan2(ysi,xsi-xf)+angle) .LE. pi/2) then
               tmp=(xsi-xf)*dcos(angle)-ysi*dsin(angle)
               ysi=(xsi-xf)*dsin(angle)+ysi*dcos(angle)
               xsi=tmp+xf
          elseif( xsi .LE. -xf .AND.
     +    abs(datan2(ysi,xsi+xf)+angle) .GE. pi/2) then
               tmp=(xsi+xf)*dcos(angle)-ysi*dsin(angle)
               ysi=(xsi+xf)*dsin(angle)+ysi*dcos(angle)
               xsi=tmp-xf
          elseif( xsi+ssl .GT. xf .AND. xsi .LT. xf .AND.
     +            ysi .GT. 0.0d0 ) then
               tangle=-(xsi+ssl-xf)/yf
               xsi=xf
               tmp=(xsi-xf)*dcos(tangle)-ysi*dsin(tangle)
               ysi=(xsi-xf)*dsin(tangle)+ysi*dcos(tangle)
               xsi=tmp+xf
          elseif( xsi-ssl .LT. -xf .AND. xsi .GT. -xf .AND.
     +             ysi .LT. 0.0d0 ) then
               tangle=(xsi-ssl+xf)/yf
               xsi=-xf
               tmp=(xsi+xf)*dcos(tangle)-ysi*dsin(tangle)
               ysi=(xsi+xf)*dsin(tangle)+ysi*dcos(tangle)
               xsi=tmp-xf
          elseif( xsi .GE. xf ) then
               tsl=yf*(abs(datan2(ysi,xsi-xf)+angle)-pi/2)
               ysi=-yf
               xsi=xf-tsl
          elseif( xsi .LE. -xf ) then
               tsl=yf*(pi/2-abs(datan2(ysi,xsi+xf)+angle))
               ysi=yf
               xsi=-xf+tsl
          endif
          alfa(ns-m)=(ns-m)*dalfa
          xs0(ns-m,l)=xsi
          ys0(ns-m,l)=ysi         
          ENDIF

          IF(jectob(l).eq.2) THEN
          alfa(m)=dble(m)*dalfa
          angle=alfa(m)
          xf=a(l)*dcos(angle)
          yf=b(l)*dsin(angle)
          if(abs(xf).gt.1.0d-3) then
            tangle=yf/xf
            tmp=(a(l)*a(l)*tangle/b(l))**2.0d0
            xsi=dsqrt((dsqrt(tmp*tmp+4.0d0*a(l)**4.0d0)-tmp)/2.0d0)
            ysi=abs(tangle)*xsi
          else
            xsi=0.0d0
            ysi=b(l)
          endif
          if(angle.le.0.5d0*pi) then
            xs0(m,l)=xsi
            ys0(m,l)=ysi
          endif
          if(angle.gt.0.5d0*pi.and.angle.le.pi) then
            xs0(m,l)=-xsi
            ys0(m,l)=ysi
          endif
          if(angle.gt.pi.and.angle.le.1.5d0*pi) then
            xs0(m,l)=-xsi
            ys0(m,l)=-ysi
          endif
          if(angle.gt.1.5d0*pi.and.angle.le.2.0d0*pi) then
            xs0(m,l)=xsi
            ys0(m,l)=-ysi
          endif
          ENDIF

          IF(jectob(l).eq.3) THEN
          alfa(m)=dble(m)*dalfa
          angle=alfa(m)
          xf=a(l)*dcos(angle)
          yf=b(l)*dsin(angle)
          xs0(m,l)=xf
          ys0(m,l)=yf
          ENDIF
        enddo
        xs0(ns,l)=xs0(0,l)
        ys0(ns,l)=ys0(0,l)

        do m=0,ns-1
          r(1,m+1)=xs0(m,l)
          r(2,m+1)=ys0(m,l)
        enddo

        call vrfftf(2,ns,r,w,2,wsave)
        do m=64,ns
          r(1,m)=0.0d0
          r(2,m)=0.0d0
        enddo
        call vrfftb(2,ns,r,w,2,wsave)

        do m=0,ns-1
          xs0(m,l)=r(1,m+1)
          ys0(m,l)=r(2,m+1)
        enddo
        xs0(ns,l)=xs0(0,l)
        ys0(ns,l)=ys0(0,l)

        do m=0,ns
          xs(m,l)=xsc(l)+
     +             xs0(m,l)*dcos(theta(l))-ys0(m,l)*dsin(theta(l))
          ys(m,l)=ysc(l)+
     +             xs0(m,l)*dsin(theta(l))+ys0(m,l)*dcos(theta(l))
        enddo
        xs(ns,l)=xs(0,l)
        ys(ns,l)=ys(0,l)

        do m=0,ns-1
          xsi=xs0(m,l)
          ysi=ys0(m,l)
          xf=xs0(m+1,l)
          yf=ys0(m+1,l)
          ang1=datan2(ysi,xsi)+pi
          ang2=datan2(yf,xf)+pi
          ra1=dsqrt(xsi*xsi+ysi*ysi)
          ra2=dsqrt(xf*xf+yf*yf)
          nd=int(ang1/dar)+1
          nu=int(ang2/dar)
          if(ang2.gt.ang1) then
            do n=nd,nu
              angle=dble(n)*dar
              ra(n,l)=(ra1*(ang2-angle)+ra2*(angle-ang1))/(ang2-ang1)
            enddo
          else
            nn=4*ns
            do n=nd,nn
              ra(n,l)=ra1
            enddo
            nn=0
            do n=nn,nu
              ra(n,l)=ra2
            enddo
          endif
        enddo

        do m=0,ns
          write(13,100)alfa(m),xs(m,l),ys(m,l)
        enddo

        do n=0,4*ns
          write(14,100)dble(n)*dar,ra(n,l)
        enddo

      ENDDO

      close(13)
      close(14)
100   format(1x,30e16.7)

      return
      end


c-----------------------------------------------------------------------
