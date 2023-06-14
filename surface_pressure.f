c-----------------------------------------------------------------------
c
      subroutine surface_pressure
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer many,nany,ic,jc,id,jd
      parameter(many=3,nany=3)
      real*8 ds,xx,yy,xn,yn,gacobi,foo,pn,rhop,rhom
      real*8 pp(nany),xa(many),ya(many),xb(many),yb(many)
      real*8 r(1,ns),w(1,ns) 
      real*8 dfx0p,dfx0m
      real*8 a1,a2,a3,a4,xt,yt

      ds=1.01d0*sqrt(dx*dx+dy*dy)
      rhop=rho_pos
      rhom=rho_neg

      DO l=1,ms
 
      do m=0,ns-1

       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
       !                  OUTSIDE: p+                          !
       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

        xn=taoy(m,l)
        yn=-taox(m,l)
        gacobi=dsqrt(xn*xn+yn*yn)
        xn=xn/gacobi
        yn=yn/gacobi
        do n=1,nany
          xx=xs(m,l)+dble(n)*ds*xn
          yy=ys(m,l)+dble(n)*ds*yn
          i=int((xx-x0)/hdx)
          j=int((yy-y0)/hdy)
          if(mod(i,2).eq.0) then
            ic=i/2
          else
            ic=(i+1)/2
          endif
          if(mod(j,2).eq.0) then
            jc=j/2
          else
            jc=(j+1)/2
          endif
          id=int(sign(1.0,xn))
          jd=int(sign(1.0,yn))
          if(id.lt.0) then
            ic=ic+1
          endif
          if(jd.lt.0) then
            jc=jc+1
          endif
          do i=1,many
            do j=1,many
              xa(j)=yc(jc+jd*(j-1))
              ya(j)=p(ic+id*(i-1),jc+jd*(j-1))
            enddo
            xb(i)=xc(ic+id*(i-1))
            call interpolate(xa,ya,many,yy,yb(i),foo)
          enddo
          call interpolate(xb,yb,many,xx,pp(n),foo)
        enddo

        !--------------------------------!
        !    interpolation: p+           !
        !--------------------------------!

        if(nany.eq.2) then
          dfx0p=2.0d0*pp(1)-1.0d0*pp(2)
        endif
        if(nany.eq.3) then
          dfx0p=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
        endif

       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
       !                   INSIDE: p-                          !
       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

        xn=-xn
        yn=-yn
        do n=1,nany
          xx=xs(m,l)+dble(n)*ds*xn
          yy=ys(m,l)+dble(n)*ds*yn
          i=int((xx-x0)/hdx)
          j=int((yy-y0)/hdy)
          if(mod(i,2).eq.0) then
            ic=i/2
          else
            ic=(i+1)/2
          endif
          if(mod(j,2).eq.0) then
            jc=j/2
          else
            jc=(j+1)/2
          endif
          id=int(sign(1.0,xn))
          jd=int(sign(1.0,yn))
          if(id.lt.0) then
            ic=ic+1
          endif
          if(jd.lt.0) then
            jc=jc+ 1
          endif
          do i=1,many
            do j=1,many
              xa(j)=yc(jc+jd*(j-1))
              ya(j)=p(ic+id*(i-1),jc+jd*(j-1))
            enddo
            xb(i)=xc(ic+id*(i-1))
            call interpolate(xa,ya,many,yy,yb(i),foo)
          enddo
          call interpolate(xb,yb,many,xx,pp(n),foo)
        enddo

        !--------------------------------!
        !    interpolation: p-           !
        !--------------------------------!

        if(nany.eq.2) then
          dfx0m = 2.0d0*pp(1)-1.0d0*pp(2)
        endif
        if(nany.eq.3) then
          dfx0m = 3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
        endif

       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
       !                     RHS GMRES                         !
       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
       
       !--------------------------------!
       !   JUMP: d=[p]                  !
       !--------------------------------!

        ! One-sided formula
        xt=taox(m,l)
        yt=taoy(m,l)
        gacobi=dsqrt(xt*xt+yt*yt)
        if (rhop > rhom) then
          dfx0(m,l)=rhop*fy(m,l)+(rhop-rhom)*dfx0m
        else
          dfx0(m,l)=rhom*fy(m,l)+(rhop-rhom)*dfx0p
        endif

        ! Both-sided formula

        !dfx0(m,l)=dfx0p-dfx0m

       !--------------------------------!
       !           Saving d0            !
       !--------------------------------!
        fx0(m,l) = dfx0(m,l)

       !--------------------------------!
       !       rhs_e = []_e - d0        !
       !--------------------------------!

        dfx0(m,l) = Prin_jc(m,l,7)-dfx0(m,l)
      enddo
      dfx0(ns,l)= dfx0(0,l)
      fx0(ns,l) = fx0(0,l)

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      !                     SMOOTHING                         !
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      do m=0,ns-1
        r(1,m+1)=dfx0(m,l)
      enddo
      call vrfftf(1,ns,r,w,1,wsave)
      do m=ns/2,ns
c        r(1,m)=0.0d0
      enddo
      call vrfftb(1,ns,r,w,1,wsave)
      do m=0,ns-1
        dfx0(m,l)=r(1,m+1)
      enddo
      dfx0(ns,l)=dfx0(0,l)

      ENDDO
   
      return
      end


c-----------------------------------------------------------------------
