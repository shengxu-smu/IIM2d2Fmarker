c-----------------------------------------------------------------------
c
      subroutine jc_p
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer many,nany,ic,jc,id,jd
      parameter(many=3,nany=3)
      real*8 ds,xx,yy,xn,yn,gacobi,foo,pn
      real*8 rhop,rhom,rhoj,ridpdn
      real*8 pp(nany),xa(many),ya(many),xb(many),yb(many)
    
      ds=1.01d0*sqrt(dx*dx+dy*dy)
      rhop=rho_pos
      rhom=rho_neg
      rhoj=rhop-rhom

      DO l=1,ms
 
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                  ---  Outside  ---                    !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      do m=0,ns-1
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
        !  Jump: [1/rho*dp/dn] (CHANGE)  !
        !--------------------------------!

        ridpdn = Prin_jc(m,l,8)

        !--------------------------------!
        ! One-sided jump: [dpdn]         !
        !--------------------------------!

        if(nany.eq.2) then
          if(rhop.gt.rhom) then
            fy0(m,l)=rhom*ridpdn+(rhoj/rhop)*
     .               (-1.0d0*pp(1)+1.0d0*pp(2))/ds
          endif
        endif
        if(nany.eq.3) then
          if(rhop.gt.rhom) then
            fy0(m,l)=rhom*ridpdn+(rhoj/rhop)*
     .        (-2.5d0*pp(1)+4.0d0*pp(2)-1.5d0*pp(3))/ds
          endif
        endif

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                   ---  Inside  ---                    !    
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
        ! One-sided jump: [dpdn]         !
        !--------------------------------!

        if(nany.eq.2) then
          if(rhop.le.rhom) then
            fy0(m,l)=rhop*ridpdn-(rhoj/rhom)*
     .               (-1.0d0*pp(1)+1.0d0*pp(2))/ds
          endif
        endif
        if(nany.eq.3) then
          if(rhop.le.rhom) then
            fy0(m,l)=rhop*ridpdn-(rhoj/rhom)*
     .        (-2.5d0*pp(1)+4.0d0*pp(2)-1.5d0*pp(3))/ds
          endif
        endif

      enddo
      fy0(ns,l)=fy0(0,l)

      ENDDO


      return
      end


c-----------------------------------------------------------------------
