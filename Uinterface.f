c
      !****************************************************************!
      !       ---  APPROXIMATE VELOCITY AT THE INTERFACE  ---          !
      !                          Sept 25,2010                          !
      !****************************************************************!
      !                                                                !
      !   This program calculates the interface velocity using the     ! 
      !   current field velocity at the grid by interpolation.         !
      !                                                                !
      !    In variables:                                               !
      !                       - curr_u    (current velocity u)         !
      !                       - curr_v    (current velocity v)         !       
      !    Out variables:                                              !
      !                       - Uinter    (at the Lagragrian points)   !
      !                       - Vinter    (at the Lagragrian points)   !
      !                                                                !
      !    Programs:                                                   !
      !                       - interpolate                            !
      !                       - dUdalfa                                !
      !                                                                !
      !    Remark: We write the approximation in the data files:       !
      !                       - N_Uxsys.dat                            !
      !                       - N_Vxsys.dat                            !
      !                                                                !
      !****************************************************************!
c
      subroutine Uinterface(curr_u,curr_v,Uinter,Vinter)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,many,nany,ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
      parameter(many=2,nany=4)
      real*8 fac,ds,xx,yy,xt,yt,xn,yn,gacobi,usi,vsi,ut,un,et,en,ed
      real*8 foo,uu(nany,2),vv(nany,2)
      real*8 xa(many),ya(many),xb(many),yb(many)
      real*8 r(2,ns),w(2,ns)
      real*8 curr_u(0:nx,0:ny+1), curr_v(0:nx+1,0:ny)
      real*8 uOut, vOut, uIn, vIn
      real*8 uI_before, uI_after, vI_before, vI_after, dalphaa
      real*8 u_initial(0:ns), v_initial(0:ns)
      real*8 Uinter(0:ns,ms), Vinter(0:ns,ms) 
      real*8 dudalpha(0:ns), dvdalpha(0:ns), dudtao, dvdtao
      real*8 Jump_mududn, Jump_mudvdn
      real*8 xso, yso, SolUout, SolVout
      integer Ix,Iy
      real*8 dudnp,dudnm,dvdnp,dvdnm,Jdudn,Jdvdn
      real*8 du,dv,alphau,alphav,betau,betav  

      ds=1.01d0*sqrt(dx*dx+dy*dy)

      DO l=1,ms
        do m=0,ns-1   
     
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !          --- Tangent and normal vector (outside) ---        !
        !                 k=1 outside, k=2 inside                     !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        DO k=1,2

          xt=taox(m,l)
          yt=taoy(m,l)
        
          if(k.eq.1) then 
            xn=yt
            yn=-xt
          else
            xn=-yt
            yn=xt
          endif
       
          gacobi=dsqrt(xn*xn+yn*yn)
          xn=xn/gacobi
          yn=yn/gacobi
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !    ---  Finding the velocity at the out/inside points ---    !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
                
          do n=1,nany
            
           !------------------------------------------!
           !            Out/Inside points             !
           !------------------------------------------!

            xx=xs(m,l)+dble(n-1)*ds*xn
            yy=ys(m,l)+dble(n-1)*ds*yn
            i=int((xx-x0)/hdx)
            j=int((yy-y0)/hdy)
            if(mod(i,2).eq.0) then
              ie=i/2
              ic=ie
            else
              ic=(i+1)/2
              ie=ic-1
            endif
            if(mod(j,2).eq.0) then
              je=j/2
              jc=je
            else
              jc=(j+1)/2
              je=jc-1
            endif
            iu=ie
            ju=jc
            iv=ic
            jv=je
            io=ic
            jo=jc

            id=int(sign(1.0,xn))
            jd=int(sign(1.0,yn))
            if(id.lt.0.0d0) then
              iu=iu+1
              iv=iv+1
            endif
            if(jd.lt.0.0d0) then
              ju=ju+1
              jv=jv+1
            endif
          
           !------------------------------------------!
           !        Velocity by interpolation         !
           !------------------------------------------!

            do j=1,many
              do i=1,many
                xa(i)=xe(iu+id*(i-1))
                ya(i)=curr_u(iu+id*(i-1),ju+jd*(j-1))
              enddo
              xb(j)=yc(ju+jd*(j-1))
              call interpolate(xa,ya,many,xx,yb(j),foo)
            enddo
            call interpolate(xb,yb,many,yy,uu(n,k),foo)

            if (n==1)then
              du = dsqrt((xx-xa(1))**2+(yy-xb(2))**2)
              alphau = abs(xx-xa(1))
              betau  = abs(yy-xb(2))
            endif
         
            do i=1,many
              do j=1,many
                xa(j)=ye(jv+jd*(j-1))
                ya(j)=curr_v(iv+id*(i-1),jv+jd*(j-1))
              enddo
              xb(i)=xc(iv+id*(i-1))
              call interpolate(xa,ya,many,yy,yb(i),foo)
            enddo
            call interpolate(xb,yb,many,xx,vv(n,k),foo)

            if (n==1)then
              dv = dsqrt((xx-xb(1))**2+(yy-xa(2))**2)
              alphav = abs(xx-xb(1))
              betav  = abs(yy-xa(2))
            endif

          enddo !n (nany)      
      
        ENDDO !k (InOut) 
         
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !       --- Jump [du/dn] using 3 out/inside points ---         !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

         dudnp = (1/ds)*(-2.5d0*uu(2,1)+4.0d0*uu(3,1)-1.5d0*uu(4,1))
         dudnm =-(1/ds)*(-2.5d0*uu(2,2)+4.0d0*uu(3,2)-1.5d0*uu(4,2))
         dvdnp = (1/ds)*(-2.5d0*vv(2,1)+4.0d0*vv(3,1)-1.5d0*vv(4,1))
         dvdnm =-(1/ds)*(-2.5d0*vv(2,2)+4.0d0*vv(3,2)-1.5d0*vv(4,2))
         Jdudn = dudnp-dudnm
         Jdvdn = dvdnp-dvdnm

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !              ---  Interpolation Formula  ---                 !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

         Uinter(m,l) = uu(1,1)! - betau*(1-alphau)*du*Jdudn 
         Vinter(m,l) = vv(1,1)! - betav*(1-alphav)*dv*Jdvdn
      
        enddo  !m (Lagragian)

        Uinter(ns,l) = Uinter(0,l)
        Vinter(ns,l) = Vinter(0,l)
      
      ENDDO ! Objects

      DO l=1,ms

        do m=0,ns-1
          r(1,m+1)=Uinter(m,l)
          r(2,m+1)=Vinter(m,l)
        enddo
        call vrfftf(2,ns,r,w,2,wsave)
        do m=ns/2,ns
          r(1,m)=0.0d0
        enddo
        call vrfftb(2,ns,r,w,2,wsave)
        do m=0,ns-1
          Uinter(m,l)=r(1,m+1)
          Vinter(m,l)=r(2,m+1)
        enddo
        Uinter(ns,l)=Uinter(0,l)
        Vinter(ns,l)=Vinter(0,l)

      ENDDO

 
      !----------------------------------------------------------------!
      !      Writing Numerical interface velocity  at (xs,ys)          !
      !----------------------------------------------------------------!

400   format(1x,2000e16.6)

      open(unit=17,file='DAT/Uinterface/N_Uxsys.dat',status='unknown')
      open(unit=27,file='DAT/Uinterface/N_Vxsys.dat',status='unknown') 

      do m=0,ns-1
        write(17,400)(Uinter(m,l),l=1,ms)
        write(27,400)(Vinter(m,l),l=1,ms)
      enddo
       
      close(17)
      close(27)
            
      !----------------------------------------------------------------!
      
      return
      end    
