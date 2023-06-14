c
      !****************************************************************!
      !          --- FIRST AND SECOND NORMAL DERIVATIVES ---           !
      !                          Sept 25,2010                          !
      !****************************************************************!
      !                                                                !
      !   This program calculates the derivative of the velocity       ! 
      !   in the normal direction outside and inside using the cur-    !
      !   rent velocity field, the interface velocity and interpo-     !
      !   lation.                                                      !
      !    In variables:                                               !
      !                       - curr_u    (current velocity u)         !
      !                       - curr_v    (current velocity v)         !
      !                       - Uxsys     (interface velocity u)       ! 
      !                       - Vxsys     (interface velocity v)     
      !    Out variables:                                              !
      !                       - dudn_pos  (du/dn)+                     !
      !                       - dudn_neg  (du/dn)-                     !
      !                       - dvdn_pos  (dv/dn)+                     !
      !                       - dvdn_neg  (dv/dn)-                     !
      !                                                                !
      !    Programs:                                                   !
      !                       - interpolate                            !
      !                                                                !
      !    Remark: We write the approximation in the data files:       !
      !                       - N_dudnp.dat                            !
      !                       - N_dudnm.dat                            !
      !                       - N_dvdnp.dat                            !
      !                       - N_dvdnm.dat                            !
      !                                                                !
      !                       - N_ddudnp.dat                           !
      !                       - N_ddudnm.dat                           !
      !                       - N_ddvdnp.dat                           !
      !                       - N_ddvdnm.dat                           !
      !                                                                !
      !****************************************************************!
c
      subroutine NormalDerivatives(cur_u,cur_v,Uxsys,Vxsys,dudnp,dudnm,
     .                          dvdnp,dvdnm,ddudnp,ddudnm,ddvdnp,ddvdnm)

      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,many,nany,ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
      parameter(many=3,nany=5)
      real*8 fac,ds,xx,yy,xt,yt,xn,yn,gacobi,usi,vsi,ut,un,et,en,ed
      real*8 foo,uu(nany,2),vv(nany,2), su, sv, so, sd
      real*8 xa(many),ya(many),xb(many),yb(many)
      real*8 r(nany,ns),w(nany,ns)
      real*8 cur_u(0:nx,0:ny+1), cur_v(0:nx+1,0:ny)
      real*8 dudtao, dvdtao, dudalpha(0:ns), dvdalpha(0:ns)
      real*8 dudtau(0:ns,ms), dvdtau(0:ns,ms)       
      real*8 dudnp(0:ns,ms), dudnm(0:ns,ms)
      real*8 dvdnp(0:ns,ms), dvdnm(0:ns,ms)
      real*8 ddudnp(0:ns,ms),ddudnm(0:ns,ms)
      real*8 ddvdnp(0:ns,ms),ddvdnm(0:ns,ms)
      integer Choose
      real*8 Uxsys(0:ns,ms),Vxsys(0:ns,ms)
      real*8 a0,a1,a2,a3,a4,a5

      ds=1.01d0*sqrt(dx*dx+dy*dy)

      DO l=1,ms
        do m=0,ns-1   
     
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !          --- Tangent and normal vector (outside) ---      !
        !                 k=1 outside, k=2 inside                   !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

          do k=1,2

          xt=taox(m,l)
          yt=taoy(m,l)
        
          if(k==1) then 
            xn=yt
            yn=-xt
          else
            xn=-yt
            yn=xt
          endif
       
          gacobi=dsqrt(xn*xn+yn*yn)
          xn=xn/gacobi
          yn=yn/gacobi
          su=us(m,l)
          sv=vs(m,l)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !   ---  Finding the velocity at the out/inside points ---   !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
                
          do n=1,nany
            
           !------------------------------------------!
           !            Out/Inside points             !
           !------------------------------------------!

            xx=xs(m,l)+dble(n)*ds*xn
            yy=ys(m,l)+dble(n)*ds*yn
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
                ya(i)=cur_u(iu+id*(i-1),ju+jd*(j-1))
              enddo
              xb(j)=yc(ju+jd*(j-1))
              call interpolate(xa,ya,many,xx,yb(j),foo)
            enddo
            call interpolate(xb,yb,many,yy,uu(n,k),foo)
         
            do i=1,many
              do j=1,many
                xa(j)=ye(jv+jd*(j-1))
                ya(j)=cur_v(iv+id*(i-1),jv+jd*(j-1))
              enddo
              xb(i)=xc(iv+id*(i-1))
              call interpolate(xa,ya,many,yy,yb(i),foo)
            enddo
            call interpolate(xb,yb,many,xx,vv(n,k),foo)
     
          enddo !n (nany)      
      
          enddo !k (InOut) 
   
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !       One-sided du/dn+ and du/dn- at each Lagr point         !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        
          Choose=1

         !-----------------------------------------------------!
         ! Option 1 : Using 3 points (without interface vel)   !
         !-----------------------------------------------------!
         
          if (Choose==1) then

         dudnp(m,l)= (1/ds)*(-2.5d0*uu(1,1)+4.0d0*uu(2,1)-1.5d0*uu(3,1))
         dudnm(m,l)=-(1/ds)*(-2.5d0*uu(1,2)+4.0d0*uu(2,2)-1.5d0*uu(3,2))
         dvdnp(m,l)= (1/ds)*(-2.5d0*vv(1,1)+4.0d0*vv(2,1)-1.5d0*vv(3,1))
         dvdnm(m,l)=-(1/ds)*(-2.5d0*vv(1,2)+4.0d0*vv(2,2)-1.5d0*vv(3,2))
         
         ddudnp(m,l)= (uu(1,1)-2.0d0*uu(2,1)+uu(3,1))/(ds**2)
         ddudnm(m,l)= (uu(1,2)-2.0d0*uu(2,2)+uu(3,2))/(ds**2)
         ddvdnp(m,l)= (vv(1,1)-2.0d0*vv(2,1)+vv(3,1))/(ds**2)
         ddvdnm(m,l)= (vv(1,2)-2.0d0*vv(2,2)+vv(3,2))/(ds**2)

         !-----------------------------------------------------!
         ! Option 2 : Using 4 points (without interface vel)   !
         !-----------------------------------------------------!

         elseif (Choose==2) then

         a1 = -26.0d0/6.0d0
         a2 =  57.0d0/6.0d0
         a3 = -42.0d0/6.0d0
         a4 =  11.0d0/6.0d0

         dudnp(m,l)= (1/ds)*(a1*uu(1,1)+a2*uu(2,1) 
     .                      +a3*uu(3,1)+a4*uu(4,1))
         dudnm(m,l)=-(1/ds)*(a1*uu(1,2)+a2*uu(2,2)
     .                      +a3*uu(3,2)+a4*uu(4,2))
         dvdnp(m,l)= (1/ds)*(a1*vv(1,1)+a2*vv(2,1)
     .                      +a3*vv(3,1)+a4*vv(4,1))
         dvdnm(m,l)=-(1/ds)*(a1*vv(1,2)+a2*vv(2,2)
     .                      +a3*vv(3,2)+a4*vv(4,2))

         a1 =  3.0d0
         a2 = -8.0d0
         a3 =  7.0d0
         a4 = -2.0d0

         ddudnp(m,l)= (1/ds**2)*(a1*uu(1,1)+a2*uu(2,1) 
     .                          +a3*uu(3,1)+a4*uu(4,1))
         ddudnm(m,l)= (1/ds**2)*(a1*uu(1,2)+a2*uu(2,2)
     .                          +a3*uu(3,2)+a4*uu(4,2))
         ddvdnp(m,l)= (1/ds**2)*(a1*vv(1,1)+a2*vv(2,1)
     .                          +a3*vv(3,1)+a4*vv(4,1))
         ddvdnm(m,l)= (1/ds**2)*(a1*vv(1,2)+a2*vv(2,2)
     .                          +a3*vv(3,2)+a4*vv(4,2))


         !-----------------------------------------------------!
         ! Option 3 : Using 5 points (without interface vel)   !
         !-----------------------------------------------------!

         elseif (Choose==3) then

         a1 = -77.0d0/12.0d0
         a2 =  107.0d0/6.0d0
         a3 = -39.0d0/2.0d0
         a4 =  61.0d0/6.0d0
         a5 = -25.0d0/12.0d0

         dudnp(m,l)= (1/ds)*(a1*uu(1,1)+a2*uu(2,1) 
     .                      +a3*uu(3,1)+a4*uu(4,1)+a5*uu(5,1))
         dudnm(m,l)=-(1/ds)*(a1*uu(1,2)+a2*uu(2,2)
     .                      +a3*uu(3,2)+a4*uu(4,2)+a5*uu(5,2))
         dvdnp(m,l)= (1/ds)*(a1*vv(1,1)+a2*vv(2,1)
     .                      +a3*vv(3,1)+a4*vv(4,1)+a5*vv(5,1))
         dvdnm(m,l)=-(1/ds)*(a1*vv(1,2)+a2*vv(2,2)
     .                      +a3*vv(3,2)+a4*vv(4,2)+a5*vv(5,2))

         a1 =  71.0d0/12.0d0
         a2 = -59.0d0/3.0d0
         a3 =  49.0d0/2.0d0
         a4 = -41.0d0/3.0d0
         a5 =  35.0d0/12.0d0 

         ddudnp(m,l)= (1/ds**2)*(a1*uu(1,1)+a2*uu(2,1) 
     .                          +a3*uu(3,1)+a4*uu(4,1)+a5*uu(5,1))
         ddudnm(m,l)= (1/ds**2)*(a1*uu(1,2)+a2*uu(2,2)
     .                          +a3*uu(3,2)+a4*uu(4,2)+a5*uu(5,2))
         ddvdnp(m,l)= (1/ds**2)*(a1*vv(1,1)+a2*vv(2,1)
     .                          +a3*vv(3,1)+a4*vv(4,1)+a5*vv(5,1))
         ddvdnm(m,l)= (1/ds**2)*(a1*vv(1,2)+a2*vv(2,2)
     .                          +a3*vv(3,2)+a4*vv(4,2)+a5*vv(5,2))

         !-----------------------------------------------------!
         ! Option 4 : Using the interface velocity (3 points)  !
         !-----------------------------------------------------!
         
         elseif (Choose==4) then

         dudnp(m,l)= (1/ds)*(-1.5d0*Uxsys(m,l)+2.0d0*uu(1,1)
     .                       -0.5d0*uu(2,1))
         dudnm(m,l)=-(1/ds)*(-1.5d0*Uxsys(m,l)+2.0d0*uu(1,2)
     .                       -0.5d0*uu(2,2))
         dvdnp(m,l)= (1/ds)*(-1.5d0*Vxsys(m,l)+2.0d0*vv(1,1)
     .                       -0.5d0*vv(2,1))
         dvdnm(m,l)=-(1/ds)*(-1.5d0*Vxsys(m,l)+2.0d0*vv(1,2)
     .                       -0.5d0*vv(2,2))
         
         ddudnp(m,l)= (1/ds**2)*(Uxsys(m,l)-2.0d0*uu(1,1)+uu(2,1))
         ddudnm(m,l)= (1/ds**2)*(Uxsys(m,l)-2.0d0*uu(1,2)+uu(2,2))
         ddvdnp(m,l)= (1/ds**2)*(Vxsys(m,l)-2.0d0*vv(1,1)+vv(2,1))
         ddvdnm(m,l)= (1/ds**2)*(Vxsys(m,l)-2.0d0*vv(1,2)+vv(2,2))


         !-----------------------------------------------------!
         ! Option 5 : Using the interface velocity (5 points)  !
         !-----------------------------------------------------!
         
         elseif (Choose==5) then

         a0 = -25.0d0/12.0d0
         a1 =  48.0d0/12.0d0
         a2 = -36.0d0/12.0d0
         a3 =  16.0d0/12.0d0
         a4 = -3.0d0/12.0d0

         dudnp(m,l)= (1/ds)*(a0*Uxsys(m,l)+a1*uu(1,1)
     .                      +a2*uu(2,1)+a3*uu(3,1)+a4*uu(4,1))
         dudnm(m,l)=-(1/ds)*(a0*Uxsys(m,l)+a1*uu(1,2)
     .                      +a2*uu(2,2)+a3*uu(3,2)+a4*uu(4,2))
         dvdnp(m,l)= (1/ds)*(a0*Vxsys(m,l)+a1*vv(1,1)
     .                      +a2*vv(2,1)+a3*vv(3,1)+a4*vv(4,1))
         dvdnm(m,l)=-(1/ds)*(a0*Vxsys(m,l)+a1*vv(1,2)
     .                      +a2*vv(2,2)+a3*vv(3,2)+a4*vv(4,2))

         a0 =  35.0d0/12.0d0
         a1 = -104.0d0/12.0d0
         a2 =  114.0d0/12.0d0
         a3 = -56.0d0/12.0d0
         a4 =  11.0d0/12.0d0

         ddudnp(m,l)= (1/ds**2)*(a0*Uxsys(m,l)+a1*uu(1,1)+a2*uu(2,1)
     .                          +a3*uu(3,1)+a4*uu(4,1))
         ddudnm(m,l)= (1/ds**2)*(a0*Uxsys(m,l)+a1*uu(1,2)+a2*uu(2,2)
     .                          +a3*uu(3,2)+a4*uu(4,2))
         ddvdnp(m,l)= (1/ds**2)*(a0*Vxsys(m,l)+a1*vv(1,1)+a2*vv(2,1)
     .                          +a3*vv(3,1)+a4*vv(4,1))
         ddvdnm(m,l)= (1/ds**2)*(a0*Vxsys(m,l)+a1*vv(1,2)+a2*vv(2,2)
     .                          +a3*vv(3,2)+a4*vv(4,2))

         endif      

        enddo  !m (Lagragian)
         
        dudnp(ns,l)  = dudnp(0,l)  
        dudnm(ns,l)  = dudnm(0,l)
        dvdnp(ns,l)  = dvdnp(0,l)
        dvdnm(ns,l)  = dvdnm(0,l)
      
        ddudnp(ns,l)  = ddudnp(0,l)  
        ddudnm(ns,l)  = ddudnm(0,l)
        ddvdnp(ns,l)  = ddvdnp(0,l)
        ddvdnm(ns,l)  = ddvdnm(0,l)

      ENDDO ! Objects


c-----------------------------------------------------------------------
c     Writing the approximate normal derivatives
c-----------------------------------------------------------------------

400   format(1x,2000e16.6)

      open(unit=19,file='DAT/NormalDer/N_dudnp.dat',status='unknown')
      open(unit=29,file='DAT/NormalDer/N_dudnm.dat',status='unknown')
      open(unit=39,file='DAT/NormalDer/N_dvdnp.dat',status='unknown')  
      open(unit=49,file='DAT/NormalDer/N_dvdnm.dat',status='unknown')  
      
      do m=0,ns-1
        write(19,400)(dudnp(m,l),l=1,ms)
        write(29,400)(dudnm(m,l),l=1,ms)
        write(39,400)(dvdnp(m,l),l=1,ms)
        write(49,400)(dvdnm(m,l),l=1,ms)
      enddo
         
      close(19)
      close(29)
      close(39)
      close(49)

c-----------------------------------------------------------------------
c     Writing the approximate second normal derivatives
c-----------------------------------------------------------------------

      open(unit=86,file='DAT/NormalDer/N_ddudnp.dat',status='unknown')
      open(unit=26,file='DAT/NormalDer/N_ddudnm.dat',status='unknown')
      open(unit=36,file='DAT/NormalDer/N_ddvdnp.dat',status='unknown')  
      open(unit=46,file='DAT/NormalDer/N_ddvdnm.dat',status='unknown') 
      
      do m=0,ns-1
        write(86,400)(ddudnp(m,l),l=1,ms)
        write(26,400)(ddudnm(m,l),l=1,ms)
        write(36,400)(ddvdnp(m,l),l=1,ms)
        write(46,400)(ddvdnm(m,l),l=1,ms)
      enddo
         
      close(86)
      close(26)
      close(36)
      close(46)

c-----------------------------------------------------------------------
      
      return
      end    
