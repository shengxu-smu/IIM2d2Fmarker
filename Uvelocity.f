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
      subroutine Uvelocity(curr_u,curr_v,Uinter,Vinter)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,many,nany,ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
      parameter(many=2,nany=1)
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
      integer Ix,Iy,itera
      real*8 hx_neg,hx_pos,hy_neg,hy_pos,DDx,DDy,uuu(3,2)  

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
                ya(i)=curr_u(iu+id*(i-1),ju+jd*(j-1))
              enddo
              xb(j)=yc(ju+jd*(j-1))
              call interpolate(xa,ya,many,xx,yb(j),foo)
            enddo
            call interpolate(xb,yb,many,yy,uu(n,k),foo)
         
            do i=1,many
              do j=1,many
                xa(j)=ye(jv+jd*(j-1))
                ya(j)=curr_v(iv+id*(i-1),jv+jd*(j-1))
              enddo
              xb(i)=xc(iv+id*(i-1))
              call interpolate(xa,ya,many,yy,yb(i),foo)
            enddo
            call interpolate(xb,yb,many,xx,vv(n,k),foo)

c            !===== my approximation for u =====!
c            ! I got the same result,(I am not using it)       !
c       
c            DDx = xe(2)-xe(1)
c            DDy = yc(2)-yc(1)
c
c            do i=1,nx
c              if (xe(i)<xx) then
c                Ix = i
c              endif
c            enddo
c       
c            do j=0,ny+1
c              if (yc(j)<yy) then
c                Iy = j
c             endif 
c            enddo
c        
c            hx_neg = xx-xe(Ix)
c            hx_pos = xe(Ix+1)-xx
c            hy_neg = yy-yc(Iy)
c            hy_pos = yc(Iy+1)-yy
c        
c            uuu(n,k)= (hy_neg*hx_neg*curr_u(Ix+1,Iy+1) 
c     .              + hy_neg*hx_pos*curr_u(Ix,Iy+1)
c     .              + hy_pos*hx_neg*curr_u(Ix+1,Iy)  
c     .              + hy_pos*hx_pos*curr_u(Ix,Iy))/(DDx*DDy)
c            
c            !=============== end ==============!
     
          enddo !n (nany)      
      
        ENDDO !k (InOut) 
         
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !       ---  Interpolation Formula (First order) ---          !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

          uOut = uu(1,1)
          vOut = vv(1,1)
          uIn  = uu(1,2)
          vIn  = vv(1,2)
 
          u_initial(m) = (mu_pos*uOut+mu_neg*uIn)/(mu_pos+mu_neg)
          v_initial(m) = (mu_pos*vOut+mu_neg*vIn)/(mu_pos+mu_neg) 
      
        enddo  !m (Lagragian)
          u_initial(ns) = u_initial(0)
          v_initial(ns) = v_initial(0)
       
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !        --- Interpolation Formula (Second order) ---         !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        
        call dUdalfa(u_initial,v_initial,dudalpha,dvdalpha)      
        
        do m=0,ns-1
          
           !------------------------------------------!
           !   Derivatives du/dtau and dv/dtau        !
           !------------------------------------------!
           

c            !========= my approximation for dudtau ===========!
c            ! I got the same result,(I am not using it)       !
c            if (m==0) then
c              uI_before = u_initial(ns-1)
c              uI_after  = u_initial(1)
c              vI_before = v_initial(ns-1)
c              vI_after  = v_initial(1)
c            elseif (i==ns-1)then
c              uI_before = u_initial(ns-2)
c              uI_after  = u_initial(0)
c              vI_before = v_initial(ns-2)
c              vI_after  = v_initial(0) 
c           else 
c              uI_before = u_initial(m-1)
c              uI_after  = u_initial(m+1)
c              vI_before = v_initial(m-1)
c              vI_after  = v_initial(m+1)
c           endif

c           dalphaa = 2.0*pi/ns 
c           dudtao =(1/gacobi)*(uI_after-uI_before)/(2*dalphaa)
c           dvdtao =(1/gacobi)*(vI_after-vI_before)/(2*dalphaa)           
c           !===================================================!           
     
           dudtao = (1/gacobi)*dudalpha(m)
           dvdtao = (1/gacobi)*dvdalpha(m)

           !------------------------------------------!
           ! Jump formulas: [mu*du/dn] and [mu*dv/dn] !
           !------------------------------------------!
         
           xt = taox(m,l)
           yt = taoy(m,l)
           gacobi=dsqrt(xt**2+yt**2)
           xt = 1/gacobi*xt
           yt = 1/gacobi*yt        
           xn =  yt
           yn = -xt 
     
           Jump_mududn = -(mu_pos-mu_neg)*((dudtao*xn+dvdtao*yn)*xt +
     .                                     (dudtao*xt+dvdtao*yt)*xn )
           Jump_mudvdn = -(mu_pos-mu_neg)*((dudtao*xn+dvdtao*yn)*yt +
     .                                     (dudtao*xt+dvdtao*yt)*yn )
           
           !------------------------------------------!
           !           Final Interpolation            !
           !------------------------------------------!
           
           Uinter(m,l)= u_initial(m)-(ds/(mu_pos+mu_neg))*Jump_mududn
           Vinter(m,l)= v_initial(m)-(ds/(mu_pos+mu_neg))*Jump_mudvdn
 
        enddo
        Uinter(ns,l) = Uinter(0,l)
        Vinter(ns,l) = Vinter(0,l)
      
       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
       do k =1,1
         call dUdalfa(Uinter,Vinter,dudalpha,dvdalpha)      
         do m=0,ns-1
            dudtao = (1/gacobi)*dudalpha(m)
            dvdtao = (1/gacobi)*dvdalpha(m)
            xt = taox(m,l)
            yt = taoy(m,l)
            gacobi=dsqrt(xt**2+yt**2)
            xt = 1/gacobi*xt
            yt = 1/gacobi*yt        
            xn =  yt
            yn = -xt 
            Jump_mududn = -(mu_pos-mu_neg)*((dudtao*xn+dvdtao*yn)*xt +
     .                                      (dudtao*xt+dvdtao*yt)*xn )
            Jump_mudvdn = -(mu_pos-mu_neg)*((dudtao*xn+dvdtao*yn)*yt +
     .                                      (dudtao*xt+dvdtao*yt)*yn )
            Uinter(m,l)= u_initial(m)-(ds/(mu_pos+mu_neg))*Jump_mududn
            Vinter(m,l)= v_initial(m)-(ds/(mu_pos+mu_neg))*Jump_mudvdn
          enddo
          Uinter(ns,l) = Uinter(0,l)
          Vinter(ns,l) = Vinter(0,l)
       enddo
       !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
        
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
