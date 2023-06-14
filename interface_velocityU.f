c-----------------------------------------------------------------------
c
      subroutine interface_velocityU
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      real*8 tao_x, tao_y, nor_x, nor_y, norm_tao, Delta_n
      real*8 x_star, y_star
      real*8 hxu_pos, hxu_neg, hyu_pos, hyu_neg
      real*8 hxv_pos, hxv_neg, hyv_pos, hyv_neg
      real*8 pa1, pa2, pa3, pa4
      real*8 u_s(1:2), v_s(1:2)
      real*8 u_Outside, v_Outside
      real*8 u_Inside, v_Inside
      real*8 uI_before, uI_after, vI_before, vI_after
      real*8 Jjacobi, dalphaa
      real*8 u_initial(0:ns-1), v_initial(0:ns-1)
      real*8 Du_Dtau(0:ns-1), Dv_Dtau(0:ns-1)
      integer Iu, Ju, Iv, Jv
     
c----------------------------------------------------------------------
c    *  0ver all the objects l=1,ms
c    *  Over all the Lagrangian points m=0,ns-1
c----------------------------------------------------------------------     

      DO l=1,ms

        do m=0,ns-1

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !     ---  Finding the outside and inside points ---         !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! unitary tangent and normal vector(outside) 
        
        norm_tao = dsqrt(taox(m,l)**2+taoy(m,l)**2)
        tao_x = taox(m,l)/norm_tao
        tao_y = taoy(m,l)/norm_tao
        nor_x =  tao_y    
        nor_y = -tao_x 

        ! outside and inside points

        Delta_n = 1.01d0*dsqrt(dx**2+dy**2)               
        xsout(m,l) = xs(m,l) + Delta_n*nor_x
        ysout(m,l) = ys(m,l) + Delta_n*nor_y
        xsin(m,l)  = xs(m,l) - Delta_n*nor_x
        ysin(m,l)  = ys(m,l) - Delta_n*nor_y
     
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !           ---  Interpolation first option ---             !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do k=1,2
           
           ! Considering the outside and inside points

           if (k==1) then
              x_star = xsout(m,l)
              y_star = ysout(m,l)
           else
              x_star = xsin(m,l)
              y_star = ysin(m,l)
           endif       
           
           ! Finding the indices of the grid points for u an v           
           
           i=1
           do while (xe(i)<x_star)
              Iu = i
              i  = i+1 
           end do
           j=1
           do while (yc(j)<y_star)
              Ju = j
              j   = j+1
           end do
           i=1
           do while (xc(i)<x_star)
              Iv = i
               i = i+1
           end do
           j=1
           do while(ye(j)<y_star)
              Jv = j
               j = j+1
           end do
            
           ! Distance from star to grid points 
           
           hxu_neg = x_star  - xe(Iu)
           hxu_pos = xe(Iu+1)- x_star
           hyu_neg = y_star  - yc(Ju)
           hyu_pos = yc(Ju+1)- y_star 
          
           hxv_neg = x_star  - xc(Iv)
           hxv_pos = xc(Iv+1)- x_star
           hyv_neg = y_star  - ye(Jv)
           hyv_pos = ye(Jv+1)- y_star
           
           ! Interpolation Formulas
            
           pa1 = hyu_neg*hxu_neg*u(Iu+1,Ju+1)
           pa2 = hyu_neg*hxu_pos*u(Iu,Ju+1)
           pa3 = hyu_pos*hxu_neg*u(Iu+1,Ju)  
           pa4 = hyu_pos*hxu_pos*u(Iu,Ju)
           u_s(k) = (pa1 + pa2 + pa3 + pa4)/(dx*dy) 
           
           pa1 = hyv_neg*hxv_neg*v(Iv+1,Jv+1)
           pa2 = hyv_neg*hxv_pos*v(Iv,Jv+1)
           pa3 = hyv_pos*hxv_neg*v(Iv+1,Jv)  
           pa4 = hyv_pos*hxv_pos*v(Iv,Jv)
           v_s(k) =(pa1 + pa2 + pa3 + pa4)/(dx*dy)
        
        enddo !(end of k)
       
        !   The outside and inside values of the velocity
        
        u_Outside = u_s(1)
        v_Outside = v_s(1)
        u_Inside  = u_s(2)
        v_Inside  = v_s(2)

        ! The final interpolation ( first option )

        u_initial(m)=(mu_pos*u_Outside+mu_neg*u_Inside)/(mu_pos+mu_neg)
        v_initial(m)=(mu_pos*v_Outside+mu_neg*v_Inside)/(mu_pos+mu_neg)
        
        enddo !(end of m)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !        ---   Interpolation Second option ---              !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do m=0,ns-1

           ! Derivatives: du/dtau and dv/dtau
            
           if (m==0) then
              uI_before = u_initial(ns-1)
              uI_after  = u_initial(1)
              vI_before = v_initial(ns-1)
              vI_after  = v_initial(1)
           elseif (i==ns-1)then
              uI_before = u_initial(ns-2)
              uI_after  = u_initial(0)
              vI_before = v_initial(ns-2)
              vI_after  = v_initial(0) 
           else 
              uI_before = u_initial(m-1)
              uI_after  = u_initial(m+1)
              vI_before = v_initial(m-1)
              vI_after  = v_initial(m+1)
           endif
          
           Jjacobi = dsqrt(taox(m,l**2)+taoy(m,l)**2)
           dalphaa = 2.0*pi/ns 
           Du_Dtau =(1/Jjacobi)*(uI_after-uI_before)/(2*dalphaa)
           Dv_Dtau =(1/Jjacobi)*(vI_after-vI_before)/(2*dalphaa)
            
            
        enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
