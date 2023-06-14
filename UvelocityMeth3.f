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
      subroutine UvelocityMeth3(curr_u,curr_v,Uinter,Vinter)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,many,nany,ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
      parameter(many=3,nany=1)
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
      real*8 MatrixA(1:2*ns,1:2*ns),vectorB(1:2*ns,1)
      real*8 MatrixAT(1:2*ns,1:2*ns)
      real*8 cc1(1:ns),cc2(1:ns),cc3(1:ns), cgamma 
      integer info, ipiv(2*ns)
      real*8  u_sys(1:2*ns,1)
      
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
        !                   --- Big matrix ---                         !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        cgamma = ds/(mu_pos+mu_neg)*(mu_pos-mu_neg)*1.0d0/(2*pi/ns)
 
        ! Constants for the Big Matrix

        do i=1,ns    
          xt=taox(i,l)
          yt=taoy(i,l)
          gacobi=dsqrt(xt**2+yt**2)
          xt = xt/gacobi
          yt = yt/gacobi
          xn= yt
          yn= -xt       
          cc1(ns-(i-1)) = -cgamma*(xn*xt + xn*xt)   
          cc2(ns-(i-1)) =  cgamma*(yn*xt + xn*yt)  
          cc3(ns-(i-1)) = -cgamma*(yn*yt + yn*yt)  
        enddo
  
        ! Matrix
   
        do i=1,ns 
          if (i==1) then
            MatrixA(i,i)   =   1
            MatrixA(i,i+1) = -cc1(i)  
            MatrixA(i,ns)  =  cc1(i) 
        
            MatrixA(i,ns + (i+1)) = -cc2(i) 
            MatrixA(i,ns + (ns))  =  cc2(i) 
   
            vectorB(i,1)  = u_initial(i-1) 
          elseif (i==ns) then
            MatrixA(i,i-1) =  cc1(i) 
            MatrixA(i,i)   =   1 
            MatrixA(i,1)   = -cc1(i) 
        
            MatrixA(i,ns + (i-1)) =  cc2(i) 
            MatrixA(i,ns + 1)     = -cc2(i) 
   
            vectorB(i,1)  = u_initial(i-1) 
          else
            MatrixA(i,i-1) =  cc1(i) 
            MatrixA(i,i)   =   1 
            MatrixA(i,i+1) = -cc1(i) 
        
            MatrixA(i,ns + i-1) =  cc2(i) 
            MatrixA(i,ns + i+1) = -cc2(i) 
        
            vectorB(i,1)  = u_initial(i-1) 
          endif
        enddo

        do i=1,ns    
          if (i==1) then
            MatrixA(ns+i, ns) =  cc2(i) 
            MatrixA(ns+i,i)   =   0 
            MatrixA(ns+i,i+1) = -cc2(i) 
        
            MatrixA(ns+i,ns +  ns) =  cc3(i) 
            MatrixA(ns+i,ns +   i) =   1 
            MatrixA(ns+i,ns + i+1) = -cc3(i) 
   
            vectorB(ns+i,1)  = v_initial(i-1) 
          elseif (i==ns) then
            MatrixA(ns+i,i-1) =  cc2(i) 
            MatrixA(ns+i,i)   =   0 
            MatrixA(ns+i,1)   = -cc2(i) 
       
            MatrixA(ns+i,ns + i-1) =  cc3(i) 
            MatrixA(ns+i,ns +   i) =   1 
            MatrixA(ns+i,ns +   1) = -cc3(i) 
      
            vectorB(ns+i,1)  = v_initial(i-1) 
          else
            MatrixA(ns+i,i-1) =  cc2(i) 
            MatrixA(ns+i,i)   =   0 
            MatrixA(ns+i,i+1) = -cc2(i) 
       
            MatrixA(ns+i,ns + i-1) =  cc3(i) 
            MatrixA(ns+i,ns +   i) =   1 
            MatrixA(ns+i,ns + i+1) = -cc3(i) 
      
            vectorB(ns+i,1)  = v_initial(i-1) 
          endif
        enddo

      ! Solve the linear system Au=B

       do i=1,2*ns
         do j=1,2*ns
           MatrixAT(i,j)=MatrixA(j,i)
         enddo
       enddo
     
       call dgesv(2*ns,1,MatrixAT,2*ns,ipiv,vectorB,2*ns,info)  
       u_sys = vectorB

      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF

      ! Assign the solution

      do  m=0,ns-1
         Uinter(m,l) = u_sys(m+1,1)
         Vinter(m,l) = u_sys(ns+m+1,1)
      enddo 
         Uinter(ns,l) = Uinter(0,l)
         Vinter(ns,l) = Vinter(0,l)   

      ENDDO !l Objects


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                   --- Smoothing ---                          !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      DO l=1,ms
        do m=0,ns-1
          r(1,m+1)=Uinter(m,l)
          r(2,m+1)=Vinter(m,l)
        enddo
        call vrfftf(2,ns,r,w,2,wsave)
        do m=16,ns
          r(1,m)=0.0d0
          r(2,m)=0.0d0
        enddo
        call vrfftb(2,ns,r,w,2,wsave)
        do m=0,ns-1
          Uinter(m,l)=r(1,m+1)
          Vinter(m,l)=r(2,m+1)
          us(m,l)=Uinter(m,l)
          vs(m,l)=Vinter(m,l)
        enddo
        Uinter(ns,l)=Uinter(0,l)
        Vinter(ns,l)=Vinter(0,l)
        us(ns,l)=us(0,l)
        vs(ns,l)=vs(0,l)
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
