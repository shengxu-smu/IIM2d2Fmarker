c
      !****************************************************************!
      !          ***   CALCULATE (du/dn)+ and (du/dn)-   ***           !
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
      !                       - dudn_pos (at the Lagragrian points)    !
      !                       - dudn_neg (at the Lagragrian points)    !
      !                       - dvdn_pos                               !
      !                       - dvdn_neg                               !
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
      !****************************************************************!
c
      subroutine dUdn(curr_u,curr_v,Uxsys,Vxsys,
     .                dudn_pos,dudn_neg,dvdn_pos,dvdn_neg)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,many,nany,ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
      parameter(many=3,nany=3)
      real*8 fac,ds,xx,yy,xt,yt,xn,yn,gacobi,usi,vsi,ut,un,et,en,ed
      real*8 foo,uu(3,2),vv(3,2), su, sv, so, sd
      real*8 xa(many),ya(many),xb(many),yb(many)
      real*8 r(3,ns),w(3,ns)
      real*8 curr_u(0:nx,0:ny+1), curr_v(0:nx+1,0:ny)
      real*8 dudtao, dvdtao     
      real*8 dudn_pos(0:ns,ms), dudn_neg(0:ns,ms)
      real*8 dvdn_pos(0:ns,ms), dvdn_neg(0:ns,ms)
      real*8 Jump_mududn, Jump_mudvdn
      integer TheOption
      real*8 Uxsys(0:ns,ms),Vxsys(0:ns,ms)

      ds=1.01d0*sqrt(dx*dx+dy*dy)

      DO l=1,ms
        do m=0,ns-1   
     
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !           --- Tangent and normal vector (outside) ---       !
        !                  k=1 outside, k=2 inside                    !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

          do k=1,2

          xt=taox(m,l)
          yt=taoy(m,l)
        
          if(k.eq.1) then 
            xn=yt
            yn=-xt
          else
            xn=-yt
            yn=xt
          endif
       
c          gacobi=dsqrt(xn*xn+yn*yn)
c          xn=xn/gacobi
c          yn=yn/gacobi
          su=us(m,l)
          sv=vs(m,l)
        
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
      
          enddo !k (InOut) 
         
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !   --- One-sided du/dn+ and du/dn- at each Lagr point ---     !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        
          TheOption=1

          !-------------------------------------------!
          ! Option 1 : without the interface velocity !
          !-------------------------------------------!
         
          if (TheOption==1) then

            dudn_pos(m,l) =  (1/ds)*( -2.5d0*uu(1,1) + 4.0d0*uu(2,1)
     .                                -1.5d0*uu(3,1) )
            dudn_neg(m,l) = -(1/ds)*( -2.5d0*uu(1,2) + 4.0d0*uu(2,2)
     .                                -1.5d0*uu(3,2) )
            dvdn_pos(m,l) =  (1/ds)*( -2.5d0*vv(1,1) + 4.0d0*vv(2,1)
     .                                -1.5d0*vv(3,1) )
            dvdn_neg(m,l) = -(1/ds)*( -2.5d0*vv(1,2) + 4.0d0*vv(2,2)
     .                                -1.5d0*vv(3,2) )
         
          !-------------------------------------------!
          ! Option 2 : Using the intrface velocity    !
          !-------------------------------------------!
         
          elseif (TheOption==2) then

            dudn_pos(m,l) =  (1/ds)*( -1.5d0*Uxsys(m,l) + 2.0d0*uu(1,1)
     .                                -0.5d0*uu(2,1) )
            dudn_neg(m,l) = -(1/ds)*( -1.5d0*Uxsys(m,l) + 2.0d0*uu(1,2)
     .                                -0.5d0*uu(2,2) )
            dvdn_pos(m,l) =  (1/ds)*( -1.5d0*Vxsys(m,l) + 2.0d0*vv(1,1)
     .                                -0.5d0*vv(2,1) )
            dvdn_neg(m,l) = -(1/ds)*( -1.5d0*Vxsys(m,l) + 2.0d0*vv(1,2)
     .                                -0.5d0*vv(2,2) )
          endif
  
        enddo  !m (Lagragian)
           
      ENDDO ! Objects
      !----------------------------------------------------------------!
      !      Writing Numerical interface velocity  at (xs,ys)          !
      !----------------------------------------------------------------!

400   format(1x,2000e16.6)

      open(unit=17,file='DAT/N_dudnp.dat',status='unknown')
      open(unit=27,file='DAT/N_dudnm.dat',status='unknown') 
      open(unit=37,file='DAT/N_dvdnp.dat',status='unknown')
      open(unit=47,file='DAT/N_dvdnm.dat',status='unknown')     

      do m=0,ns-1
        write(17,400)(dudn_pos(m,l),l=1,ms)
        write(27,400)(dudn_neg(m,l),l=1,ms)
        write(37,400)(dvdn_pos(m,l),l=1,ms)
        write(47,400)(dvdn_neg(m,l),l=1,ms)
      enddo
       
      close(17)
      close(27)
      close(37)
      close(47)
            
      !----------------------------------------------------------------!
      
      return
      end    
