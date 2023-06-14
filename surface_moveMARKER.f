c-----------------------------------------------------------------------  
c
      subroutine surface_moveMARKER(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'old.inc'
      include 'field.inc'
      integer krk
      real*8 fac
      real*8 gacobi,xn,yn,xt,yt,xtt,ytt
      real*8 Fspeed
      real*8 Uxsys(0:ns,ms),Vxsys(0:ns,ms)
      real*8 r(2,ns),w(2,ns)

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                  ---  MARKER METHOD ---                       !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 

      call UvelocityMeth3(u,v,Uxsys,Vxsys)

      DO l=1,ms
        do m=0,ns-1
        
          !---------------------------------|
          !         Derivatives             |
          !---------------------------------|
          
          xt  = taox(m,l)
          yt  = taoy(m,l)
          gacobi = dsqrt(xt**2+yt**2)
          xn = yt/gacobi
          yn = -xt/gacobi
          xtt = dtaox(m,l)
          ytt = dtaoy(m,l)
         
          !---------------------------------|
          !            Speed F              |
          !---------------------------------|
       
          !Fspeed = 1.0d0
           Fspeed = Uxsys(m,l)*xn +Vxsys(m,l)*yn
          !---------------------------------|
          !            Update               |
          !---------------------------------|
    
           xsrk(m,l,krk) = Uxsys(m,l)!Fspeed*xn 
           ysrk(m,l,krk) = Vxsys(m,l)!Fspeed*yn

        enddo
        xsrk(ns,l,krk)=xsrk(0,l,krk)
        ysrk(ns,l,krk)=ysrk(0,l,krk)
      ENDDO  
      
      !---------------------------------|
      !        Update the xs, ys        |
      !---------------------------------|
    
      if (krk.le.3) then
      DO l=1,ms
        do m=0,ns
          xs(m,l)=xsn(m,l)+fac*dt*xsrk(m,l,krk)
          ys(m,l)=ysn(m,l)+fac*dt*ysrk(m,l,krk)
        enddo
      ENDDO
      endif

      if (krk==4) then
        do l=1,ms
          do m=0,ns
            xs(m,l)=xsn(m,l)+(dt/6.0d0)*(xsrk(m,l,1)+
     .                   2.0d0*(xsrk(m,l,2)+xsrk(m,l,3))+xsrk(m,l,4))
            ys(m,l)=ysn(m,l)+(dt/6.0d0)*(ysrk(m,l,1)+
     .                   2.0d0*(ysrk(m,l,2)+ysrk(m,l,3))+ysrk(m,l,4))
          enddo
        enddo
      endif

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                   --- Smoothing ---                          !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      DO l=1,ms
        do m=0,ns-1
          r(1,m+1)=xs(m,l)
          r(2,m+1)=ys(m,l)
        enddo
        call vrfftf(2,ns,r,w,2,wsave)
        do m=16,ns
          r(1,m)=0.0d0
          r(2,m)=0.0d0
        enddo
        call vrfftb(2,ns,r,w,2,wsave)
        do m=0,ns-1
          xs(m,l)=r(1,m+1)
          ys(m,l)=r(2,m+1)
        enddo
        xs(ns,l)=xs(0,l)
        ys(ns,l)=ys(0,l)
      ENDDO

      return
      end


c-----------------------------------------------------------------------
