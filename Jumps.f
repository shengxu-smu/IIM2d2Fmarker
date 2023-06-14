c-----------------------------------------------------------------------
c
      subroutine Jumps
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      real*8 nada, nada1
      real*8 Exa_Uxsys(0:ns,ms),Exa_Vxsys(0:ns,ms)
      real*8 Num_Uxsys(0:ns,ms),Num_Vxsys(0:ns,ms)
      real*8 Edudnp(0:ns,ms),Edudnm(0:ns,ms)
      real*8 Edvdnp(0:ns,ms),Edvdnm(0:ns,ms)
      real*8 dudnp(0:ns,ms),dudnm(0:ns,ms)
      real*8 dvdnp(0:ns,ms),dvdnm(0:ns,ms)
      real*8 ddudnp(0:ns,ms),ddudnm(0:ns,ms)
      real*8 ddvdnp(0:ns,ms),ddvdnm(0:ns,ms)
      real*8 Eddudnp(0:ns,ms),Eddudnm(0:ns,ms)
      real*8 Eddvdnp(0:ns,ms),Eddvdnm(0:ns,ms)
      real*8 dudnpnorm,dudnmnorm,Uxsysnorm,Vxsysnorm
      integer DisplayError

c      call ExactUvelocity(Exa_Uxsys,Exa_Vxsys)      
c      call ExactNormalDer(Edudnp,Edudnm,Edvdnp,Edvdnm,
c     .                    Eddudnp,Eddudnm,Eddvdnp,Eddvdnm)

      !------------------------------------------!
      ! Calculating the interface velocity       !
      !------------------------------------------!
     
c      call Uvelocity(u,v,Num_Uxsys,Num_Vxsys)
      call UvelocityMeth3(u,v,Num_Uxsys,Num_Vxsys)
c      call Uinterface(u,v,Num_Uxsys,Num_Vxsys)

      !------------------------------------------!
      !    Calling:  Normal derivatives (-/+)    !
      !------------------------------------------!

      call NormalDerivatives(u,v,Num_Uxsys,Num_Vxsys,
     .                       dudnp,dudnm,dvdnp,dvdnm,
     .                       ddudnp,ddudnm,ddvdnp,ddvdnm)
        
      !------------------------------------------!
      !   Calculating the Principal Jumps        !
      !------------------------------------------!

      call Principal_Jumps(Num_Uxsys,Num_Vxsys,
     .                     dudnp,dudnm,dvdnp,dvdnm,
     .                     ddudnp,ddudnm,ddvdnp,ddvdnm)
      
      !------------------------------------------!
      !  Exact: *  Interface velocity            !
      !         *  Normal derivatives            !
      !         *  Principal jumps               !
      !------------------------------------------!

      call ExactUvelocity(Exa_Uxsys,Exa_Vxsys)
      call ExactNormalDer(Edudnp,Edudnm,Edvdnp,Edvdnm,
     .                    Eddudnp,Eddudnm,Eddvdnp,Eddvdnm)
      call ExactJumps

      !------------------------------------------!
      !     ERRORS: dudn+ and U interface        !
      !------------------------------------------!

      DisplayError=0
      if (DisplayError==1) then
      dudnpnorm=0.0d0
      dudnmnorm=0.0d0
      Uxsysnorm=0.0d0
      Vxsysnorm=0.0d0
      DO l=1,ms
        do m=0,ns
           dudnpnorm = max(dudnpnorm,abs(dudnp(m,l)-Edudnp(m,l)))
           dudnmnorm = max(dudnmnorm,abs(dudnm(m,l)-Edudnm(m,l)))
           Uxsysnorm = max(Uxsysnorm,abs(Exa_Uxsys(m,l)-Num_Uxsys(m,l)))
           Vxsysnorm = max(Vxsysnorm,abs(Exa_Vxsys(m,l)-Num_Vxsys(m,l)))
        enddo
      enddo

      write(*,*),'_______________________________________'
      write(*,*),'Error dudn+  =', dudnpnorm
      write(*,*),'Error dudn-  =', dudnmnorm
      write(*,*),'Error Uxsys  =', Uxsysnorm
      write(*,*),'Error Vxsys  =', Vxsysnorm
      write(*,*),'_______________________________________'
      endif

      !------------------------------------------!

      return
      end


c-----------------------------------------------------------------------
