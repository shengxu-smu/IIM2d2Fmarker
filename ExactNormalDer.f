c
c     !****************************************************************!
c     !        ---  EXACT SOLUTION FOR THE NORMAL DERIVATIVES  ---     !
c     !                          Sept 25,2010                          !
c     !****************************************************************!
c     !                                                                !
c     !    Remark: It creates the following data files:                !
c     !                       - E_dudnp                                !  
c     !                       - E_dudnm                                !
c     !                       - E_dvdnp                                !
c     !                       - E_dvdnm                                !
c     !                                                                !
c     !****************************************************************!
c
      subroutine  ExactNormalDer(Edudnp,Edudnm,Edvdnp,Edvdnm,
     .                           Eddudnp,Eddudnm,Eddvdnp,Eddvdnm)
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'

      real*8 xi,yj
      real*8 au,av,a_p
      real*8 aa,bb,dd, Omega1,Omega2, xxc, yyc
      real*8 Q, Apos, Bpos, Aneg, Bneg      
      real*8 rr, Ud
      real*8 dudxp, dudxm, dudyp, dudym 
      real*8 dvdxp, dvdxm, dvdyp, dvdym
      real*8 xn, yn, rhop, rhom, xt, yt
      real*8 ppos, pneg
      real*8 dpdxp, dpdxm, dpdyp, dpdym
      real*8 ddpdxp, ddpdxm, ddpdyp, ddpdym
      real*8 Edudnp(0:ns,ms), Edudnm(0:ns,ms)
      real*8 Edvdnp(0:ns,ms), Edvdnm(0:ns,ms)
      real*8 Eddudnp(0:ns,ms),Eddudnm(0:ns,ms)
      real*8 Eddvdnp(0:ns,ms),Eddvdnm(0:ns,ms)
      real*8 dudtau, dvdtau     
      integer OptionProblem

      OptionProblem = 2 

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                          Parameters                            !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      aa        = 0.1d0   ! Inside radius
      bb        = 3.0d0   ! Outside radius
      dd        = 0.5d0   ! interface radius
      Omega1    = 1.0d0   ! Inside velocity
      Omega2    = 2.0d0   ! Outside velocity
      xxc       = 0.d0    ! center (xxc,yyc)  
      yyc       = 0.d0  

      !------------------------------!
      !  Option 1: Two cylinders     !
      !------------------------------!
    
      If (OptionProblem==1) then
     
      Q    = 1/(mu_neg*(1/bb**2-1/dd**2)+ mu_pos*(1/dd**2-1/aa**2))
      Apos = Q*(mu_neg*(Omega1/bb**2-Omega2/dd**2) 
     .          + mu_pos*Omega2*(1/dd**2-1/aa**2)) 
      Bpos = Q*((Omega2-Omega1)*mu_neg)
      Aneg = Q*(mu_neg*Omega1*(1/bb**2-1/dd**2)
     .          +mu_pos*(Omega1/dd**2-Omega2/aa**2))
      Bneg = Q*((Omega2-Omega1)*mu_pos)
 
      !------------------------------!
      ! Option 2: One cylinder       !
      !------------------------------!

      elseif (OptionProblem==2) then

      Apos = 1.0d0!Omega2
      Aneg = 1.0d0!Omega2
      Bpos = 0.0d0
      Bneg = 0.0d0

      endif

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !        Calculating the derivatives of the velocity             !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      Ud = Aneg*dd + Bneg/dd

      DO l=1,ms
        DO m=0,ns

          !-----------------------------!
          ! Unitary tangential vector   !
          ! Unitary normal vector (+/-) !
          !-----------------------------!

          xt = -dsin(alfa(m))
          yt =  dcos(alfa(m))
          xn =  dcos(alfa(m))
          yn =  dsin(alfa(m))
              
          xi = xs(m,l) - xsc(l)
          yj = ys(m,l) - ysc(l)
          rr = dsqrt(xi**2 + yj**2)
         
          !------------------------------!
          !  Cartesian derivatives       !
          !------------------------------!

          dudxm =  2*Bneg*xi*yj/(rr**4)
          dudym = -(Aneg + Bneg/(rr**2)) + 2*Bneg*yj**2/(rr**4) 
          dvdxm =  (Aneg + Bneg/(rr**2)) - 2*Bneg*xi**2/(rr**4)
          dvdym = -2*Bneg*xi*yj/(rr**4)
          dudxp =  2*Bpos*xi*yj/(rr**4)
          dudyp = -(Apos + Bpos/(rr**2)) + 2*Bpos*yj**2/(rr**4) 
          dvdxp =  (Apos + Bpos/(rr**2)) - 2*Bpos*xi**2/(rr**4)
          dvdyp = -2*Bpos*xi*yj/(rr**4)
          dpdxm = Aneg**2*xi + Bneg**2*xi/(rr**4)+2*Aneg*Bneg*xi/(rr**2)
          dpdym = Aneg**2*yj + Bneg**2*yj/(rr**4)+2*Aneg*Bneg*yj/(rr**2)
          dpdxp = Apos**2*xi + Bpos**2*xi/(rr**4)+2*Apos*Bpos*xi/(rr**2)
          dpdyp = Apos**2*yj + Bpos**2*yj/(rr**4)+2*Apos*Bpos*yj/(rr**2)
      
          !------------------------------!
          !      du/dtau &  dv/dtau      !
          !------------------------------!

          dudtau = Ud/dd*(-dcos(alfa(m)))
          dvdtau = Ud/dd*(-dsin(alfa(m)))

          !------------------------------!
          !    du/dn & dv/dn  (+/-)      !
          !------------------------------!

          Edudnp(m,l) = -dsin(alfa(m))*(Apos-Bpos/(dd**2))
          Edudnm(m,l) = -dsin(alfa(m))*(Aneg-Bneg/(dd**2))
          Edvdnp(m,l) =  dcos(alfa(m))*(Apos-Bpos/(dd**2))
          Edvdnm(m,l) =  dcos(alfa(m))*(Aneg-Bneg/(dd**2)) 
           
          !------------------------------!
          ! Second Cartesian derivatives !
          !------------------------------!

          ddpdxm=rho_neg*(Aneg**2-4*(Bneg*xi)**2/(rr**6)+Bneg**2/(rr**4)
     .                + 2*Aneg*Bneg/(rr**2) - 4*Aneg*Bneg*xi**2/(rr**4))
          ddpdym=rho_neg*(Aneg**2-4*(Bneg*yj)**2/(rr**6)+Bneg**2/(rr**4)
     .                + 2*Aneg*Bneg/(rr**2) - 4*Aneg*Bneg*yj**2/(rr**4))
          ddpdxp=rho_pos*(Apos**2-4*(Bpos*xi)**2/(rr**6)+Bpos**2/(rr**4)
     .                + 2*Apos*Bpos/(rr**2) - 4*Apos*Bpos*xi**2/(rr**4))
          ddpdyp=rho_pos*(Apos**2-4*(Bpos*yj)**2/(rr**6)+Bpos**2/(rr**4)
     .                + 2*Apos*Bpos/(rr**2) - 4*Apos*Bpos*yj**2/(rr**4))

          !--------------------------------!
          ! d^2(u)/dn^2 & d^2(v)/dn^2 (+/-)!
          !--------------------------------!
          ! this is not the solution       !

          Eddudnp(m,l) = -dsin(alfa(m))*(2.0d0*Bpos/(dd**3))
          Eddudnm(m,l) = -dsin(alfa(m))*(2.0d0*Bneg/(dd**3))
          Eddvdnp(m,l) =  dcos(alfa(m))*(2.0d0*Bpos/(dd**3))
          Eddvdnm(m,l) =  dcos(alfa(m))*(2.0d0*Bneg/(dd**3))
      
        enddo
        Edudnp(ns,l) = Edudnp(0,l) 
        Edudnm(ns,l) = Edudnm(0,l)
        Edvdnp(ns,l) = Edvdnp(0,l)
        Edvdnm(ns,l) = Edvdnm(0,l)                  

      ENDDO   

c-----------------------------------------------------------------------
c     Writing the exact Normal Derivatives
c-----------------------------------------------------------------------

400   format(1x,2000e16.6)
      open(unit=18,file='DAT/NormalDer/E_dudnp.dat',status='unknown')
      open(unit=28,file='DAT/NormalDer/E_dudnm.dat',status='unknown')
      open(unit=38,file='DAT/NormalDer/E_dvdnp.dat',status='unknown')   
      open(unit=48,file='DAT/NormalDer/E_dvdnm.dat',status='unknown')  
      
      do m=0,ns-1
        write(18,400)(Edudnp(m,l),l=1,ms)
        write(28,400)(Edudnm(m,l),l=1,ms)
        write(38,400)(Edvdnp(m,l),l=1,ms)
        write(48,400)(Edvdnm(m,l),l=1,ms)
      enddo
         
      close(18)
      close(28)
      close(38)
      close(48)

c-----------------------------------------------------------------------
c     Writing the exact Normal Derivatives
c-----------------------------------------------------------------------

      open(unit=17,file='DAT/NormalDer/E_ddudnp.dat',status='unknown')
      open(unit=27,file='DAT/NormalDer/E_ddudnm.dat',status='unknown')
      open(unit=37,file='DAT/NormalDer/E_ddvdnp.dat',status='unknown') 
      open(unit=47,file='DAT/NormalDer/E_ddvdnm.dat',status='unknown') 
      
      do m=0,ns-1
        write(17,400)(Eddudnp(m,l),l=1,ms)
        write(27,400)(Eddudnm(m,l),l=1,ms)
        write(37,400)(Eddvdnp(m,l),l=1,ms)
        write(47,400)(Eddvdnm(m,l),l=1,ms)
      enddo
         
      close(17)
      close(27)
      close(37)
      close(47)

c-----------------------------------------------------------------------

      return
      end

