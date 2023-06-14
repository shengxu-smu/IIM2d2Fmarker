c
      !****************************************************************!
      !                       --- EXACT JUMPS ---                      !
      !                            Sept,2010                           !
      !****************************************************************!
      !                                                                !
      !   Remarks:                                                     ! 
      !   (1) We can choose between the two problems:                  !
      !          * Two cylinders (TheOption=1)                         !
      !          * One cylinder  (TheOption=2)                         !
      !                                                                !
      !    (2) We write the solution in the data files:                !
      !                       - Ejc_u.dat                              !
      !                       - Ejc_dudn.dat                           !
      !                       - Ejc_Lu.dat                             !
      !                       - Ejc_v.dat                              !
      !                       - Ejc_dvdn.dat                           !
      !                       - Ejc_Lv.dat                             !
      !                       - Ejc_p.dat                              !
      !                       - Ejc_1rhodpdn.dat                       !
      !                       - Ejc_Lp.dat                             !
      !                                                                !
      !****************************************************************!
c
      subroutine  ExactJumps
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'

      real*8 xi,yj
      real*8 au,av,a_p
      real*8 aa,bb,dd, Omega1,Omega2, xxc, yyc
      real*8 Q, Apos, Bpos, Aneg, Bneg      
      real*8 rr
      real*8 dudnp, dudnm, dvdnp, dvdnm 
      real*8 dudxp, dudxm, dudyp, dudym 
      real*8 dvdxp, dvdxm, dvdyp, dvdym
      real*8 xn, yn, rhop, rhom, xt, yt
      real*8 ppos, pneg
      real*8 dpdxp, dpdxm, dpdyp, dpdym
      real*8 ddpdxp, ddpdxm, ddpdyp, ddpdym
      real*8 ddudxp, ddudxm, ddudyp, ddudym
      real*8 ddvdxp, ddvdxm, ddvdyp, ddvdym
      real*8 Ejc_u(0:ns,ms), Ejc_dudn(0:ns,ms), Ejc_Lu(0:ns,ms)
      real*8 Ejc_v(0:ns,ms), Ejc_dvdn(0:ns,ms), Ejc_Lv(0:ns,ms)
      real*8 Ejc_p(0:ns,ms), Ejc_1rhodpdn(0:ns,ms), Ejc_Lp(0:ns,ms)
      real*8 dudtau, dvdtau     
      real*8 Ejc_mududn, Ejc_mudvdn, Ud
      real*8 ddudnp, ddudnm, ddvdnp, ddvdnm
      real*8 Ejc_nududn, Ejc_nudvdn, fun1
      integer ProblemOption

      ProblemOption = 2 

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
    
      If (ProblemOption==1) then
     
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

      elseif (ProblemOption==2) then

      Apos = 1.0d0!Omega2
      Aneg = 1.0d0!Omega2
      Bpos = 0.0d0
      Bneg = 0.0d0

      endif

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !        Calculating the derivatives of the velocity             !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      Ud = Aneg*dd + Bneg/dd

      DO l=1,ms
        DO m=0,ns

          !-----------------------------!
          ! Unitary tangential vector   !
          ! Unitary normal vector (+/-) !
          !-----------------------------!

          xt  = -dsin(alfa(m))
          yt  = dcos(alfa(m))
          xn = dcos(alfa(m))
          yn = dsin(alfa(m))
              
          xi = xs(m,l) - xsc(l)
          yj = ys(m,l) - ysc(l)
          rr = dsqrt(xi**2 + yj**2)
  
          !------------------------------!
          !         pressure             !
          !------------------------------!

          pneg =rho_neg*(0.5d0*(Aneg*rr)**2 - 0.5d0*(Bneg/rr)**2
     .           + Aneg*Bneg*dlog(rr**2))
          ppos =rho_pos*(0.5d0*(Apos*rr)**2 - 0.5d0*(Bpos/rr)**2
     .           + Apos*Bpos*dlog(rr**2))
          
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
          !dudtau = dudxp*xt + dudyp*yt
          !dvdtau = dvdxp*xt + dvdyp*yt  

          dudtau = Ud/dd*(-dcos(alfa(m)))
          dvdtau = Ud/dd*(-dsin(alfa(m)))

          !------------------------------!
          !    du/dn & dv/dn  (+/-)      !
          !------------------------------!
          !dudnp = dudxp*xn + dudyp*yn
          !dudnm = dudxm*xn + dudym*yn
          !dvdnp = dvdxp*xn + dvdyp*yn
          !dvdnm = dvdxm*xn + dvdym*yn

          dudnp = -dsin(alfa(m))*(Apos-Bpos/(dd**2))
          dudnm = -dsin(alfa(m))*(Aneg-Bneg/(dd**2))
          dvdnp =  dcos(alfa(m))*(Apos-Bpos/(dd**2))
          dvdnm =  dcos(alfa(m))*(Aneg-Bneg/(dd**2)) 

   
          !------------------------------!
          !   [mu*du/dn] & [mu*dvdn]     !
          !------------------------------!
          !Ejc_mududn = mu_pos*dudnp - mu_neg*dudnm
          !Ejc_mududn = mu_pos*dvdnp - mu_neg*dvdnm

          Ejc_mududn = -dsin(alfa(m))*(mu_pos-mu_neg)*Ud/dd 
          Ejc_mudvdn =  dcos(alfa(m))*(mu_pos-mu_neg)*Ud/dd    
           
          !------------------------------!
          !         [nu*dudn]            !
          !------------------------------!
          
          Ejc_nududn = (mu_pos/rho_pos)*dudnp-(mu_neg/rho_neg)*dudnm
          Ejc_nudvdn = (mu_pos/rho_pos)*dvdnp-(mu_neg/rho_neg)*dvdnm    
 
          !------------------------------!
          !            fun1              !
          !------------------------------!

          fun1 = -xt*Ejc_nududn-yt*Ejc_nudvdn 
     .          +(mu_pos/rho_pos-mu_neg/rho_neg)*(dudtau*xn+dvdtau*yn)
         
          !------------------------------!
          !    d^2/dn^2 (+/-)            !
          !------------------------------!

          ddudnp = -dsin(alfa(m))*(2.0d0*Bpos/(dd**3))
          ddudnm = -dsin(alfa(m))*(2.0d0*Bneg/(dd**3))
          ddvdnp =  dcos(alfa(m))*(2.0d0*Bpos/(dd**3))
          ddvdnm =  dcos(alfa(m))*(2.0d0*Bneg/(dd**3))

          !------------------------------!
          ! Second Cartesian derivatives !
          !------------------------------!

          ! u velocity

          ddudxm=-8.0d0*Bneg*(xi**2)*yj/(rr**6)+2.0d0*Bneg*yj/(rr**4)
          ddudym=-8.0d0*Bneg*(yj**3)/(rr**6)+6.0d0*Bneg*yj/(rr**4)
          ddudxp=-8.0d0*Bpos*(xi**2)*yj/(rr**6)+2.0d0*Bpos*yj/(rr**4)
          ddudyp=-8.0d0*Bpos*(yj**3)/(rr**6)+6.0d0*Bpos*yj/(rr**4) 

          ! v velocity

          ddvdxm=8.0d0*Bneg*(xi**3)/(rr**6)-6.0d0*Bneg*xi/(rr**4)
          ddvdym=8.0d0*Bneg*(yj**2)*xi/(rr**6)-2.0d0*Bneg*xi/(rr**4)
          ddvdxp=8.0d0*Bpos*(xi**3)/(rr**6)-6.0d0*Bpos*xi/(rr**4)
          ddvdyp=8.0d0*Bpos*(yj**2)*xi/(rr**6)-2.0d0*Bpos*xi/(rr**4)
  
          ! pressure
  
          ddpdxm=rho_neg*(Aneg**2-4*(Bneg*xi)**2/(rr**6)+Bneg**2/(rr**4)
     .                + 2*Aneg*Bneg/(rr**2) - 4*Aneg*Bneg*xi**2/(rr**4))
          ddpdym=rho_neg*(Aneg**2-4*(Bneg*yj)**2/(rr**6)+Bneg**2/(rr**4)
     .                + 2*Aneg*Bneg/(rr**2) - 4*Aneg*Bneg*yj**2/(rr**4))
          ddpdxp=rho_pos*(Apos**2-4*(Bpos*xi)**2/(rr**6)+Bpos**2/(rr**4)
     .                + 2*Apos*Bpos/(rr**2) - 4*Apos*Bpos*xi**2/(rr**4))
          ddpdyp=rho_pos*(Apos**2-4*(Bpos*yj)**2/(rr**6)+Bpos**2/(rr**4)
     .                + 2*Apos*Bpos/(rr**2) - 4*Apos*Bpos*yj**2/(rr**4))

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                  Calculating the jumps                         !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
         
          !------------------------------!
          !  Option 1: Two cylinders     !
          !------------------------------!
    
          If (ProblemOption==1) then
          Ejc_u(m,l)       = 0.0d0 
          Ejc_dudn(m,l)    = dudnp  - dudnm
          Ejc_Lu(m,l)      = (ddudxp+ddudyp)-(ddudxm+ddudym)
          !Ejc_Lu(m,l)      = (ddudnp-ddudnm)
          Ejc_v(m,l)       = 0.0d0
          Ejc_dvdn(m,l)    = dvdnp  - dvdnm
          Ejc_Lv(m,l)      = (ddvdxp+ddvdyp)-(ddvdxm+ddvdym)
          !Ejc_Lv(m,l)      = (ddvdnp-ddvdnm)
          Ejc_p(m,l)       = ppos-pneg
          Ejc_1rhodpdn(m,l)= 0.0d0 
          Ejc_Lp(m,l)      = (ddpdxp+ddpdyp)-(ddpdxm+ddpdym)
 
          !------------------------------!
          ! Option 2: One cylinder       !
          !------------------------------!

          elseif (ProblemOption==2) then
          Ejc_u(m,l)        =  0.0d0
          Ejc_dudn(m,l)     =  0.0d0
          Ejc_Lu(m,l)       =  0.0d0
          Ejc_v(m,l)        =  0.0d0
          Ejc_dvdn(m,l)     =  0.0d0
          Ejc_Lv(m,l)       =  0.0d0
          Ejc_p(m,l)        =  0.0d0
          Ejc_1rhodpdn(m,l) =  0.0d0
          Ejc_Lp(m,l)       =  2*(rho_pos-rho_neg)
          
          !------------------------------!
          ! Option 3: u=sin(x)sin(y)     !
          !------------------------------!

          elseif (ProblemOption==3) then
          Ejc_u(m,l)        =  2*(xi*dsin(xi)*dcos(yj)
     .                           -yj*dcos(xi)*dsin(yj))
          Ejc_dudn(m,l)     =  dcos(xi)*dsin(yj)*xn+
     .                         dsin(xi)*dcos(yj)*yn
          Ejc_Lu(m,l)       =  -2*dsin(xi)*dsin(yj)
          Ejc_v(m,l)        =  dsin(xi)*dsin(yj)
          Ejc_dvdn(m,l)     =  dcos(xi)*dsin(yj)*xn+
     .                         dsin(xi)*dcos(yj)*yn
          Ejc_Lv(m,l)       =  -2*dsin(xi)*dsin(yj)
          Ejc_p(m,l)        = -(mu_pos-mu_neg)*2*
     .             ((xi*dsin(xi)*dcos(yj)-yj*dcos(xi)*dsin(yj))*(xt+yt))
          Ejc_1rhodpdn(m,l) = -(xt+yt)*(mu_pos/rho_pos)*
     .             (dcos(xi)*dsin(yj)*xn+dsin(xi)*dcos(yj)*yn)+
     .             ((mu_pos/rho_pos)-(mu_neg/rho_neg))*2*
     .             ((xi*dsin(xi)*dcos(yj)-yj*dcos(xi)*dsin(yj))*(xn+yn))
          Ejc_Lp(m,l)       =  0.0d0
 
          endif        

        enddo
        Ejc_u(ns,l)        = Ejc_u(0,l)
        Ejc_dudn(ns,l)     = Ejc_dudn(0,l)
        Ejc_Lu(ns,l)       = Ejc_Lu(0,l)
        Ejc_v(ns,l)        = Ejc_v(0,l) 
        Ejc_dvdn(ns,l)     = Ejc_dvdn(0,l)
        Ejc_Lv(ns,l)       = Ejc_Lv(0,l)
        Ejc_p(ns,l)        = Ejc_p(0,l)
        Ejc_1rhodpdn(ns,l) = Ejc_1rhodpdn(0,l)
        Ejc_Lp(ns,l)       = Ejc_Lp(0,l)                  

      ENDDO   

c-----------------------------------------------------------------------
c     Writing the exact jump conditions
c-----------------------------------------------------------------------

400   format(1x,2000e16.6)
      open(unit=18,file='DAT/Jumps/Ejc_u.dat',status='unknown')     
      open(unit=28,file='DAT/Jumps/Ejc_dudn.dat',status='unknown')     
      open(unit=38,file='DAT/Jumps/Ejc_Lu.dat',status='unknown')   
      open(unit=48,file='DAT/Jumps/Ejc_v.dat',status='unknown')     
      open(unit=58,file='DAT/Jumps/Ejc_dvdn.dat',status='unknown')     
      open(unit=68,file='DAT/Jumps/Ejc_Lv.dat',status='unknown')   
      open(unit=78,file='DAT/Jumps/Ejc_p.dat',status='unknown')     
      open(unit=88,file='DAT/Jumps/Ejc_1rhodpdn.dat',status='unknown')  
      open(unit=98,file='DAT/Jumps/Ejc_Lp.dat',status='unknown')       
      
      do m=0,ns-1
        write(18,400)(Ejc_u(m,l),l=1,ms)
        write(28,400)(Ejc_dudn(m,l),l=1,ms)
        write(38,400)(Ejc_Lu(m,l),l=1,ms)
        write(48,400)(Ejc_v(m,l),l=1,ms)
        write(58,400)(Ejc_dvdn(m,l),l=1,ms)
        write(68,400)(Ejc_Lv(m,l),l=1,ms)
        write(78,400)(Ejc_p(m,l),l=1,ms)
        write(88,400)(Ejc_1rhodpdn(m,l),l=1,ms)
        write(98,400)(Ejc_Lp(m,l),l=1,ms)
      enddo
         
      close(18)
      close(28)
      close(38)
      close(48)
      close(58)
      close(68)
      close(78)
      close(88)
      close(98)     

c-----------------------------------------------------------------------
      
      return
      end

