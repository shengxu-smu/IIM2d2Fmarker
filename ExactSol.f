c
c     !****************************************************************!
c     !        ---  EXACT SOLUTION u,v,p FOR EACH POINT  ---           !
c     !                          Sept 25,2010                          !     
c     !****************************************************************!
c     !    In variables:                                               !
c     !                       - xi  (x coordinate)                     !  
c     !                       - yj  (y coordinate)                     !
c     !    Out variables:                                              !
c     !                       - au  (analytical u)                     !
c     !                       - av  (analytical v)                     !
c     !                       - a_p (analytical p)                     !
c     !                                                                !
c     !****************************************************************!
c
      subroutine  ExactSol(xi,yj,au,av,a_p)
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'

      real*8 xi,yj,au,av,a_p
      real*8 aa,bb,dd, Omega1,Omega2, x_c, y_c, rr
      real*8 Q, Apos, Bpos, Aneg, Bneg      
      integer ProblemOption
      
      ProblemOption = 2 

      !----------------------------------------------------------------!
      ! Option 1: Two cylinders and the interface (cuadratic behavior) !
      !----------------------------------------------------------------!
      ! It has dimensions but I only need this to check jump conditions

      if (ProblemOption==1) then
 
        aa        = 0.1d0   ! Inside radius
        bb        = 3.0d0   ! Outside radius
        dd        = 0.5d0   ! interface radius
        Omega1    = 1.0d0   ! Inside velocity
        Omega2    = 2.0d0   ! Outside velocity
        x_c       = 0.d0    ! center (x_c,y_c)  
        y_c       = 0.d0      

        Q    = 1/(mu_neg*(1/bb**2-1/dd**2)+ mu_pos*(1/dd**2-1/aa**2))
        Apos = Q*(mu_neg*(Omega1/bb**2-Omega2/dd**2) 
     .            + mu_pos*Omega2*(1/dd**2-1/aa**2)) 
        Bpos = Q*((Omega2-Omega1)*mu_neg)
        Aneg = Q*(mu_neg*Omega1*(1/bb**2-1/dd**2)
     .            +mu_pos*(Omega1/dd**2-Omega2/aa**2))
        Bneg = Q*((Omega2-Omega1)*mu_pos)

        rr = dsqrt((xi-x_c)**2 + (yj-y_c)**2)
 
        if (rr<dd) then
          au = -(Aneg + Bneg/(rr**2))*(yj-y_c) 
          av =  (Aneg + Bneg/(rr**2))*(xi-x_c)
          a_p =  rho_neg*(0.5d0*(Aneg*rr)**2 - 0.5d0*(Bneg/rr)**2
     .        + Aneg*Bneg*log(rr**2))
        else 
          au = -(Apos + Bpos/(rr**2))*(yj-y_c) 
          av =  (Apos + Bpos/(rr**2))*(xi-x_c)
          a_p = rho_pos*(0.5d0*(Apos*rr)**2 - 0.5d0*(Bpos/rr)**2
     .         + Apos*Bpos*log(rr**2))
        endif

      !----------------------------------------------------------------!
      ! Option 2: One cylinder and the interface (linear behavior)     !
      !----------------------------------------------------------------!      
      ! Dimensionless

      elseif (ProblemOption==2) then

        rr = dsqrt(xi**2 + yj**2)
 
        au = -yj
        av =  xi
        if (rr >= 0.5) then
          a_p=  0.5d0*(rr**2)*rho_pos 
        else 
          a_p = 0.5d0*(rr**2)*rho_neg
        endif
 
      !----------------------------------------------------------------!
      ! Option 2: One cylinder and the interface (linear behavior)     !
      !----------------------------------------------------------------!      

      elseif (ProblemOption==3) then

        x_c = 0.d0     
        y_c = 0.d0
        rr  = dsqrt((xi-x_c)**2 + (yj-y_c)**2)

        if (rr<dd) then
          au = (2.0d0*rr)**2*dsin(xi)*dsin(yj)
          av = 2.0d0*rr*sin(xi)*sin(yj);
          a_p = 0;
        else
          au = dsin(xi)*dsin(yj)
          av = dsin(xi)*dsin(yj)
          a_p = 0.0d0
        endif

      endif
     
      return
      end

c!*********************************************************************|

