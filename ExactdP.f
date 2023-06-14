c-----------------------------------------------------------------------
c
      subroutine  ExactdP(xi,yj,dpdx,dpdy)
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'

      real*8 xi,yj,dpdx,dpdy
      real*8 aa,bb,dd, Omega1,Omega2, xxc, yyc
      real*8 Q, Apos, Bpos, Aneg, Bneg      
      real*8 rr
      integer ProblemOption

      ProblemOption = 2

      !-----------------------------------!
      !            Parameters             !
      !-----------------------------------!

      aa        = 0.2d0   ! Inside radius
      bb        = 3.0d0   ! Outside radius
      dd        = 0.5d0   ! interface radius
      Omega1    = 1.0d0
      Omega2    = 2.0d0
      xxc       = 0.d0    ! center (xxc,yyc)  
      yyc       = 0.d0      

      !-----------------------------------!
      !  Option 1: Two cylinders          !
      !-----------------------------------!

      If (ProblemOption==1) then

      Q    = 1/(mu_neg*(1/bb**2-1/dd**2)+ mu_pos*(1/dd**2-1/aa**2))
      Apos = Q*(mu_neg*(Omega1/bb**2-Omega2/dd**2) 
     .          + mu_pos*Omega2*(1/dd**2-1/aa**2)) 
      Bpos = Q*((Omega2-Omega1)*mu_neg)
      Aneg = Q*(mu_neg*Omega1*(1/bb**2-1/dd**2)
     .          +mu_pos*(Omega1/dd**2-Omega2/aa**2))
      Bneg = Q*((Omega2-Omega1)*mu_pos)
      
      xi = xi - xxc
      yj = yj - yyc
      rr = dsqrt(xi**2+yj**2)
    
      dpdx = ((Apos**2)*xi+(Bpos**2)*xi/(rr**4)+2*Apos*Bpos*xi/(rr**2))*
     .        rho_pos
      dpdy = ((Apos**2)*yj+(Bpos**2)*yj/(rr**4)+2*Apos*Bpos*yj/(rr**2))*
     .        rho_pos
    
      !-----------------------------------!
      ! Option 2: One cylinder            !
      !-----------------------------------!

      elseif (ProblemOption==2) then

c      dpdx = (Omega2**2)*xi*rho_pos
c      dpdy = (Omega2**2)*yj*rho_pos

       dpdx = xi*rho_pos
       dpdy = yj*rho_pos
      
      !-----------------------------------!
      ! Option 3: p(X,Y)=sinX*sinY        !
      !-----------------------------------!

      elseif (ProblemOption==3) then
        
        dpdx = dcos(xi)*dsin(yj)
        dpdy = dsin(xi)*dcos(yj)
       
      endif

      return
      end

c-----------------------------------------------------------------------
