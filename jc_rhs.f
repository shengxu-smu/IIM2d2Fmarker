!-----------------------------------------------------------------------
!     Given a function (in this case fy and fx) the program            !
!     calculates the derivative with respect to the variable           !
!     alpha.                                                           !
!                   dfy  = d/dalfa(fy)                                 !
!                   ddfy = d^2/dalpha^2(fy)                            !
!                   dfx  = J*[dp/dn])        			       !
!                                                                      !
!-----------------------------------------------------------------------  
c
      subroutine jc_rhs
      include 'parameter.inc'
      include 'surface.inc'
      real*8 foo
      real*8 tmp(2,0:ns),ccc(2,0:ns)
      real*8 gacobi
      

      DO l=1,ms

      do m=0,ns

        !=========================================!
        !   Modification for 2-fluids pressure    !
        !-----------------------------------------!
        ! We will find the derivatives:           !
        !  d[p]/d(alfa)     = dfy                 !
        !  d^2[p]/d(alfa^2) = ddfy                !
        !-----------------------------------------!

        tmp(1,m)=fy(m,l)
        tmp(2,m)=0.0d0
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,2)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,dfy(m,l),1)
      enddo
      dfy(ns,l)=dfy(0,l)
        !==========================================!
  
      do m=0,ns
        gacobi=sqrt(taox(m,l)*taox(m,l)+taoy(m,l)*taoy(m,l))
        dfx(m,l)=gacobi*Prin_jc(m,l,8)
      enddo

      do m=0,ns
        tmp(1,m)=dfx(m,l)
        tmp(2,m)=dfy(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,2)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,ddfx(m,l),1)
        call fdf(alfa(m),alfa(m+1),tmp(2,m),tmp(2,m+1),
     .           ccc(2,m),ccc(2,m+1),alfa(m),foo,ddfy(m,l),1)
      enddo
      ddfx(ns,l)=ddfx(0,l)
      ddfy(ns,l)=ddfy(0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
