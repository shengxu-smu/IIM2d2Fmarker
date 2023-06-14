      subroutine gmres_pressure(fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
c      include "mkl_rci.fi"
      integer izes,nknown,ks,itercount,irc_request
      parameter (izes=128,nknown=kks*ms)
      integer ipar(izes)
      real*8 dvar,xt,yt,gacobi,xx0,xx1,yy0,yy1,cc0,cc1,foo
      real*8 dpar(izes)
      real*8 tmp(nknown*(2*nknown+1)+(nknown*(nknown+9))/2+1)
      real*8 rhs(nknown),rhsb(nknown),solution(nknown),residual(nknown)
      real*8 xx(0:kks),yy(1,0:kks),cc(1,0:kks)
      real*8 r(1,ns),w(1,ns)
      real*8 fac
c
c---------------------------------------------------------------------------
c External BLAS functions are taken from MKL BLAS to use
c with the RCI (P)FGMRES solver
c---------------------------------------------------------------------------
      real*8 dnrm2
      external dnrm2

c---------------------------------------------------------------------------
c Initialize the right hand side
c---------------------------------------------------------------------------
  
      DO l=1,ms
      ks=(l-1)*kks

      do k=1,kks
        m=kips+(k-1)*kips
        xt=taox(m,l)
        yt=taoy(m,l)
        gacobi=dsqrt(xt*xt+yt*yt)
        rhs(ks+k)=dfx0(m,l)
        ! Initial guess [p]
        solution(ks+k)= fy0(m,l)
      enddo
 
      ENDDO

c---------------------------------------------------------------------------
c Save the right-hand side in vector rhsb for future use
c---------------------------------------------------------------------------
      call dcopy(nknown,rhs,1,rhsb,1)

c---------------------------------------------------------------------------
c Initialize the solver
c---------------------------------------------------------------------------
      call dfgmres_init(nknown,solution,rhs,irc_request,ipar,dpar,tmp)
      if(irc_request.ne.0) goto 999

c---------------------------------------------------------------------------
c Set the desired parameters:
c do the restart after 5 iterations
c LOGICAL parameters:
c do not do the stopping test for the maximal number of iterations
c do not do the Preconditioned iterations of FGMRES method
c DOUBLE PRECISION parameters
c set the relative tolerance to 1.0d-3 instead of default value 1.0d-6
c---------------------------------------------------------------------------
      ipar(15)=5
      ipar(5)=100
      ipar(8)=1
      ipar(11)=0
      dpar(1)=1.0d-3
c
c---------------------------------------------------------------------------
c Check the correctness and consistency of the newly set parameters
c---------------------------------------------------------------------------
      call dfgmres_check(nknown,solution,rhs,irc_request,ipar,dpar,tmp)
      if(irc_request.ne.0) goto 999

c---------------------------------------------------------------------------
c Compute the solution by RCI (P)FGMRES solver with
c Reverse Communication starts here
c---------------------------------------------------------------------------
1     call dfgmres(nknown,solution,rhs,irc_request,ipar,dpar,tmp)

c---------------------------------------------------------------------------
c If irc_request=0, then the solution was found with the required precision
c---------------------------------------------------------------------------
      if(irc_request.eq.0) goto 3

c---------------------------------------------------------------------------
c If irc_request=1, then compute the matrix-vector prodct A*tmp(ipar(22))
c and put the result in vector tmp(ipar(23))
c---------------------------------------------------------------------------
      if(irc_request.eq.1) then

      DO l=1,ms

      ks=(l-1)*kks
      do k=1,kks
        m=kips+(k-1)*kips
        xx(k)=alfa(m)
        yy(1,k)=tmp(ipar(22)+ks+k-1)
      enddo
      xx(0)=xx(kks)-sl
      yy(1,0)=yy(1,kks)
      call cubic_spline(xx,yy,cc,kks,kks,1)

      do m=0,ns-1
        k=int(m/kips)
        xx0=xx(k)
        yy0=yy(1,k)
        cc0=cc(1,k)
        xx1=xx(k+1)
        yy1=yy(1,k+1)
        cc1=cc(1,k+1)
        call fdf(xx0,xx1,yy0,yy1,cc0,cc1,alfa(m),fy0(m,l),foo,1)
      enddo
      fy0(ns,l)=fy0(0,l)

      do m=0,ns-1
        r(1,m+1)=fy0(m,l)
      enddo
      call vrfftf(1,ns,r,w,1,wsave)
      do m=ns/2,ns
c        r(1,m)=0.0d0
      enddo
      call vrfftb(1,ns,r,w,1,wsave)
      do m=0,ns-1
        fy0(m,l)=r(1,m+1)
      enddo
      fy0(ns,l)=fy0(0,l)

      ENDDO

      call mvproduct_gmres(fac)

      DO l=1,ms
      ks=(l-1)*kks

      do k=1,kks
        m=kips+(k-1)*kips
        xt=taox(m,l)
        yt=taoy(m,l)
        gacobi=dsqrt(xt*xt+yt*yt)
        tmp(ipar(23)+ks+k-1)=dfx0(m,l)
      enddo

      ENDDO

      goto 1

      endif

c---------------------------------------------------------------------------
c If irc_request=2, then do the user-defined stopping test
c The residual stopping test for the computed solution is performed here
c---------------------------------------------------------------------------
c NOTE: from this point vector rhsb is no longer containing the right-hand
c side of the problem! It contains the current FGMRES approximation to the
c solution. If you need to keep the right-hand side, save it in some other
c vector before the call to DFGMRES routine. Here we saved it in vector
c rhs. The vector rhsb is used instead of rhs to preserve the original
c right-hand side of the problem and guarantee the proper restart of FGMRES
c method. Vector rhsb will be altered when computing the residual stopping
c criterion!
c---------------------------------------------------------------------------
      if(irc_request.eq.2) then
c Request to the DFGMRES_GET routine to put the solution into rhsb via ipar(13)
        ipar(13)=1
c Get the current FGMRES solution in the vector rhsb
        call dfgmres_get(nknown,solution,rhsb,irc_request,ipar,
     .                   dpar,tmp,itercount)
c Compute the current true residual
      DO l=1,ms

      ks=(l-1)*kks
      do k=1,kks
        m=kips+(k-1)*kips
        xx(k)=alfa(m)
        yy(1,k)=rhsb(ks+k)
      enddo
      xx(0)=xx(kks)-sl
      yy(1,0)=yy(1,kks)
      call cubic_spline(xx,yy,cc,kks,kks,1)

      do m=0,ns-1
        k=int(m/kips)
        xx0=xx(k)
        yy0=yy(1,k)
        cc0=cc(1,k)
        xx1=xx(k+1)
        yy1=yy(1,k+1)
        cc1=cc(1,k+1)
        call fdf(xx0,xx1,yy0,yy1,cc0,cc1,alfa(m),fy0(m,l),foo,1)
      enddo
      fy0(ns,l)=fy0(0,l)

      do m=0,ns-1
        r(1,m+1)=fy0(m,l)
      enddo
      call vrfftf(1,ns,r,w,1,wsave)
      do m=ns/2,ns
c        r(1,m)=0.0d0
      enddo
      call vrfftb(1,ns,r,w,1,wsave)
      do m=0,ns-1
        fy0(m,l)=r(1,m+1)
      enddo
      fy0(ns,l)=fy0(0,l)

      ENDDO
      
      !=================================
      !write(*,*),'Second fy0=',fy0(0:8,1) 
      !================================= 

      call mvproduct_gmres(fac)

      DO l=1,ms
      ks=(l-1)*kks

      do k=1,kks
        m=kips+(k-1)*kips
        xt=taox(m,l)
        yt=taoy(m,l)
        gacobi=dsqrt(xt*xt+yt*yt)
        residual(ks+k)=dfx0(m,l)
      enddo

      ENDDO

      call daxpy(nknown,-1.0d0,rhs,1,residual,1)
      dvar=dnrm2(nknown,residual,1)
      if(dvar.lt.1.0d-3) then
        goto 3
      else
        goto 1
      endif

      endif

C---------------------------------------------------------------------------
C If irc_request=4, then check if the norm of the next generated vector is
C not zero up to rounding and computational errors. The norm is contained
C in dpar(7) parameter
C---------------------------------------------------------------------------
      if(irc_request.eq.4) then
        if(dpar(7).lt.1.0d-12) then
          goto 3
        else
          goto 1
        endif
      else
c---------------------------------------------------------------------------
c If irc_request=anything else, then DFGMRES subroutine failed
c to compute the solution vector: solution(nknown)
c---------------------------------------------------------------------------
        goto 999
      endif
c
c---------------------------------------------------------------------------
c Reverse Communication ends here
c Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
c call DFGMRES_GET routine as computed_solution is still containing
c the initial guess!). Request to DFGMRES_GET to put the solution into
c vector solution(nknown) via ipar(13)
c---------------------------------------------------------------------------
3     ipar(13)=0
      call dfgmres_get(nknown,solution,rhs,irc_request,ipar,
     .                 dpar,tmp,itercount)
      print *,' Number of GMRES iterations: ',itercount
      goto 1000

999   print *,'The solver has returned the ERROR code ', irc_request

1000  continue

      DO l=1,ms

      ks=(l-1)*kks
      do k=1,kks
        m=kips+(k-1)*kips
        xx(k)=alfa(m)
        yy(1,k)=solution(ks+k)
      enddo
      xx(0)=xx(kks)-sl
      yy(1,0)=yy(1,kks)
      call cubic_spline(xx,yy,cc,kks,kks,1)

      do m=0,ns-1
        k=int(m/kips)
        xx0=xx(k)
        yy0=yy(1,k)
        cc0=cc(1,k)
        xx1=xx(k+1)
        yy1=yy(1,k+1)
        cc1=cc(1,k+1)
        call fdf(xx0,xx1,yy0,yy1,cc0,cc1,alfa(m),fy0(m,l),foo,1)
      enddo
      fy0(ns,l)=fy0(0,l)

      do m=0,ns-1
        r(1,m+1)=fy0(m,l)
      enddo
      call vrfftf(1,ns,r,w,1,wsave)
      do m=ns/2,ns
c        r(1,m)=0.0d0
      enddo
      call vrfftb(1,ns,r,w,1,wsave)
      do m=0,ns-1
        fy0(m,l)=r(1,m+1)
      enddo
      fy0(ns,l)=fy0(0,l)

      ENDDO
c-----------------------------------------------------------------------

      END
