c  parameter.inc
      implicit none

      integer i,j,k,l,m,n,nstart,nend,nstep,nx,ny,ms,ns,kips,kks

      integer icfl,iread,iwrite,iplot,itout,ianimation
      integer isingular,move,ifxy,ipfft,ichorin,itemporal
      integer lw_ubc,le_ubc,ls_ubc,ln_ubc
      integer lw_vbc,le_vbc,ls_vbc,ln_vbc
      integer lw_pbc,le_pbc,ls_pbc,ln_pbc
      integer NumN

      real*8 pi,cflc,cflv,dt,dtcfl,dt0,dt1,dt2,t,tout,t0
      real*8 x0,y0,xl,yl,sl,dx,dy,dalfa,dar,dx1,dx2,dy1,dy2,hdx,hdy
      real*8 mu_pos,mu_neg,rho_pos,rho_neg,Lrepre,Urepre,rho0
      real*8 mupos,muneg,rhopos,rhoneg,Re_pos,Re_neg
      real*8 Rb,Rey,rhor,mur,Bo 
 
      parameter(pi=
     .3.1415926535897932384626433832795028841971693993751058209749446d0)

      !---------------------------------------------------------!
      !         --- Dimensionless parameters ---                !
      !---------------------------------------------------------!

      parameter(Bo      = 2.5d1)
      parameter(Re_pos  = 3.5d1)           
      parameter(Re_neg  = 3.5d0)           
      parameter(rho_pos = 1.0d0)           
      parameter(rho_neg = 1.0d-3)
      parameter(mu_pos  = rho_pos/Re_pos)
      parameter(mu_neg  = rho_neg/Re_neg)

      !---------------------------------------------------------!

      parameter(xl=1.0d0,yl=2.0d0,sl=2.0d0*pi)
      parameter(x0=-0.5d0,y0=-0.5d0)
      parameter(NumN=64)
      parameter(nx=NumN+1,ny=NumN+1,ms=1,ns=4*NumN)
      parameter(kips=1,kks=ns/kips)

      parameter(isingular=1,move=1,ifxy=1,ipfft=1,ichorin=0,itemporal=0)
      parameter(lw_ubc=1,le_ubc=1,ls_ubc=1,ln_ubc=1)
      parameter(lw_vbc=2,le_vbc=2,ls_vbc=1,ln_vbc=1)
      parameter(lw_pbc=2,le_pbc=2,ls_pbc=2,ln_pbc=2)

      common/run1/nstep,nstart,nend,
     .           iread,iwrite,itout,iplot,ianimation,icfl
      common/run2/cflc,cflv,dt,dtcfl,dt0,dt1,dt2,t,tout,t0,
     .           dx,dy,dx1,dy1,dx2,dy2,hdx,hdy,dalfa,dar
