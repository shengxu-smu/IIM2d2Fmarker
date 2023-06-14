c  surface.inc
      integer jectob(ms),ncxe(ms),ncxc(ms),ncye(ms),ncyc(ms)
      integer ncxen(ms),ncxcn(ms),ncyen(ms),ncycn(ms)

      integer ixe(0:2*nx,ms),jexe(0:2*nx,ms),jcxe(0:2*nx,ms)
      integer ixc(0:2*nx,ms),jcxc(0:2*nx,ms),jexc(0:2*nx,ms)
      integer jye(0:2*ny,ms),ieye(0:2*ny,ms),icye(0:2*ny,ms)
      integer jyc(0:2*ny,ms),icyc(0:2*ny,ms),ieyc(0:2*ny,ms)

      integer ixen(0:2*nx,ms),jexen(0:2*nx,ms),jcxen(0:2*nx,ms)
      integer ixcn(0:2*nx,ms),jcxcn(0:2*nx,ms),jexcn(0:2*nx,ms)
      integer jyen(0:2*ny,ms),ieyen(0:2*ny,ms),icyen(0:2*ny,ms)
      integer jycn(0:2*ny,ms),icycn(0:2*ny,ms),ieycn(0:2*ny,ms)

      integer juyexe(0:2*nx,ms),juyexc(0:2*nx,ms),iuxcyc(0:2*ny,ms)
      integer ivxeye(0:2*ny,ms),jvycxe(0:2*nx,ms),jvycxc(0:2*nx,ms)

      real*8 a(ms),b(ms),xsc(ms),ysc(ms),theta(ms)
      real*8 xsct(ms),ysct(ms),thetat(ms)
      real*8 xsctt(ms),ysctt(ms),thetatt(ms),area(ms)
      real*8 alfa(0:ns),ra(0:4*ns,ms)
      real*8 xs(0:ns,ms),ys(0:ns,ms),xs0(0:ns,ms),ys0(0:ns,ms)
      real*8 taox(0:ns,ms),taoy(0:ns,ms),dtaox(0:ns,ms),dtaoy(0:ns,ms)
      real*8 us(0:ns,ms),vs(0:ns,ms),fx(0:ns,ms),fy(0:ns,ms)
      real*8 dfx(0:ns,ms),dfy(0:ns,ms),ddfx(0:ns,ms),ddfy(0:ns,ms)
      real*8 us0(0:ns,ms),vs0(0:ns,ms),fx0(0:ns,ms),fy0(0:ns,ms)
      real*8 dfx0(0:ns,ms),dfy0(0:ns,ms),ddfx0(0:ns,ms),ddfy0(0:ns,ms)

c-------------------------------------------------------
c    Parameters used in interface_velocityU.f
c-------------------------------------------------------   

      real*8 xsout(0:ns,ms), ysout(0:ns,ms) 
      real*8 xsin(0:ns,ms),  ysin(0:ns,ms)

c-------------------------------------------------------

      real*8 falfaxe(0:14,0:2*nx,ms),falfaxc(0:14,0:2*nx,ms)
      real*8 falfaye(0:14,0:2*ny,ms),falfayc(0:14,0:2*ny,ms)

      real*8 falfaxen(3,0:2*nx,ms),falfaxcn(3,0:2*nx,ms)
      real*8 falfayen(3,0:2*ny,ms),falfaycn(3,0:2*ny,ms)

      real*8 ujc(6,0:ns,ms),vjc(6,0:ns,ms),pjc(8,0:ns,ms)

      real*8 ujcxe(6,0:2*nx,ms),vjcxe(6,0:2*nx,ms)
      real*8 ujcye(6,0:2*ny,ms),vjcye(6,0:2*ny,ms)
      real*8 ujcxc(6,0:2*nx,ms),vjcxc(6,0:2*nx,ms)
      real*8 ujcyc(6,0:2*ny,ms),vjcyc(6,0:2*ny,ms)

      real*8 pjcxc(8,0:2*nx,ms),pjcyc(8,0:2*ny,ms)

      real*8 ujcxen(2,0:2*nx,ms),ujcxcn(2,0:2*nx,ms)
      real*8 vjcxen(2,0:2*nx,ms),vjcxcn(2,0:2*nx,ms)
      real*8 ujcyen(2,0:2*ny,ms),ujcycn(2,0:2*ny,ms)
      real*8 vjcyen(2,0:2*ny,ms),vjcycn(2,0:2*ny,ms)

      real*8 uciyexe(0:2*nx,ms),uciyexc(0:2*nx,ms),ucixcyc(0:2*ny,ms)
      real*8 vcixeye(0:2*ny,ms),vciycxe(0:2*nx,ms),vciycxc(0:2*nx,ms)

      real*8 ucudx(0:2*ny,ms),ucpdx(0:2*ny,ms),ucdy(0:2*nx,ms)
      real*8 ucdxx(2,0:2*ny,ms),ucdyy(2,0:2*nx,ms)
      real*8 vcdx(0:2*ny,ms),vcvdy(0:2*nx,ms),vcpdy(0:2*nx,ms)
      real*8 vcdxx(2,0:2*ny,ms),vcdyy(2,0:2*nx,ms)
      real*8 pcudx(0:2*ny,ms),pcvdx(0:2*ny,ms)
      real*8 pcudy(0:2*nx,ms),pcvdy(0:2*nx,ms)
      real*8 pcdxx(2,0:2*ny,ms),pcdyy(2,0:2*nx,ms)

      real*8 wsave(3*ns+15)

c-------------------------------------------------------
c    Definitions of the Principal jump conditions
c-------------------------------------------------------
   
       real*8 Prin_jc(0:ns,ms,9) 
c-------------------------------------------------------


      common/ncount/jectob,ncxe,ncxc,ncye,ncyc,ncxen,ncxcn,ncyen,ncycn
      common/index/ixe,jexe,jcxe,ixc,jcxc,jexc,
     .             jye,ieye,icye,jyc,icyc,ieyc,
     .             juyexe,juyexc,iuxcyc,ivxeye,jvycxe,jvycxc,
     .             ixen,jexen,jcxen,ixcn,jcxcn,jexcn,
     .             jyen,ieyen,icyen,jycn,icycn,ieycn
      common/surface/a,b,xsc,ysc,theta,xsct,ysct,thetat,
     .               xsctt,ysctt,thetatt,area,alfa,ra,
     .               xs,ys,xs0,ys0,taox,taoy,dtaox,dtaoy,us,vs,
     .               us0,vs0,falfaxe,falfaxc,falfaye,falfayc,
     .               falfaxen,falfaxcn,falfayen,falfaycn,
     .               xsout,ysout,xsin,ysin

      common/force/fx,fy,dfx,dfy,ddfx,ddfy,
     .             fx0,fy0,dfx0,dfy0,ddfx0,ddfy0
      common/jc/ujc,vjc,pjc,
     .          ujcxe,vjcxe,ujcye,vjcye,pjcxc,
     .          ujcxc,vjcxc,ujcyc,vjcyc,pjcyc,
     .          ujcxen,vjcxen,ujcyen,vjcyen,
     .          ujcxcn,vjcxcn,ujcycn,vjcycn,
     .          Prin_jc
      common/ci/uciyexe,uciyexc,ucixcyc,vcixeye,vciycxe,vciycxc
      common/cd/ucudx,ucpdx,ucdy,ucdxx,ucdyy,
     .          vcdx,vcvdy,vcpdy,vcdxx,vcdyy,
     .          pcudx,pcvdx,pcudy,pcvdy,pcdxx,pcdyy
      common/fftsave/wsave
