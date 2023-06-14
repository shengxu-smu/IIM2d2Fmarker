c field.inc
      real*8 x(-2:2*nx+2),y(-2:2*ny+2)
      real*8 xe(-1:nx+1),ye(-1:ny+1),xc(0:nx+1),yc(0:ny+1)
      real*8 wsavei(3*nx+15),wsavej(3*ny+15)

      real*8 u(0:nx,0:ny+1),v(0:nx+1,0:ny),p(0:nx+1,0:ny+1),
     .       ap(0:nx+1,0:ny+1)
      real*8 funHp(0:nx+1,0:ny+1),funHu(0:nx+1,0:ny+1)
      real*8 funHv(0:nx+1,0:ny+1)  
      real*8 d(0:nx+1,0:ny+1),o(0:nx+1,0:ny+1),ph(0:nx+1,0:ny+1)
      real*8 uee(0:nx,0:ny),uce(nx,0:ny),ucc(nx,0:ny+1)
      real*8 vee(0:nx,0:ny),vec(0:nx,ny),vcc(0:nx+1,ny)

      common/grid/x,y,xe,ye,xc,yc
      common/field/u,v,p,ap,d,o,ph,uee,uce,ucc,vee,vec,vcc,funHp,
     .             funHu,funHv
      common/fftfield/wsavei,wsavej
