c-----------------------------------------------------------------------
c
      subroutine pbc_fft(pw,pe,ps,pn)
      include 'parameter.inc'
      include 'field.inc'
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      real*8 o1,o2,r1,r2,a,b,foo,rr
      real*8 nada,nada1,pww

c  periodic
      do j=1,ny
        pw(j)=0.0d0
        pe(j)=0.0d0
      enddo
      do i=1,nx
        ps(i)=0.0d0
        pn(i)=0.0d0
      enddo

c  dirichlet 
      if(lw_pbc.eq.1) then
        do j=1,ny
          pw(j)=-(0.125d0*(yc(j)*yc(j)-xc(1)*xc(1))+
     .             0.25d0*(yc(j)*xc(1)))/
     .           (yc(j)*yc(j)+xc(1)*xc(1))**2.0d0
          pw(j)=0.0d0
        enddo
      endif

      if(le_pbc.eq.1) then
        do j=1,ny
          pe(j)=0.0d0
        enddo
      endif

      if(ls_pbc.eq.1) then
        do i=1,nx
          ps(i)=0.0d0
        enddo
      endif

      if(ln_pbc.eq.1) then
        do i=1,nx
          pn(i)=0.0d0
        enddo
      endif

c  newmann
      if(lw_pbc.eq.2) then
        do j=1,ny
          pw(j)=((dx2/Re_pos)*(ucc(2,j)-ucc(1,j)
     .         +(dx/dy)*(vee(0,j)-vee(0,j-1)))
     .         +(dy2/Re_pos)*(ucc(1,j+1)-2.0d0*ucc(1,j)+ucc(1,j-1))
     .         -dx1*(u(1,j)*u(1,j)-u(0,j)*u(0,j))
     .         -dy1*(uce(1,j)*v(1,j)-uce(1,j-1)*v(1,j-1)))

           !====================================!
           !     -- for 2-fluid pressure --     !
           !------------------------------------!

c            call ExactdP(xc(1),yc(j),pw(j),nada)
           !====================================!
        enddo
      endif

      if(le_pbc.eq.2) then
        do j=1,ny
          pe(j)=((dx2/Re_pos)*(ucc(nx-1,j)-ucc(nx,j)
     .         -(dx/dy)*(vee(nx,j)-vee(nx,j-1)))
     .         +(dy2/Re_pos)*(ucc(nx,j+1)-2.0d0*ucc(nx,j)+ucc(nx,j-1))
     .         -dx1*(u(nx,j)*u(nx,j)-u(nx-1,j)*u(nx-1,j))
     .         -dy1*(uce(nx,j)*v(nx,j)-uce(nx,j-1)*v(nx,j-1)))

           !====================================!
           !    -- for 2-fluid pressure --      !
           !------------------------------------!

c            call ExactdP(xc(nx),yc(j),pe(j),nada)

           !====================================!

        enddo
      endif

      if(ls_pbc.eq.2) then
        do i=1,nx
          ps(i)=((dy2/Re_pos)*(vcc(i,2)-vcc(i,1)
     .         +(dy/dx)*(uee(i,0)-uee(i-1,0)))
     .         +(dx2/Re_pos)*(vcc(i+1,1)-2.0d0*vcc(i,1)+vcc(i-1,1))
     .         -dy1*(v(i,1)*v(i,1)-v(i,0)*v(i,0))
     .         -dx1*(vec(i,1)*u(i,1)-vec(i-1,1)*u(i-1,1)))-1.0d0
           
           !====================================!
           !    -- for 2-fluid pressure --      !
           !------------------------------------!
          
c           call ExactdP(xc(i),yc(1),nada,ps(i))          

           !====================================!
        enddo
      endif

      if(ln_pbc.eq.2) then
        do i=1,nx
          pn(i)=((dy2/Re_pos)*(vcc(i,ny-1)-vcc(i,ny)
     .         -(dy/dx)*(uee(i,ny)-uee(i-1,ny)))
     .         +(dx2/Re_pos)*(vcc(i+1,ny)-2.0d0*vcc(i,ny)+vcc(i-1,ny))
     .         -dy1*(v(i,ny)*v(i,ny)-v(i,ny-1)*v(i,ny-1))
     .         -dx1*(vec(i,ny)*u(i,ny)-vec(i-1,ny)*u(i-1,ny)))-1.0d0

           !====================================!
           !    -- for 2-fluid pressure --      !
           !------------------------------------!

c            call ExactdP(xc(i),yc(ny),nada,pn(i))

           !====================================!
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
