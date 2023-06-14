c-----------------------------------------------------------------------
c
      subroutine rhs(krk)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer krk,i1,j1
      real*8 conv,grad,visc
      real*8 rr, Reu, Rev, Rep
      real*8 Reone, Retwo
      real*8 signx, signy, maxx,dd,xx,yy
      real*8 forcegx,forcegy

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                ---  Body forces  ---                  !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      forcegx = 0.0d0
      forcegy = -1.0d0 

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                  ---  rhs_u  ---                      !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      do j=1,ny
        do i=1,nx-1
          Reu = Re_neg*funHu(i,j)+Re_pos*(1.0d0-funHu(i,j))
          conv=-dx1*(ucc(i+1,j)*ucc(i+1,j)-ucc(i,j)*ucc(i,j))
     .         -dy1*(uee(i,j)*vee(i,j)-uee(i,j-1)*vee(i,j-1))
          visc=(dx2)*(u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))+
     .         (dy2)*(u(i,j+1)-2.0d0*u(i,j)+u(i,j-1))
          urk(i,j,krk) = conv + (1/Reu)*visc + forcegx
        enddo
      enddo

      !  Adding the jump contributions

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxe(l)
            urk(ixe(i,l),jexe(i,l)+1,krk)=urk(ixe(i,l),jexe(i,l)+1,krk)
     .                                   -ucdy(i,l)
            i1 = ixe(i,l)
            j1 = jcxe(i,l)
            Reone = Re_neg*funHu(i1,j1)+Re_pos*(1.0d0-funHu(i1,j1))
            urk(ixe(i,l),jcxe(i,l),krk)=urk(ixe(i,l),jcxe(i,l),krk)
     .                                  +ucdyy(1,i,l)*(1/Reone)
            i1 = ixe(i,l)
            j1 = jcxe(i,l)+1
            Retwo = Re_neg*funHu(i1,j1)+Re_pos*(1.0d0-funHu(i1,j1))
            urk(ixe(i,l),jcxe(i,l)+1,krk)=urk(ixe(i,l),jcxe(i,l)+1,krk)
     .                                    +ucdyy(2,i,l)*(1/Retwo)
          enddo

          do j=0,ncyc(l)
            i1 = icyc(j,l)
            j1 = jyc(j,l)
            urk(icyc(j,l),jyc(j,l),krk)=urk(icyc(j,l),jyc(j,l),krk)
     .                                -ucudx(j,l)-ucpdx(j,l)
            i1 = ieyc(j,l)
            j1 = jyc(j,l)
            Reone = Re_neg*funHu(i1,j1)+Re_pos*(1.0d0-funHu(i1,j1))
            urk(ieyc(j,l),jyc(j,l),krk)=urk(ieyc(j,l),jyc(j,l),krk)
     .                                  +ucdxx(1,j,l)*(1/Reone)
            i1 = ieyc(j,l)+1
            j1 = jyc(j,l)
            Retwo = Re_neg*funHu(i1,j1)+Re_pos*(1.0d0-funHu(i1,j1))
            urk(ieyc(j,l)+1,jyc(j,l),krk)=urk(ieyc(j,l)+1,jyc(j,l),krk)
     .                                    +ucdxx(2,j,l)*(1/Retwo)
          enddo
        enddo
      endif

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww! 
      !                  ---  rhs_v  ---                      !    
      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      do j=1,ny-1
        do i=1,nx
          Rev = Re_neg*funHv(i,j) + Re_pos*(1.0d0-funHv(i,j))
          conv=-dy1*(vcc(i,j+1)*vcc(i,j+1)-vcc(i,j)*vcc(i,j))
     .         -dx1*(uee(i,j)*vee(i,j)-uee(i-1,j)*vee(i-1,j))
          visc=(dx2)*(v(i+1,j)-2.0d0*v(i,j)+v(i-1,j))+
     .         (dy2)*(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))
          vrk(i,j,krk) = conv + (1/Rev)*visc + forcegy
        enddo
      enddo

      !  Adding the jump contributions

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            i1 = ixc(i,l)
            j1 = jcxc(i,l)
            vrk(ixc(i,l),jcxc(i,l),krk)=vrk(ixc(i,l),jcxc(i,l),krk)
     .                                 -vcvdy(i,l)-vcpdy(i,l)
            i1 = ixc(i,l)
            j1 = jexc(i,l)
            Reone = Re_neg*funHv(i1,j1)+Re_pos*(1.0d0-funHv(i1,j1))
            vrk(ixc(i,l),jexc(i,l),krk)=vrk(ixc(i,l),jexc(i,l),krk)
     .                                  +vcdyy(1,i,l)*(1/Reone)
            i1 = ixc(i,l)
            j1 = jexc(i,l)+1
            Retwo = Re_neg*funHv(i1,j1)+Re_pos*(1.0d0-funHv(i1,j1))
            vrk(ixc(i,l),jexc(i,l)+1,krk)=vrk(ixc(i,l),jexc(i,l)+1,krk)
     .                                  +vcdyy(2,i,l)*(1/Retwo)
          enddo

          do j=0,ncye(l)
            vrk(ieye(j,l)+1,jye(j,l),krk)=vrk(ieye(j,l)+1,jye(j,l),krk)
     .                                    -vcdx(j,l)
            i1 = icye(j,l)
            j1 = jye(j,l)
            Reone = Re_neg*funHv(i1,j1)+Re_pos*(1.0d0-funHv(i1,j1))
            vrk(icye(j,l),jye(j,l),krk)=vrk(icye(j,l),jye(j,l),krk)
     .                                  +vcdxx(1,j,l)*(1/Reone)
            i1 = icye(j,l)+1
            j1 = jye(j,l)
            Retwo = Re_neg*funHv(i1,j1)+Re_pos*(1.0d0-funHv(i1,j1))
            vrk(icye(j,l)+1,jye(j,l),krk)=vrk(icye(j,l)+1,jye(j,l),krk)
     .                                    +vcdxx(2,j,l)*(1/Retwo)
          enddo
        enddo
      endif
      
      return
      end


c-----------------------------------------------------------------------
