c
      !****************************************************************!
      !          --- CALCULATE THE PRINCIPAL JUMP CONDITIONS ---       !
      !                            Sept,2010                           !
      !****************************************************************!
      !                                                                !
      !  In variables:                                                 !
      !                       - Uxsys                                  !
      !                       - Vxsys                                  !
      !                       - dudnp   - dudnm                        !
      !                       - dvdnp   - dvdnm                        !
      !                       - ddudnp  - ddudnm                       !
      !                       - ddvdnp  - ddvdnm                       !
      !                                                                !
      !   Remark:   We write the solution in the data files:           !
      !                       - Njc_u.dat                              !
      !                       - Njc_dudn.dat                           !
      !                       - Njc_Lu.dat                             !
      !                       - Njc_v.dat                              !
      !                       - Njc_dvdn.dat                           !
      !                       - Njc_Lv.dat                             !
      !                       - Njc_p.dat                              !
      !                       - Njc_1rhodpdn.dat                       !
      !                       - Njc_Lp.dat                             !
      !                                                                !
      !****************************************************************!
c
      subroutine Principal_Jumps(Uxsys,Vxsys,dudnp,dudnm,dvdnp,dvdnm,
     .                            ddudnp,ddudnm,ddvdnp,ddvdnm)
      
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      
      real*8 r(1,ns),w(1,ns)

      real*8 Uxsys(0:ns,ms), Vxsys(0:ns,ms)
      real*8 dudnp(0:ns,ms), dudnm(0:ns,ms)
      real*8 dvdnp(0:ns,ms), dvdnm(0:ns,ms)
      real*8 ddudnp(0:ns,ms), ddudnm(0:ns,ms)
      real*8 ddvdnp(0:ns,ms), ddvdnm(0:ns,ms)
      
      real*8 xt,yt,xn,yn,gacobi
      real*8 coord(2,0:ns),ccc(2,0:ns)
      real*8 coord1(1,0:ns),cc1(1,0:ns)
      real*8 dualfa, dvalfa 
      real*8 dutau(0:ns,ms), dvtau(0:ns,ms)
      real*8 mup, mum, rhop, rhom
      real*8 dfun1alfa,foo,determ
      real*8 fun1(0:ns,ms)    
      real*8 dudxp, dudyp, dudxm, dudym  
      real*8 dvdxp, dvdyp, dvdxm, dvdym   
      
      real*8 Jump_dudn, Jump_dvdn
      real*8 Jump_mududn, Jump_mudvdn
      real*8 Jump_Lu, Jump_Lv
      real*8 Jump_p
      real*8 Jump_nududn, Jump_nudvdn
      real*8 Jump_Laplace_p
      real*8 xnp, xnm, ynp, ynm
      real*8 dtx, dty
      real*8 con1,con2
      real*8 alphaa(0:ns)
      real*8 Kcur,forceft,forcefn

      mup  = mu_pos
      mum  = mu_neg
      rhop = rho_pos      
      rhom = rho_neg

      DO l=1,ms   
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !               ---   JUMP CONDITIONS  ---                    !
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

          !------------------------------------------!
          !           du/dtau and dv/dtau            !
          !------------------------------------------!
  
        do m=0,ns
          coord(1,m) = Uxsys(m,l)
          coord(2,m) = Vxsys(m,l)
          alphaa(m)= alfa(m)
        enddo

        call cubic_spline(alphaa,coord,ccc,ns,ns,2)
        do m=0,ns-1

          !  Unitary tangential and normal vector 
          xt = taox(m,l)
          yt = taoy(m,l)
          dtx= dtaox(m,l)
          dty= dtaoy(m,l)
          xn = yt 
          yn =-xt
          gacobi = dsqrt(xn*xn+yn*yn)
          xt = xt/gacobi
          yt = yt/gacobi
          xn = xn/gacobi
          yn = yn/gacobi
          dtx= dtx/(gacobi**2)
          dty= dty/(gacobi**2)

          !------------------------------------------!
          !       curvature & forces: fn,ft          !
          !------------------------------------------!

          Kcur = -(dty*xt-dtx*yt)
          forcefn = (1.0d0/Bo)*Kcur
          forceft = 0.0d0

          !------------------------------------------!
          !              dutau & dvtau               !
          !------------------------------------------!

          call fdf(alphaa(m),alphaa(m+1),coord(1,m),coord(1,m+1),
     .             ccc(1,m),ccc(1,m+1),alphaa(m),foo,dualfa,1)
          
          call fdf(alphaa(m),alphaa(m+1),coord(2,m),coord(2,m+1),
     .             ccc(2,m),ccc(2,m+1),alphaa(m),foo,dvalfa,1)
         
          dutau(m,l)=dualfa/gacobi
          dvtau(m,l)=dvalfa/gacobi
       
          !------------------------------------------!
          !        [mu*du/dn] and [mu*dv/dn]         !
          !------------------------------------------!
  
          Jump_mududn = -(mup-mum)*
     .                   ((dutau(m,l)*xn+dvtau(m,l)*yn)*xt
     .                   +(dutau(m,l)*xt+dvtau(m,l)*yt)*xn )
          Jump_mudvdn = -(mup-mum)*
     .                   ((dutau(m,l)*xn+dvtau(m,l)*yn)*yt
     .                   +(dutau(m,l)*xt+dvtau(m,l)*yt)*yn )
          
          !------------------------------------------!
          !           [du/dn] and [dv/dn]            !
          !------------------------------------------!
         
          if (mup>mum) then
            Jump_dudn = (1.0d0/mup)*(Jump_mududn - (mup-mum)*dudnm(m,l))
            Jump_dvdn = (1.0d0/mup)*(Jump_mudvdn - (mup-mum)*dvdnm(m,l))
          else
            Jump_dudn = (1.0d0/mum)*(Jump_mududn - (mup-mum)*dudnp(m,l))
            Jump_dvdn = (1.0d0/mum)*(Jump_mudvdn - (mup-mum)*dvdnp(m,l))
          endif   
         
          !------------------------------------------!
          !             [Lu] and [Lv]                !
          !------------------------------------------!          
            
          Jump_Lu = ddudnp(m,l)-ddudnm(m,l)+Kcur*Jump_dudn
          Jump_Lv = ddvdnp(m,l)-ddvdnm(m,l)+Kcur*Jump_dvdn

          !------------------------------------------!
          !                  [p]                     !
          !------------------------------------------!
         
          Jump_p =forcefn-2.0d0*(mup-mum)*(dutau(m,l)*xt+dvtau(m,l)*yt)

          !------------------------------------------!
          !       [nu*dudn]  and [nu*dvdn]           !
          !------------------------------------------!
 
          if (mum<mup) then
            Jump_nududn =(1.0d0/rhop)*Jump_mududn
     .                   -(mum)*(rhop-rhom)/(rhop*rhom)*dudnm(m,l)
            Jump_nudvdn =(1.0d0/rhop)*Jump_mudvdn
     .                   -(mum)*(rhop-rhom)/(rhop*rhom)*dvdnm(m,l)
          else
            Jump_nududn =(1.0d0/rhom)*Jump_mududn
     .                   -(mup)*(rhop-rhom)/(rhop*rhom)*dudnp(m,l)
            Jump_nudvdn =(1.0d0/rhom)*Jump_mudvdn
     .                   -(mup)*(rhop-rhom)/(rhop*rhom)*dvdnp(m,l)
          endif

          !------------------------------------------!
          !               Making fun1                !
          !------------------------------------------!

          fun1(m,l) = - xt*Jump_nududn-yt*Jump_nudvdn 
     .                +(mup/rhop-mum/rhom)*(dutau(m,l)*xn+dvtau(m,l)*yn)
          
          !------------------------------------------!
          !                 Jumps                    !
          !------------------------------------------!

          Prin_jc(m,l,1) = 0.0d0
          Prin_jc(m,l,2) = Jump_dudn
          Prin_jc(m,l,3) = Jump_Lu
          Prin_jc(m,l,4) = 0.0d0
          Prin_jc(m,l,5) = Jump_dvdn
          Prin_jc(m,l,6) = Jump_Lv
          Prin_jc(m,l,7) = Jump_p

        enddo

        Prin_jc(ns,l,1) = Prin_jc(0,l,1)
        Prin_jc(ns,l,2) = Prin_jc(0,l,2)
        Prin_jc(ns,l,3) = Prin_jc(0,l,3)
        Prin_jc(ns,l,4) = Prin_jc(0,l,4)
        Prin_jc(ns,l,5) = Prin_jc(0,l,5)
        Prin_jc(ns,l,6) = Prin_jc(0,l,6)
        Prin_jc(ns,l,7) = Prin_jc(0,l,7)
        fun1(ns,l)      = fun1(0,l)

          !------------------------------------------!
          !       [1/rho*dp/dn] = d(fun1)/dalfa      !
          !------------------------------------------!

        do m=0,ns
          coord1(1,m) = fun1(m,l)
        enddo
        call cubic_spline(alphaa,coord1,cc1,ns,ns,1)
        do m=0,ns-1
          call fdf(alphaa(m),alphaa(m+1),coord1(1,m),coord1(1,m+1),
     .             cc1(1,m),cc1(1,m+1),alphaa(m),foo,dfun1alfa,1)
          
          Prin_jc(m,l,8) = dfun1alfa
          
          !------------------------------------------!
          !              [Laplace p]                 !
          !------------------------------------------!

          ! Finding du/dx+ du/dx- and the others
         
          xt = taox(m,l)
          yt = taoy(m,l)
          xn = yt
          yn = -xt
          gacobi = dsqrt(xt*xt+yt*yt)
          xt  = xt/gacobi
          yt  = yt/gacobi 
          xn = xn/gacobi
          yn = yn/gacobi

          determ = xt*yn-yt*xn
          dudxp = (1/determ)*(yn*dutau(m,l)-yt*dudnp(m,l)) 
          dudyp = (1/determ)*(xt*dudnp(m,l)-xn*dutau(m,l))

          dudxm = (1/determ)*(yn*dutau(m,l)-yt*dudnm(m,l))
          dudym = (1/determ)*(xt*dudnm(m,l)-xn*dutau(m,l))
         
          dvdxp = (1/determ)*(yn*dvtau(m,l)-yt*dvdnp(m,l)) 
          dvdyp = (1/determ)*(xt*dvdnp(m,l)-xn*dvtau(m,l))
       
          dvdxm = (1/determ)*(yn*dvtau(m,l)-yt*dvdnm(m,l))
          dvdym = (1/determ)*(xt*dvdnm(m,l)-xn*dvtau(m,l))  

          Jump_Laplace_p  = 2*(dudxp*dvdyp-dudyp*dvdxp- 
     .                        (dudxm*dvdym-dudym*dvdxm))
          Prin_jc(m,l,9) = Jump_Laplace_p

        enddo

        Prin_jc(ns,l,8) = Prin_jc(0,l,8)
        Prin_jc(ns,l,9) = Prin_jc(0,l,9)
  
      ENDDO 

      DO l=1,ms
        do m=0,ns-1
          r(1,m+1)=Prin_jc(m,l,8)
        enddo
        call vrfftf(1,ns,r,w,1,wsave)
        do m=ns/2,ns
          r(1,m)=0.0d0
        enddo
        call vrfftb(1,ns,r,w,1,wsave)
        do m=0,ns-1
          Prin_jc(m,l,8)=r(1,m+1)
        enddo
        Prin_jc(ns,l,8)=Prin_jc(0,l,8)

      ENDDO


c-----------------------------------------------------------------------
c     Writing the exact jump conditions
c-----------------------------------------------------------------------

400   format(1x,2000e16.6)

      open(unit=18,file='DAT/Jumps/Njc_u.dat',status='unknown')
      open(unit=28,file='DAT/Jumps/Njc_dudn.dat',status='unknown')
      open(unit=38,file='DAT/Jumps/Njc_Lu.dat',status='unknown')
      open(unit=48,file='DAT/Jumps/Njc_v.dat',status='unknown')
      open(unit=58,file='DAT/Jumps/Njc_dvdn.dat',status='unknown')
      open(unit=68,file='DAT/Jumps/Njc_Lv.dat',status='unknown')
      open(unit=78,file='DAT/Jumps/Njc_p.dat',status='unknown')
      open(unit=88,file='DAT/Jumps/Njc_1rhodpdn.dat',status='unknown')
      open(unit=98,file='DAT/Jumps/Njc_Lp.dat',status='unknown')

      do m=0,ns-1
        write(18,400)(Prin_jc(m,l,1),l=1,ms)
        write(28,400)(Prin_jc(m,l,2),l=1,ms)
        write(38,400)(Prin_jc(m,l,3),l=1,ms)
        write(48,400)(Prin_jc(m,l,4),l=1,ms)
        write(58,400)(Prin_jc(m,l,5),l=1,ms)
        write(68,400)(Prin_jc(m,l,6),l=1,ms)
        write(78,400)(Prin_jc(m,l,7),l=1,ms)
        write(88,400)(Prin_jc(m,l,8),l=1,ms)
        write(98,400)(Prin_jc(m,l,9),l=1,ms)
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
