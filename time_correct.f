c-----------------------------------------------------------------------
c
      subroutine time_correct(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,io,jo,ie,je,id,jd,is,js,ii,jj
      real*8 xa(2),ya(2),foo
      real*8 fac,ta,tb,ts
      real*8 uxjc,uyjc,utjc,vxjc,vyjc,vtjc,c1,c2,c3,c4

      ta=t-dt
      tb=ta+fac*dt

      DO l=1,ms

      do io=0,ncxen(l)
        ii=ixen(io,l)
        call index_search(io,ii,ie,l,ixe,2*nx,ms,ncxe(l))
        if(ie.ne.-1) then
          jd=jcxe(ie,l)-jcxen(io,l)
          if(jd.ge.1) then
            do js=jcxen(io,l)+1,jcxe(ie,l)
              xa(1)=falfaxen(1,io,l)
              xa(2)=falfaxe(2,ie,l)
              ya(1)=ta
              ya(2)=tb
              call interpolate(xa,ya,2,yc(js),ts,foo)
              ya(1)=ujcxen(1,io,l)
              ya(2)=ujcxe(1,ie,l)
              call interpolate(xa,ya,2,yc(js),uxjc,foo)
              ya(1)=ujcxen(2,io,l)
              ya(2)=ujcxe(2,ie,l)
              call interpolate(xa,ya,2,yc(js),uyjc,foo)
              utjc=(falfaxen(2,io,l)*uxjc+falfaxen(3,io,l)*uyjc)
              if(krk.eq.1) then
                u(ii,js)=u(ii,js)-(-utjc*(tb-ts))
              elseif(krk.eq.2) then
                u(ii,js)=u(ii,js)-(-utjc*(ta-ts))
              elseif(krk.eq.3) then
                if(ts.ge.0.5d0*(ta+tb)) then
                  u(ii,js)=u(ii,js)-(-utjc*(tb-ts))
                else
                  u(ii,js)=u(ii,js)-(-utjc*(ta-ts))
                endif
              elseif(krk.eq.4) then
                c1=-dt1*(utjc*(tb-ts))
                c4=-dt1*(utjc*(ta-ts))
                if(ts.ge.0.5*(ta+tb)) then
                  c2=-dt1*(utjc*(tb-ts))
                else
                  c2=-dt1*(utjc*(ta-ts))
                endif
                c3=c2
                u(ii,js)=u(ii,js)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
              endif
            enddo
          elseif(jd.le.-1) then
            do js=jcxen(io,l),jcxe(ie,l)+1,-1
              xa(1)=falfaxen(1,io,l)
              xa(2)=falfaxe(2,ie,l)
              ya(1)=ta
              ya(2)=tb
              call interpolate(xa,ya,2,yc(js),ts,foo)
              ya(1)=ujcxen(1,io,l)
              ya(2)=ujcxe(1,ie,l)
              call interpolate(xa,ya,2,yc(js),uxjc,foo)
              ya(1)=ujcxen(2,io,l)
              ya(2)=ujcxe(2,ie,l)
              call interpolate(xa,ya,2,yc(js),uyjc,foo)
              utjc=-(falfaxen(2,io,l)*uxjc+falfaxen(3,io,l)*uyjc)
              if(krk.eq.1) then
                u(ii,js)=u(ii,js)-(-utjc*(tb-ts))
              elseif(krk.eq.2) then
                u(ii,js)=u(ii,js)-(-utjc*(ta-ts))
              elseif(krk.eq.3) then
                if(ts.ge.0.5d0*(ta+tb)) then
                  u(ii,js)=u(ii,js)-(-utjc*(tb-ts))
                else
                  u(ii,js)=u(ii,js)-(-utjc*(ta-ts))
                endif
              elseif(krk.eq.4) then
                c1=-dt1*(utjc*(tb-ts))
                c4=-dt1*(utjc*(ta-ts))
                if(ts.ge.0.5d0*(ta+tb)) then
                  c2=-dt1*(utjc*(tb-ts))
                else
                  c2=-dt1*(utjc*(ta-ts))
                endif
                c3=c2
                u(ii,js)=u(ii,js)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
              endif
            enddo
          endif
        endif
      enddo
c
      do jo=0,ncycn(l)
        jj=jycn(jo,l)
        call index_search(jo,jj,je,l,jyc,2*ny,ms,ncyc(l))
        if(je.ne.-1) then
          id=ieyc(je,l)-ieycn(jo,l)
          if(id.ge.1) then
            do is=ieycn(jo,l)+1,ieyc(je,l)
              call index_search(0,is,io,l,ixen,2*nx,ms,ncxen(l))
              call index_search(0,is,ie,l,ixe,2*nx,ms,ncxe(l))
              if(io.eq.-1.or.ie.eq.-1) then
                xa(1)=falfaycn(1,jo,l)
                xa(2)=falfayc(1,je,l)
                ya(1)=ta
                ya(2)=tb
                call interpolate(xa,ya,2,xe(is),ts,foo)
                ya(1)=ujcycn(1,jo,l)
                ya(2)=ujcyc(1,je,l)
                call interpolate(xa,ya,2,xe(is),uxjc,foo)
                ya(1)=ujcycn(2,jo,l)
                ya(2)=ujcyc(2,je,l)
                call interpolate(xa,ya,2,xe(is),uyjc,foo)
                utjc=(falfaycn(2,jo,l)*uxjc+falfaycn(3,jo,l)*uyjc)
                if(krk.eq.1) then
                  u(is,jj)=u(is,jj)-(-utjc*(tb-ts))
                elseif(krk.eq.2) then
                  u(is,jj)=u(is,jj)-(-utjc*(ta-ts))
                elseif(krk.eq.3) then
                  if(ts.ge.0.5d0*(ta+tb)) then
                    u(is,jj)=u(is,jj)-(-utjc*(tb-ts))
                  else
                    u(is,jj)=u(is,jj)-(-utjc*(ta-ts))
                  endif
                elseif(krk.eq.4) then
                  c1=-dt1*(utjc*(tb-ts))
                  c4=-dt1*(utjc*(ta-ts))
                  if(ts.ge.0.5d0*(ta+tb)) then
                    c2=-dt1*(utjc*(tb-ts))
                  else
                    c2=-dt1*(utjc*(ta-ts))
                  endif
                  c3=c2
                  u(is,jj)=u(is,jj)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
                endif
              endif
            enddo
          elseif(id.le.-1) then
            do is=ieycn(jo,l),ieyc(je,l)+1,-1
              call index_search(0,is,io,l,ixen,2*nx,ms,ncxen)
              call index_search(0,is,ie,l,ixe,2*nx,ms,ncxe)
              if(io.eq.-1.or.ie.eq.-1) then
                xa(1)=falfaycn(1,jo,l)
                xa(2)=falfayc(1,je,l)
                ya(1)=ta
                ya(2)=tb
                call interpolate(xa,ya,2,xe(is),ts,foo)
                ya(1)=ujcycn(1,jo,l)
                ya(2)=ujcyc(1,je,l)
                call interpolate(xa,ya,2,xe(is),uxjc,foo)
                ya(1)=ujcycn(2,jo,l)
                ya(2)=ujcyc(2,je,l)
                call interpolate(xa,ya,2,xe(is),uyjc,foo)
                utjc=-(falfaycn(2,jo,l)*uxjc+falfaycn(3,jo,l)*uyjc)
                if(krk.eq.1) then
                  u(is,jj)=u(is,jj)-(-utjc*(tb-ts))
                elseif(krk.eq.2) then
                  u(is,jj)=u(is,jj)-(-utjc*(ta-ts))
                elseif(krk.eq.3) then
                  if(ts.ge.0.5d0*(ta+tb)) then
                    u(is,jj)=u(is,jj)-(-utjc*(tb-ts))
                  else
                    u(is,jj)=u(is,jj)-(-utjc*(ta-ts))
                  endif
                elseif(krk.eq.4) then
                  c1=-dt1*(utjc*(tb-ts))
                  c4=-dt1*(utjc*(ta-ts))
                  if(ts.ge.0.5d0*(ta+tb)) then
                    c2=-dt1*(utjc*(tb-ts))
                  else
                    c2=-dt1*(utjc*(ta-ts))
                  endif
                  c3=c2
                  u(is,jj)=u(is,jj)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
                endif
              endif
            enddo
          endif
        endif
      enddo
c
      do io=0,ncxcn(l)
        ii=ixcn(io,l)
        call index_search(io,ii,ie,l,ixc,2*nx,ms,ncxc(l))
        if(ie.ne.-1) then
          jd=jexc(ie,l)-jexcn(io,l)
          if(jd.ge.1) then
            do js=jexcn(io,l)+1,jexc(ie,l)
              xa(1)=falfaxcn(1,io,l)
              xa(2)=falfaxc(2,ie,l)
              ya(1)=ta
              ya(2)=tb
              call interpolate(xa,ya,2,ye(js),ts,foo)
              ya(1)=vjcxcn(1,io,l)
              ya(2)=vjcxc(1,ie,l)
              call interpolate(xa,ya,2,ye(js),vxjc,foo)
              ya(1)=vjcxcn(2,io,l)
              ya(2)=vjcxc(2,ie,l)
              call interpolate(xa,ya,2,ye(js),vyjc,foo)
              vtjc=(falfaxcn(2,io,l)*vxjc+falfaxcn(3,io,l)*vyjc)
              if(krk.eq.1) then
                v(ii,js)=v(ii,js)-(-vtjc*(tb-ts))
              elseif(krk.eq.2) then
                v(ii,js)=v(ii,js)-(-vtjc*(ta-ts))
              elseif(krk.eq.3) then
                if(ts.ge.0.5d0*(ta+tb)) then
                  v(ii,js)=v(ii,js)-(-vtjc*(tb-ts))
                else
                  v(ii,js)=v(ii,js)-(-vtjc*(ta-ts))
                endif
              elseif(krk.eq.4) then
                c1=-dt1*(vtjc*(tb-ts))
                c4=-dt1*(vtjc*(ta-ts))
                if(ts.ge.0.5d0*(ta+tb)) then
                  c2=-dt1*(vtjc*(tb-ts))
                else
                  c2=-dt1*(vtjc*(ta-ts))
                endif
                c3=c2
                v(ii,js)=v(ii,js)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
              endif
            enddo
          elseif(jd.le.-1) then
            do js=jexcn(io,l),jexc(ie,l)+1,-1
              xa(1)=falfaxcn(1,io,l)
              xa(2)=falfaxc(2,ie,l)
              ya(1)=ta
              ya(2)=tb
              call interpolate(xa,ya,2,ye(js),ts,foo)
              ya(1)=vjcxcn(1,io,l)
              ya(2)=vjcxc(1,ie,l)
              call interpolate(xa,ya,2,ye(js),vxjc,foo)
              ya(1)=vjcxcn(2,io,l)
              ya(2)=vjcxc(2,ie,l)
              call interpolate(xa,ya,2,ye(js),vyjc,foo)
              vtjc=-(falfaxcn(2,io,l)*vxjc+falfaxcn(3,io,l)*vyjc)
              if(krk.eq.1) then
                v(ii,js)=v(ii,js)-(-vtjc*(tb-ts))
              elseif(krk.eq.2) then
                v(ii,js)=v(ii,js)-(-vtjc*(ta-ts))
              elseif(krk.eq.3) then
                if(ts.ge.0.5d0*(ta+tb)) then
                  v(ii,js)=v(ii,js)-(-vtjc*(tb-ts))
                else
                  v(ii,js)=v(ii,js)-(-vtjc*(ta-ts))
                endif
              elseif(krk.eq.4) then
                c1=-dt1*(vtjc*(tb-ts))
                c4=-dt1*(vtjc*(ta-ts))
                if(ts.ge.0.5d0*(ta+tb)) then
                  c2=-dt1*(vtjc*(tb-ts))
                else
                  c2=-dt1*(vtjc*(ta-ts))
                endif
                c3=c2
                v(ii,js)=v(ii,js)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
              endif
            enddo
          endif
        endif
      enddo
c
      do jo=0,ncyen(l)
        jj=jyen(jo,l)
        call index_search(jo,jj,je,l,jye,2*ny,ms,ncye(l))
        if(je.ne.-1) then
          id=icye(je,l)-icyen(jo,l)
          if(id.ge.1) then
            do is=icyen(jo,l)+1,icye(je,l)
              call index_search(0,is,io,l,ixcn,2*nx,ms,ncxcn(l))
              call index_search(0,is,ie,l,ixc,2*nx,ms,ncxc(l))
              if(io.eq.-1.or.ie.eq.-1) then
                xa(1)=falfayen(1,jo,l)
                xa(2)=falfaye(1,je,l)
                ya(1)=ta
                ya(2)=tb
                call interpolate(xa,ya,2,xc(is),ts,foo)
                ya(1)=vjcyen(1,jo,l)
                ya(2)=vjcye(1,je,l)
                call interpolate(xa,ya,2,xc(is),vxjc,foo)
                ya(1)=vjcyen(2,jo,l)
                ya(2)=vjcye(2,je,l)
                call interpolate(xa,ya,2,xc(is),vyjc,foo)
                vtjc=(falfayen(2,jo,l)*vxjc+falfayen(3,jo,l)*vyjc)
                if(krk.eq.1) then
                  v(is,jj)=v(is,jj)-(-vtjc*(tb-ts))
                elseif(krk.eq.2) then
                  v(is,jj)=v(is,jj)-(-vtjc*(ta-ts))
                elseif(krk.eq.3) then
                  if(ts.ge.0.5d0*(ta+tb)) then
                    v(is,jj)=v(is,jj)-(-vtjc*(tb-ts))
                  else
                    v(is,jj)=v(is,jj)-(-vtjc*(ta-ts))
                  endif
                elseif(krk.eq.4) then
                  c1=-dt1*(vtjc*(tb-ts))
                  c4=-dt1*(vtjc*(ta-ts))
                  if(ts.ge.0.5d0*(ta+tb)) then
                    c2=-dt1*(vtjc*(tb-ts))
                  else
                    c2=-dt1*(vtjc*(ta-ts))
                  endif
                  c3=c2
                  v(is,jj)=v(is,jj)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
                endif
              endif
            enddo
          elseif(id.le.-1) then
            do is=icyen(jo,l),icye(je,l)+1,-1
              call index_search(0,is,io,l,ixcn,2*nx,ms,ncxcn(l))
              call index_search(0,is,ie,l,ixc,2*nx,ms,ncxc(l))
              if(io.eq.-1.or.ie.eq.-1) then
                xa(1)=falfayen(1,jo,l)
                xa(2)=falfaye(1,je,l)
                ya(1)=ta
                ya(2)=tb
                call interpolate(xa,ya,2,xc(is),ts,foo)
                ya(1)=vjcyen(1,jo,l)
                ya(2)=vjcye(1,je,l)
                call interpolate(xa,ya,2,xc(is),vxjc,foo)
                ya(1)=vjcyen(2,jo,l)
                ya(2)=vjcye(2,je,l)
                call interpolate(xa,ya,2,xc(is),vyjc,foo)
                vtjc=-(falfayen(2,jo,l)*vxjc+falfayen(3,jo,l)*vyjc)
                if(krk.eq.1) then
                  v(is,jj)=v(is,jj)-(-vtjc*(tb-ts))
                elseif(krk.eq.2) then
                  v(is,jj)=v(is,jj)-(-vtjc*(ta-ts))
                elseif(krk.eq.3) then
                  if(ts.ge.0.5d0*(ta+tb)) then
                    v(is,jj)=v(is,jj)-(-vtjc*(tb-ts))
                  else
                    v(is,jj)=v(is,jj)-(-vtjc*(ta-ts))
                  endif
                elseif(krk.eq.4) then
                  c1=-dt1*(vtjc*(tb-ts))
                  c4=-dt1*(vtjc*(ta-ts))
                  if(ts.ge.0.5d0*(ta+tb)) then
                    c2=-dt1*(vtjc*(tb-ts))
                  else
                    c2=-dt1*(vtjc*(ta-ts))
                  endif
                  c3=c2
                  v(is,jj)=v(is,jj)-(dt/6.0d0)*(c1+2.0d0*(c2+c3)+c4)
                endif
              endif
            enddo
          endif
        endif
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
