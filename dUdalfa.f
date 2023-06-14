c-----------------------------------------------------------------------
c
      subroutine dUdalfa(u_Lag,v_Lag,du_dtau,dv_dtau)
      include 'parameter.inc'
      include 'surface.inc'
      real*8 foo
      real*8 coord(2,0:ns),ccc(2,0:ns)
      real*8 u_Lag(0:ns), v_Lag(0:ns)
      real*8 du_dtau(0:ns), dv_dtau(0:ns)
       
        do m=0,ns
          coord(1,m) = u_Lag(m)
          coord(2,m) = v_Lag(m)
        enddo
        call cubic_spline(alfa,coord,ccc,ns,ns,2)
        do m=0,ns-1
          call fdf(alfa(m),alfa(m+1),coord(1,m),coord(1,m+1),
     .             ccc(1,m),ccc(1,m+1),alfa(m),foo,du_dtau(m),1)
          call fdf(alfa(m),alfa(m+1),coord(2,m),coord(2,m+1),
     .             ccc(2,m),ccc(2,m+1),alfa(m),foo,dv_dtau(m),1)
        enddo
        du_dtau(ns)=du_dtau(0)
        dv_dtau(ns)=dv_dtau(0)

      return
      end
c
c-----------------------------------------------------------------------
