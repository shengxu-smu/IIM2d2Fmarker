c-----------------------------------------------------------------------
c
      subroutine surface_property
      include 'parameter.inc'
      include 'surface.inc'
      real*8 foo
      real*8 coord(2,0:ns),ccc(2,0:ns)

      !call StepFun
     
      DO l=1,ms

        do m=0,ns
          coord(1,m)=xs(m,l)
          coord(2,m)=ys(m,l)
        enddo
        call cubic_spline(alfa,coord,ccc,ns,ns,2)
        do m=0,ns-1
          call fdf(alfa(m),alfa(m+1),coord(1,m),coord(1,m+1),
     .             ccc(1,m),ccc(1,m+1),alfa(m),foo,taox(m,l),1)
          call fdf(alfa(m),alfa(m+1),coord(2,m),coord(2,m+1),
     .             ccc(2,m),ccc(2,m+1),alfa(m),foo,taoy(m,l),1)
        enddo
        taox(ns,l)=taox(0,l)
        taoy(ns,l)=taoy(0,l)

        do m=0,ns
          coord(1,m)=taox(m,l)
          coord(2,m)=taoy(m,l)
        enddo
        call cubic_spline(alfa,coord,ccc,ns,ns,2)
        do m=0,ns-1
          call fdf(alfa(m),alfa(m+1),coord(1,m),coord(1,m+1),
     .             ccc(1,m),ccc(1,m+1),alfa(m),foo,dtaox(m,l),1)
          call fdf(alfa(m),alfa(m+1),coord(2,m),coord(2,m+1),
     .             ccc(2,m),ccc(2,m+1),alfa(m),foo,dtaoy(m,l),1)
        enddo
        dtaox(ns,l)=dtaox(0,l)
        dtaoy(ns,l)=dtaoy(0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
