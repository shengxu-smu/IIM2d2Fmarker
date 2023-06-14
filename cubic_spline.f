c-----------------------------------------------------------------------
c
      subroutine cubic_spline(xx,yy,cc,ns,nn,num)
      integer ns,nn,num
      real*8 dm,dn
      real*8 xx(0:nn),yy(num,0:nn),cc(num,0:nn)
      real*8 h(0:ns-1),dyy(num,0:ns-1)
      real*8 a(ns),b(ns),d(num,ns)

c  coefficient matrix
c
      do m=0,ns-1
        h(m)=xx(m+1)-xx(m)
        do n=1,num
          dyy(n,m)=yy(n,m+1)-yy(n,m)
        enddo
      enddo

      do m=1,ns-1
        a(m)=h(m-1)/(h(m)+h(m-1))
        b(m)=1.0-a(m)
        do n=1,num
          d(n,m)=6.0d0*(dyy(n,m)/h(m)-dyy(n,m-1)/h(m-1))/(h(m)+h(m-1))
        enddo
      enddo
      a(ns)=h(0)/(h(0)+h(ns-1))
      b(ns)=1.0d0-a(ns)
      do n=1,num
        d(n,ns)=6.0d0*(dyy(n,0)/h(0)-dyy(n,ns-1)/h(ns-1))/(h(0)+h(ns-1))
      enddo

      a(1)=a(1)/2.0d0
      b(1)=b(1)/2.0d0
      do n=1,num
        d(n,1)=d(n,1)/2.0d0
      enddo
      b(ns)=b(ns)/2.0d0
      a(ns)=a(ns)/2.0d0
      do n=1,num
        d(n,ns)=d(n,ns)/2.0d0
      enddo
 
c  elimination to upper triangular matrix
c
c      
c |  1                                        -1       |
c | a(1)    1     b(1)                                 |
c |  0     a(2)    2      b(2)      		       |
c |						       |
c |                      ...		    	       |
c |						       |
c |                             a(ns-1)   2   b(ns-1)  |
c |  0     b(ns)                         a(ns)   1     |
c
c
c  ======>
c
c
c |  1                                        -1       |
c |  0      1     b(1)                       a(1)      |
c |  0      0       1      b(2)              a(2)      |
c |                                                    |
c |                      ...                           |
c |                               1   b(ns-2) a(ns-2)  |
c |                            a(ns-1)   2    b(ns-1)  |
c |  0                          b(ns)   a(ns)    1     |
c
c
c  ======>
c
c
c |  1                                        -1       |
c |  0      1     b(1)                       a(1)      |
c |  0      0       1      b(2)              a(2)      |
c |                                                    |
c |                      ...                           |
c |                               1   b(ns-2) a(ns-2)  |
c |                               0     1    a(ns-1)   |
c |  0                            0     0      1       |
c
c
      do m=2,ns-2
        dm=2.0d0-b(m-1)*a(m)
        do n=1,num
          d(n,m)=(d(n,m)-d(n,m-1)*a(m))/dm
        enddo
        a(m)=-a(m-1)*a(m)/dm
        b(m)=b(m)/dm

        dn=1.0d0-a(m-1)*b(ns)
        do n=1,num
          d(n,ns)=(d(n,ns)-d(n,m-1)*b(ns))/dn
        enddo
        b(ns)=-b(m-1)*b(ns)/dn
        a(ns)=a(ns)/dn
      enddo

      dm=2.0d0-b(ns-2)*a(ns-1)
      do n=1,num
        d(n,ns-1)=(d(n,ns-1)-d(n,ns-2)*a(ns-1))/dm
      enddo
      a(ns-1)=(b(ns-1)-a(ns-2)*a(ns-1))/dm

      dn=1.0d0-a(ns-2)*b(ns)
      do n=1,num
        d(n,ns)=(d(n,ns)-d(n,ns-2)*b(ns))/dn
      enddo
      a(ns)=(a(ns)-b(ns-2)*b(ns))/dn

      dn=1.0d0-a(ns-1)*a(ns)
      do n=1,num
        d(n,ns)=(d(n,ns)-d(n,ns-1)*a(ns))/dn
      enddo

c  back substitution
c
      do n=1,num
        cc(n,ns)=d(n,ns)
        cc(n,ns-1)=d(n,ns-1)-a(ns-1)*cc(n,ns)
        do m=ns-2,1,-1
          cc(n,m)=d(n,m)-a(m)*cc(n,ns)-b(m)*cc(n,m+1)
        enddo
        cc(n,0)=cc(n,ns)
      enddo
 
      return
      end


c-----------------------------------------------------------------------
