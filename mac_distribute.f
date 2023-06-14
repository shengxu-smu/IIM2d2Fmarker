c-----------------------------------------------------------------------
c
      subroutine mac_distribute
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      real*8 alfa0,alfa1,alpha,signx,signy
      real*8 aa(22),bb(22),f(22)
      
      DO l=1,ms

      do i=0,ncxe(l)
        alpha=falfaxe(0,i,l)
        m=int(alpha/dalfa)
        if(m.eq.-1) m=0
        if(m.eq.ns) m=ns-1

        alfa0=alfa(m)
        alfa1=alfa(m+1)
        do n=1,6
          aa(n+2)=(ujc(n,m+1,l)-ujc(n,m,l))/dalfa
          bb(n+2)=(ujc(n,m,l)*alfa1-ujc(n,m+1,l)*alfa0)/dalfa
          aa(n+8)=(vjc(n,m+1,l)-vjc(n,m,l))/dalfa
          bb(n+8)=(vjc(n,m,l)*alfa1-vjc(n,m+1,l)*alfa0)/dalfa
        enddo

        do n=3,14
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,6
          ujcxe(n,i,l)=f(n+2)
          vjcxe(n,i,l)=f(n+8)
        enddo
      enddo
       
      do i=0,ncxc(l) 
        alpha=falfaxc(0,i,l)
        m=int(alpha/dalfa)
        if(m.eq.-1) m=0
        if(m.eq.ns) m=ns-1

        alfa0=alfa(m)
        alfa1=alfa(m+1)
        do n=1,6
          aa(n+2)=(ujc(n,m+1,l)-ujc(n,m,l))/dalfa
          bb(n+2)=(ujc(n,m,l)*alfa1-ujc(n,m+1,l)*alfa0)/dalfa
          aa(n+8)=(vjc(n,m+1,l)-vjc(n,m,l))/dalfa
          bb(n+8)=(vjc(n,m,l)*alfa1-vjc(n,m+1,l)*alfa0)/dalfa
          aa(n+14)=(pjc(n,m+1,l)-pjc(n,m,l))/dalfa
          bb(n+14)=(pjc(n,m,l)*alfa1-pjc(n,m+1,l)*alfa0)/dalfa
        enddo
        aa(21)=(pjc(7,m+1,l)-pjc(7,m,l))/dalfa
        bb(21)=(pjc(7,m,l)*alfa1-pjc(7,m+1,l)*alfa0)/dalfa
        aa(22)=(pjc(8,m+1,l)-pjc(8,m,l))/dalfa
        bb(22)=(pjc(8,m,l)*alfa1-pjc(8,m+1,l)*alfa0)/dalfa

        do n=3,22
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,6
          ujcxc(n,i,l)=f(n+2)
          vjcxc(n,i,l)=f(n+8)
          pjcxc(n,i,l)=f(n+14)
        enddo
        pjcxc(7,i,l)=f(21)
        pjcxc(8,i,l)=f(22)
      enddo

      do j=0,ncye(l)
        alpha=falfaye(0,j,l)
        m=int(alpha/dalfa)
        if(m.eq.-1) m=0
        if(m.eq.ns) m=ns-1

        alfa0=alfa(m)
        alfa1=alfa(m+1)
        do n=1,6
          aa(n+2)=(ujc(n,m+1,l)-ujc(n,m,l))/dalfa
          bb(n+2)=(ujc(n,m,l)*alfa1-ujc(n,m+1,l)*alfa0)/dalfa
          aa(n+8)=(vjc(n,m+1,l)-vjc(n,m,l))/dalfa
          bb(n+8)=(vjc(n,m,l)*alfa1-vjc(n,m+1,l)*alfa0)/dalfa
        enddo

        do n=3,14
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,6
          ujcye(n,j,l)=f(n+2)
          vjcye(n,j,l)=f(n+8)
        enddo
      enddo

      do j=0,ncyc(l)
        alpha=falfayc(0,j,l)
        m=int(alpha/dalfa)
        if(m.eq.-1) m=0
        if(m.eq.ns) m=ns-1

        alfa0=alfa(m)
        alfa1=alfa(m+1)
        do n=1,6
          aa(n+2)=(ujc(n,m+1,l)-ujc(n,m,l))/dalfa
          bb(n+2)=(ujc(n,m,l)*alfa1-ujc(n,m+1,l)*alfa0)/dalfa
          aa(n+8)=(vjc(n,m+1,l)-vjc(n,m,l))/dalfa
          bb(n+8)=(vjc(n,m,l)*alfa1-vjc(n,m+1,l)*alfa0)/dalfa
          aa(n+14)=(pjc(n,m+1,l)-pjc(n,m,l))/dalfa
          bb(n+14)=(pjc(n,m,l)*alfa1-pjc(n,m+1,l)*alfa0)/dalfa
        enddo
        aa(21)=(pjc(7,m+1,l)-pjc(7,m,l))/dalfa
        bb(21)=(pjc(7,m,l)*alfa1-pjc(7,m+1,l)*alfa0)/dalfa
        aa(22)=(pjc(8,m+1,l)-pjc(8,m,l))/dalfa
        bb(22)=(pjc(8,m,l)*alfa1-pjc(8,m+1,l)*alfa0)/dalfa

        do n=3,22
          f(n)=aa(n)*alpha+bb(n)
        enddo
        do n=1,6
          ujcyc(n,j,l)=f(n+2)
          vjcyc(n,j,l)=f(n+8)
          pjcyc(n,j,l)=f(n+14)
        enddo
        pjcyc(7,j,l)=f(21)
        pjcyc(8,j,l)=f(22)
      enddo

      ENDDO

      DO l=1,ms

      do i=0,ncxe(l)
        signx=falfaxe(2,i,l)
        signy=falfaxe(3,i,l)

        ujcxe(1,i,l)=signx*ujcxe(1,i,l)
        vjcxe(1,i,l)=signx*vjcxe(1,i,l)
        ujcxe(3,i,l)=signx*ujcxe(3,i,l)
        vjcxe(3,i,l)=signx*vjcxe(3,i,l)
        ujcxe(5,i,l)=signx*ujcxe(5,i,l)
        vjcxe(5,i,l)=signx*vjcxe(5,i,l)

        ujcxe(2,i,l)=signy*ujcxe(2,i,l)
        vjcxe(2,i,l)=signy*vjcxe(2,i,l)
        ujcxe(4,i,l)=signy*ujcxe(4,i,l)
        vjcxe(4,i,l)=signy*vjcxe(4,i,l)
        ujcxe(6,i,l)=signy*ujcxe(6,i,l)
        vjcxe(6,i,l)=signy*vjcxe(6,i,l)
      enddo

      do i=0,ncxc(l)
        signx=falfaxc(2,i,l)
        signy=falfaxc(3,i,l)

        ujcxc(1,i,l)=signx*ujcxc(1,i,l)
        vjcxc(1,i,l)=signx*vjcxc(1,i,l)
        ujcxc(3,i,l)=signx*ujcxc(3,i,l)
        vjcxc(3,i,l)=signx*vjcxc(3,i,l)
        ujcxc(5,i,l)=signx*ujcxc(5,i,l)
        vjcxc(5,i,l)=signx*vjcxc(5,i,l)

        ujcxc(2,i,l)=signy*ujcxc(2,i,l)
        vjcxc(2,i,l)=signy*vjcxc(2,i,l)
        ujcxc(4,i,l)=signy*ujcxc(4,i,l)
        vjcxc(4,i,l)=signy*vjcxc(4,i,l)
        ujcxc(6,i,l)=signy*ujcxc(6,i,l)
        vjcxc(6,i,l)=signy*vjcxc(6,i,l)
      enddo

      do j=0,ncye(l)
        signx=falfaye(2,j,l)
        signy=falfaye(3,j,l)

        ujcye(1,j,l)=signx*ujcye(1,j,l)
        vjcye(1,j,l)=signx*vjcye(1,j,l)
        ujcye(3,j,l)=signx*ujcye(3,j,l)
        vjcye(3,j,l)=signx*vjcye(3,j,l)
        ujcye(5,j,l)=signx*ujcye(5,j,l)
        vjcye(5,j,l)=signx*vjcye(5,j,l)

        ujcye(2,j,l)=signy*ujcye(2,j,l)
        vjcye(2,j,l)=signy*vjcye(2,j,l)
        ujcye(4,j,l)=signy*ujcye(4,j,l)
        vjcye(4,j,l)=signy*vjcye(4,j,l)
        ujcye(6,j,l)=signy*ujcye(6,j,l)
        vjcye(6,j,l)=signy*vjcye(6,j,l)
      enddo

      do j=0,ncyc(l)
        signx=falfayc(2,j,l)
        signy=falfayc(3,j,l)

        ujcyc(1,j,l)=signx*ujcyc(1,j,l)
        vjcyc(1,j,l)=signx*vjcyc(1,j,l)
        ujcyc(3,j,l)=signx*ujcyc(3,j,l)
        vjcyc(3,j,l)=signx*vjcyc(3,j,l)
        ujcyc(5,j,l)=signx*ujcyc(5,j,l)
        vjcyc(5,j,l)=signx*vjcyc(5,j,l)

        ujcyc(2,j,l)=signy*ujcyc(2,j,l)
        vjcyc(2,j,l)=signy*vjcyc(2,j,l)
        ujcyc(4,j,l)=signy*ujcyc(4,j,l)
        vjcyc(4,j,l)=signy*vjcyc(4,j,l)
        ujcyc(6,j,l)=signy*ujcyc(5,j,l)
        vjcyc(6,j,l)=signy*vjcyc(6,j,l)
      enddo

      ENDDO


      return
      end


c-----------------------------------------------------------------------
