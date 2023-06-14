c-----------------------------------------------------------------------
c
      subroutine jc_pressure
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ic,jc,ic1,jc1
      real*8 signx,signy,r3,hp,hm,sx,sy
      real*8 duyp,duym,dvyp,dvym,duxp,duxm,dvxp,dvxm

      DO l=1,ms

      do i=0,ncxc(l)
        signx=falfaxc(2,i,l)
        signy=falfaxc(3,i,l)

        pjcxc(1,i,l)=signx*pjcxc(1,i,l)
        pjcxc(3,i,l)=signx*pjcxc(3,i,l)
        pjcxc(5,i,l)=signx*pjcxc(5,i,l)

        pjcxc(2,i,l)=signy*pjcxc(2,i,l)
        pjcxc(4,i,l)=signy*pjcxc(4,i,l)
        pjcxc(6,i,l)=signy*pjcxc(6,i,l)
      enddo

      do j=0,ncyc(l)
        signx=falfayc(2,j,l)
        signy=falfayc(3,j,l)

        pjcyc(1,j,l)=signx*pjcyc(1,j,l)
        pjcyc(3,j,l)=signx*pjcyc(3,j,l)
        pjcyc(5,j,l)=signx*pjcyc(5,j,l)

        pjcyc(2,j,l)=signy*pjcyc(2,j,l)
        pjcyc(4,j,l)=signy*pjcyc(4,j,l)
        pjcyc(6,j,l)=signy*pjcyc(6,j,l)
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
