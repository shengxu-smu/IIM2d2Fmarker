c-----------------------------------------------------------------------
c
      subroutine surface_plot
      include 'parameter.inc'
      include 'surface.inc'
      character*1 fnext
      character*16 fname

      DO l=1,ms

      write(unit=fnext,fmt='(i1)') l
      fname='DAT/surface'//fnext//'.dat'
      open(unit=69,file=fname,status='unknown')
      do m=0,ns
        write(69,300)alfa(m),xs(m,l),ys(m,l),
     .               us(m,l),vs(m,l),us0(m,l),vs0(m,l),
     .               fx(m,l),dfy(m,l),
     .               (pjc(n,m,l),n=1,6)
      enddo
      close(69)

      ENDDO

      l=1
      open(unit=68,file='DAT/falfa.dat',status='unknown')
      do j=0,ncyc(l)
        write(68,300)falfayc(0,j,l),(falfayc(2,j,l)*ujcyc(n,j,l),n=1,4)
     .              ,(vjcyc(n,j,l),falfayc(3,j,l)*pjcyc(n,j,l),n=1,4)
      enddo
      close(68)

300   format(1x,30e16.6)

      return
      end


c-----------------------------------------------------------------------
