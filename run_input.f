c-----------------------------------------------------------------------
c
      subroutine run_input
      include 'parameter.inc'
      include 'surface.inc'

      open(unit=8,file='DAT/input.run',status='old')
      rewind 8
      read(8,*)
      read(8,*)nstep,icfl,iread,iwrite,itout,iplot,ianimation
      read(8,*)
      read(8,*)
      read(8,*)cflc,cflv,dtcfl,dt0,tout
      close(8)

      open(unit=88,file='DAT/objects.run',status='old')
      rewind(88)
      read(88,*)
      do l=1,ms
        read(88,*)jectob(l),a(l),b(l),xsc(l),ysc(l),theta(l)
      enddo
      close(88)

      open(unit=16,file='DAT/st.dat',status='unknown')
      open(unit=26,file='DAT/pt.dat',status='unknown')
      open(unit=36,file='DAT/uct.dat',status='unknown')
      open(unit=46,file='DAT/vct.dat',status='unknown')
      open(unit=56,file='DAT/wot.dat',status='unknown')
      open(unit=66,file='DAT/cdcl.dat',status='unknown')
      open(unit=76,file='DAT/xsys.dat',status='unknown')

      return
      end


c-----------------------------------------------------------------------

