c
c     |****************************************************************|
c     |           ---  EXACT VELOCITY AT THE INTERFACE  ---            |
c     |                          Sept 25,2010                          |
c     |****************************************************************|
c     |    Out variables:                                              |
c     |                       - E_Uxsys                                |
c     |                       - E_Vxsys                                |
c     |                                                                |
c     |    Programs:                                                   |
c     |                       - ExactSol                               |
c     |                                                                |
c     |    Remark: We write the solution in the data files:            |
c     |                       - E_Uxsys.dat                            |
c     |                       - E_Vxsys.dat                            |
c     |                                                                |
c     |****************************************************************|
c
      subroutine ExactUvelocity(E_Uxsys,E_Vxsys)
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'
      real*8 E_Uxsys(0:ns,ms), E_Vxsys(0:ns,ms)
      real*8 foo    

      DO l=1,ms
        do m=0,ns-1
          call ExactSol(xs(m,l),ys(m,l),E_Uxsys(m,l),E_Vxsys(m,l),foo)
        enddo
        E_Uxsys(ns,l)= E_Uxsys(0,l)
        E_Vxsys(ns,l)= E_Vxsys(0,l)
      ENDDO

      !----------------------------------------------------------------!
      !      Writing Exact interface velocity  at (xs,ys)              !
      !----------------------------------------------------------------!

400   format(1x,2000e16.6)

      open(unit=69,file='DAT/Uinterface/E_Uxsys.dat',status='unknown')
      open(unit=79,file='DAT/Uinterface/E_Vxsys.dat',status='unknown') 

      do m=0,ns-1
        write(69,400)(E_Uxsys(m,l),l=1,ms)
        write(79,400)(E_Vxsys(m,l),l=1,ms)
      enddo
       
      close(69)
      close(79)
            
      !----------------------------------------------------------------!
      return
      end

