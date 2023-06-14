c-----------------------------------------------------------------------
c
      subroutine singular_force
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      !==========================================!
      !    Modification for 2-fluids pressure    !
      !------------------------------------------!
      ! In this case, we will set:               !
      ! fy = [p]                                 !
      !------------------------------------------!

      DO l=1,ms
        do m=0,ns
          fx(m,l) = 0.0d0
          fy(m,l) = 0.0d0 ! initial guess 4 GMRES
        enddo
      ENDDO
     
      !==========================================!
     
      return
      end


c-----------------------------------------------------------------------
