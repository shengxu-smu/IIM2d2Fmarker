c-----------------------------------------------------------------------
c
      subroutine rk
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer iflag,krk,modify
      real*8 fac, MethodSolver

      MethodSolver=1
      modify=1
     
      !--------------------------------------!
      !  Calculating the Step function:      !
      !  We need to do it in each time step  !
      !  because the shape of the interface  !
      !  will change.                        !
      !--------------------------------------!     
      
c  krk=1:
      krk=1
      fac=0.5d0

      if(isingular.eq.1) then
        call surface_property
        call euler_link
        call StepFun     ! <-----------------!
        call Jumps
        call singular_force
        call jc_rhs
        call jc_firstsecond
        call mac_distribute
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call correction_velocity
      else
        call field_interpolate
        call uv_strain
      endif
   
      call old_save

      !-------------------------------!
      !          Pressure             ! 
      !-------------------------------!

      IF (MethodSolver==1) THEN       
        if(isingular.eq.1) then   
          call jc_pressure
          call correction_pressure
        endif
        call pressure(krk,fac)

        if(modify.eq.1) then 
          call surface_pressure
          call gmres_pressure(fac)
        endif
      ELSEIF (MethodSolver==2) THEN
          call IterPressureSolver
      ENDIF
     
      !-------------------------------!

      call rhs(krk)
 1    call pgrad(krk)

      if(move.eq.1) then
        call surface_moveMARKER(krk,fac)
        !call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,1)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,1)
        enddo
      enddo 
      call ubc(krk,fac)
      call vbc(krk,fac)
 
c  krk=2:
      krk=2
      fac=0.5d0
      if(isingular.eq.1) then
        call surface_property
        call euler_link
        call StepFun
        call Jumps
        call singular_force
        call jc_rhs
        call jc_firstsecond
        call mac_distribute
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call correction_velocity
      else
        call field_interpolate
        call uv_strain
      endif

      !-------------------------------!
      !          Pressure             ! 
      !-------------------------------!

      IF (MethodSolver==1) THEN       
        if(isingular.eq.1) then   
          call jc_pressure
          call correction_pressure
        endif
        call pressure(krk,fac)

        if(modify.eq.1) then 
          call surface_pressure
          call gmres_pressure(fac)
        endif
      ELSEIF (MethodSolver==2) THEN
          call IterPressureSolver
      ENDIF
     
      !-------------------------------!

      call rhs(krk)
2     call pgrad(krk)
      
      if(move.eq.1) then
        call surface_moveMARKER(krk,fac)
        !call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,2)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,2)
        enddo
      enddo
      call ubc(krk,fac)
      call vbc(krk,fac)

c  krk=3:
      krk=3
      fac=1.0d0

      if(isingular.eq.1) then
        call surface_property
        call euler_link
        call StepFun
        call Jumps
        call singular_force
        call jc_rhs
        call jc_firstsecond
        call mac_distribute
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call correction_velocity
      else
        call field_interpolate
        call uv_strain
      endif

      !-------------------------------!
      !          Pressure             ! 
      !-------------------------------!

      IF (MethodSolver==1) THEN       
        if(isingular.eq.1) then   
          call jc_pressure
          call correction_pressure
        endif
        call pressure(krk,fac)

        if(modify.eq.1) then 
          call surface_pressure
          call gmres_pressure(fac)
        endif
      ELSEIF (MethodSolver==2) THEN
          call IterPressureSolver
      ENDIF
     
      !-------------------------------!

      call rhs(krk)
 3    call pgrad(krk)

      if(move.eq.1) then
        call surface_moveMARKER(krk,fac)
        !call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,3)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,3)
        enddo
      enddo
      call ubc(krk,fac)
      call vbc(krk,fac)

c  krk=4:
      krk=4
      fac=1.0d0

      if(isingular.eq.1) then
        call surface_property
        call euler_link
        call StepFun
        call Jumps
        call singular_force
        call jc_rhs
        call jc_firstsecond
        call mac_distribute
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call correction_velocity
      else
        call field_interpolate
        call uv_strain
      endif

      !-------------------------------!
      !          Pressure             ! 
      !-------------------------------!

      IF (MethodSolver==1) THEN       
        if(isingular.eq.1) then   
          call jc_pressure
          call correction_pressure
        endif
        call pressure(krk,fac)

        if(modify.eq.1) then 
          call surface_pressure
          call gmres_pressure(fac)
        endif
      ELSEIF (MethodSolver==2) THEN
          call IterPressureSolver
      ENDIF
     
      !-------------------------------!

      call rhs(krk)
 4    call pgrad(krk)

      if(move.eq.1) then
        call surface_moveMARKER(krk,fac)
        !call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+(dt/6.0d0)*(urk(i,j,1)+
     .                   2.0d0*(urk(i,j,2)+urk(i,j,3))+urk(i,j,4))
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+(dt/6.0d0)*(vrk(i,j,1)+
     .                   2.0d0*(vrk(i,j,2)+vrk(i,j,3))+vrk(i,j,4))
        enddo
      enddo
      call ubc(krk,fac)
      call vbc(krk,fac)

      return
      end


c-----------------------------------------------------------------------
