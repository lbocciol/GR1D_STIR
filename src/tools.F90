!-*-f90-*-
! various aux routines
 subroutine findnaninf(vector,length,nan,inf)

   implicit none
   real*8 vector(length)
   integer length
   integer i

   logical nan,inf

   nan = .false.
   inf = .false.


   do i=1,length
      if(vector(i).ne.vector(i)) then
         write(*,*) "NaN at index ",i
         nan = .true.
      endif
      
   enddo
   
 end subroutine findnaninf

 subroutine ApplyEOS_limits

  use GR1D_module
  use nulibtable
  implicit none

  integer :: i

  if (.not. limit_EoS_table) return

   do i=1, n1

      if ( i .le. M1_imaxradii+ghosts1-1 ) then
         if ( log10(temp(i)) .le. nulibtable_logtemp_min ) temp(i) = 10**nulibtable_logtemp_min
         if ( log10(temp(i)) .ge. nulibtable_logtemp_max ) temp(i) = 10**nulibtable_logtemp_max

         if ( ye(i) .le. nulibtable_ye_min ) ye(i) = nulibtable_ye_min
         if ( ye(i) .ge. nulibtable_ye_max ) ye(i) = nulibtable_ye_max

      else

         if ( rho(i)/rho_gf .le. EoS_rhomin ) rho(i) = EoS_rhomin*rho_gf
         if ( rho(i)/rho_gf .ge. EoS_rhomax ) rho(i) = EoS_rhomax*rho_gf

         if ( temp(i) .le. EoS_tempmin ) temp(i) = EoS_tempmin
         if ( temp(i) .ge. EoS_tempmax ) temp(i) = EoS_tempmax

         if ( ye(i) .le. EoS_yemin ) ye(i) = EoS_yemin
         if ( ye(i) .ge. EoS_yemax ) ye(i) = EoS_yemax
      endif
    enddo

 end subroutine ApplyEOS_limits


