!-*-f90-*-
program GR1D
    	 
#ifdef HAVE_OMP
  use omp_lib
#endif
  use GR1D_module
  use timers
  implicit none

  real(8) :: t1, t2
  integer :: i

  CALL GetThisTime(t_start)
  timer_code = t_start

  !Welcome to GR1D
  write(*,*) "#################################################"
  write(*,*) "#################################################"
  write(*,*) "########### GR1D SPHERICAL HYDRO v2 #############"
  write(*,*) "######### Now with Neutrino Transport ###########"
  write(*,*) "################# Nov ??, 2014 ##################"
  write(*,*) "#################################################"

#ifdef HAVE_OMP
  write(*,*) "Running with", omp_get_max_threads(), "cores"
#else
write(*,*) "Running in serial mode, you sure you don't want OMP?"
#endif
  ! Call problem setup and allocate/initialize variables 
  call start
  write(*,*) "Done with initial data :-)"

  write(*,*) "Begin time integration loop:"
  IntegrationLoop: do 
    do i=ghosts1+1,n1-ghosts1
      if (temp(i) .lt. 1.0d-5) then
          write(*,*) "temp do loop 1: ", temp(i)
          stop "temperature too low before Step"
      endif
    enddo
     call SetTimeStep
    do i=ghosts1+1,n1-ghosts1
      if (temp(i) .lt. 1.0d-5) then
          write(*,*) "temp do loop 2: ", temp(i)
          stop "temperature too low before Step"
      endif
    enddo
     call handle_output

    do i=ghosts1+1,n1-ghosts1
      if (temp(i) .lt. 1.0d-5) then
          write(*,*) "temp do loop 3: ", temp(i)
          stop "temperature too low before Step"
      endif
    enddo
!!   Integrate
     CALL GetThisTime(t1)
     call Step(dt)
     CALL GetThisTime(t2)
     timer_step = timer_step + (t2 - t1)

     call postStep_analysis
     call flush(6)

  enddo IntegrationLoop
      
  write(*,*) "Shutting down!"
  write(*,*) " "

  call PrintTimers()

end program GR1D
