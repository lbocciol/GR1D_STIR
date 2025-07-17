!-*-f90-*-
program GR1D	
  
  use OMP_LIB
  use GR1D_module
  implicit none

  real*8 :: tstart, tfin, t1, t2
  
  timer_code = omp_get_wtime()
  !Welcome to GR1D
  write(*,*) "#################################################"
  write(*,*) "#################################################"
  write(*,*) "########### GR1D SPHERICAL HYDRO v2 #############"
  write(*,*) "######### Now with Neutrino Transport ###########"
  write(*,*) "################# Nov ??, 2014 ##################"
  write(*,*) "#################################################"

  ! Call problem setup and allocate/initialize variables 
  call start
  write(*,*) "Done with initial data :-)"

  write(*,*) "Begin time integration loop:"
  IntegrationLoop: do 

     call SetTimeStep

     call handle_output

!!   Integrate
     t1 = omp_get_wtime()
     call Step(dt)
     t2 = omp_get_wtime()
     timer_step = timer_step + (t2 - t1)

     call postStep_analysis
     call flush(6)

  enddo IntegrationLoop
      
  write(*,*) "Shutting down!"
  write(*,*) " "

  call PrintTimers()

end program GR1D
