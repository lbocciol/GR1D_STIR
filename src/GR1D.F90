!-*-f90-*-
program GR1D	
    	  
  use GR1D_module
  use timers
  implicit none

  real(8) :: t1, t2

  CALL GetThisTime(timer_code)

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
