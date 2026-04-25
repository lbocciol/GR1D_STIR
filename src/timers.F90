!-*-f90-*-
module timers

  ! Add some useful timers
  real*8 :: timer_hydro   = 0.0d0
  real*8 :: timer_eosf    = 0.0d0
  real*8 :: timer_eoss    = 0.0d0
  real*8 :: timer_M1      = 0.0d0
  real*8 :: timer_fdhlle  = 0.0d0
  real*8 :: timer_rec     = 0.0d0
  real*8 :: timer_c2GR    = 0.0d0
  real*8 :: timer_c2p     = 0.0d0
  real*8 :: timer_M1_exp  = 0.0d0
  real*8 :: timer_M1_imp  = 0.0d0
  real*8 :: timer_step    = 0.0d0
  real*8 :: timer_M1_clo  = 0.0d0
  real*8 :: timer_M1_rec  = 0.0d0
  real*8 :: timer_M1_eas  = 0.0d0
  real*8 :: timer_M1cons  = 0.0d0
  real*8 :: timer_code    = 0.0d0
  real*8 :: t_start       = 0.0d0

contains

  subroutine start_timers
    implicit none

    timer_hydro   = 0.0d0
    timer_eosf    = 0.0d0
    timer_eoss    = 0.0d0
    timer_M1      = 0.0d0
    timer_fdhlle  = 0.0d0
    timer_rec     = 0.0d0
    timer_c2GR    = 0.0d0
    timer_c2p     = 0.0d0
    timer_M1_exp  = 0.0d0
    timer_M1_imp  = 0.0d0
    timer_step    = 0.0d0
    timer_M1_clo  = 0.0d0
    timer_M1_rec  = 0.0d0
    timer_M1_eas  = 0.0d0
    timer_M1cons  = 0.0d0
    timer_code    = 0.0d0

  end subroutine start_timers

#ifdef HAVE_OMP
  subroutine GetThisTime(t)

    use omp_lib
    implicit none
    real*8, intent(inout) :: t

    t = omp_get_wtime()
  
  end subroutine GetThisTime
#else
  subroutine GetThisTime(t)

    implicit none
    real*8, intent(inout) :: t

    call cpu_time(t)
  
  end subroutine GetThisTime
#endif

end module timers
