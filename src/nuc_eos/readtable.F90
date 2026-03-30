!SVN:EOSdriver Revision #7
subroutine readtable(eos_filename)
! This routine reads the table and initializes
! all variables in the module. 

  use eosmodule
  use hdf5 

  implicit none

  character(*) eos_filename

  character(len=100) message

! HDF5 vars
  integer(HID_T) file_id,dset_id,dspace_id
  integer(HSIZE_T) dims1(1), dims3(3)
  integer error,rank,accerr
  integer i,j,k

  real*8 amu_cgs_andi
  real*8 buffer1,buffer2,buffer3,buffer4
  accerr=0

  write(*,*) "Reading Nuclear EOS Table"

  call h5open_f(error)

  call h5fopen_f (trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)

  write(6,*) trim(adjustl(eos_filename))

! read scalars
  dims1(1)=1
  call h5dopen_f(file_id, "pointsrho", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
  call h5dclose_f(dset_id,error)

  if(error.ne.0) then
     stop "Could not read EOS table file"
  endif

  dims1(1)=1
  call h5dopen_f(file_id, "pointstemp", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
  call h5dclose_f(dset_id,error)

  if(error.ne.0) then
     stop "Could not read EOS table file"
  endif

  dims1(1)=1
  call h5dopen_f(file_id, "pointsye", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nye, dims1, error)
  call h5dclose_f(dset_id,error)

  if(error.ne.0) then
     stop "Could not read EOS table file"
  endif

  write(message,"(a25,i5,i5,i5)") "We have nrho ntemp nye: ", nrho,ntemp,nye
  write(*,*) message

  allocate(alltables(nrho,ntemp,nye,nvars))

  ! index variable mapping:
  !  1 -> logpress
  !  2 -> logenergy
  !  3 -> entropy
  !  4 -> munu
  !  5 -> cs2
  !  6 -> dedT
  !  7 -> dpdrhoe
  !  8 -> dpderho
  !  9 -> muhat
  ! 10 -> mu_e
  ! 11 -> mu_p
  ! 12 -> mu_n
  ! 13 -> xa
  ! 14 -> xh
  ! 15 -> xn
  ! 16 -> xp
  ! 17 -> abar
  ! 18 -> zbar


  dims3(1)=nrho
  dims3(2)=ntemp
  dims3(3)=nye
  call h5dopen_f(file_id, "logpress", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,1), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  call h5dopen_f(file_id, "logenergy", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,2), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  call h5dopen_f(file_id, "entropy", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,3), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  call h5dopen_f(file_id, "munu", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,4), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  call h5dopen_f(file_id, "cs2", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,5), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  call h5dopen_f(file_id, "dedt", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,6), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  call h5dopen_f(file_id, "dpdrhoe", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,7), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  call h5dopen_f(file_id, "dpderho", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,8), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

! chemical potentials
  call h5dopen_f(file_id, "muhat", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,9), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  call h5dopen_f(file_id, "mu_e", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,10), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  call h5dopen_f(file_id, "mu_p", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,11), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  call h5dopen_f(file_id, "mu_n", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,12), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

! compositions
  call h5dopen_f(file_id, "Xa", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,13), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  call h5dopen_f(file_id, "Xh", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,14), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  call h5dopen_f(file_id, "Xn", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,15), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  call h5dopen_f(file_id, "Xp", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,16), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error


! average nucleus
  call h5dopen_f(file_id, "Abar", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,17), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  call h5dopen_f(file_id, "Zbar", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,18), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

! Gamma
  call h5dopen_f(file_id, "gamma", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,19), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  allocate(logrho(nrho))
  dims1(1)=nrho
  call h5dopen_f(file_id, "logrho", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  allocate(logtemp(ntemp))
  dims1(1)=ntemp
  call h5dopen_f(file_id, "logtemp", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  allocate(ye(nye))
  dims1(1)=nye
  call h5dopen_f(file_id, "ye", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error


  call h5dopen_f(file_id, "energy_shift", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, energy_shift, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  if(accerr.ne.0) then
    stop "Problem reading EOS table file"
  endif


  call h5fclose_f (file_id,error)

  call h5close_f (error)

  ! set min-max values:

  eos_rhomin = 10.0d0**logrho(1)
  eos_rhomax = 10.0d0**logrho(nrho)

  eos_yemin = ye(1)
  eos_yemax = ye(nye)

  eos_tempmin = 10.0d0**logtemp(1)
  eos_tempmax = 10.0d0**logtemp(ntemp)

  write(6,*) "Done reading eos tables"


end subroutine readtable

subroutine get_extrema_EoS_table( rhomin, rhomax, &
        yemin, yemax, tmin, tmax )

  use eosmodule
  implicit none

  real*8, intent(out) :: rhomin, rhomax
  real*8, intent(out) :: yemin, yemax
  real*8, intent(out) :: tmin, tmax

  rhomin = eos_rhomin
  rhomax = eos_rhomax

  yemin = eos_yemin
  yemax = eos_yemax

  tmin = eos_tempmin
  tmax = eos_tempmax

end subroutine get_extrema_EoS_table

! Extra things to calculate BV frequency

subroutine Calculate_additional_derivs

  use eosmodule

  integer :: ipress = 1, ientropy = 3, imunu = 4, ics2 = 5
  integer :: i,j,k,iE

  real*8 :: x1, x2, f1, f2

  allocate(ynu(nrho,ntemp,nye))
  allocate(dpdye(nrho,ntemp,nye))
  allocate(dyLdye(nrho,ntemp,nye))
  allocate(dsdye(nrho,ntemp,nye))
  allocate(dyLdT(nrho,ntemp,nye))
  allocate(dpdT(nrho,ntemp,nye))
  allocate(dsdT(nrho,ntemp,nye))
  allocate(dyLdrho(nrho,ntemp,nye))
  allocate(dsdrho(nrho,ntemp,nye))
  allocate(dpdrho(nrho,ntemp,nye))


  allocate(dpdyL_s_rho(nrho,ntemp,nye))
  allocate(dpdye_s_rho(nrho,ntemp,nye))
  allocate(dpds_rho_yL(nrho,ntemp,nye))
  allocate(dpds_rho_ye(nrho,ntemp,nye))
  allocate(dpdrho_s_yL(nrho,ntemp,nye))
  allocate(dpdrho_s_ye(nrho,ntemp,nye))

  ynu = 0.0d0
  dpdye = 0.0d0
  dyLdye = 0.0d0
  dsdye = 0.0d0
  dyLdT = 0.0d0
  dpdT = 0.0d0
  dsdT = 0.0d0
  dyLdrho = 0.0d0
  dsdrho = 0.0d0
  dpdrho = 0.0d0
  dpdyL_s_rho = 0.0d0
  dpds_rho_yL = 0.0d0
  dpdrho_s_yL = 0.0d0

  do k=1,nye
    do i=1,nrho
      do j=1,ntemp

        do iE=1,number_groups_eos
          fnu = 1.0d0 / (1.0d0 + exp((nu_energies_eos(iE) - &
              alltables(i,j,k,imunu))/10.0d0**logtemp(j)))
          fanu = 1.0d0 / (1.0d0 + exp((nu_energies_eos(iE) + &
              alltables(i,j,k,imunu))/10.0d0**logtemp(j)))
          ynu(i,j,k) = ynu(i,j,k) + 4.0d0*pi / 10.0d0**logrho(i) / (hbarc_mevcm*2.0d0*pi)**3 * &
              (fnu - fanu) * nu_energies_eos(iE)**2 * nu_ewidths_eos(iE) * amu_cgs
        enddo

      enddo
    enddo
  enddo

  ! Derivatives at constant temperature
  do j=1,ntemp
     do i=1,nrho
        ! dpdye
        do k=2,nye-1
           x1 = ye(k-1)
           f1 = alltables(i,j,k-1,ipress)
           x2 = ye(k+1)
           f2 = alltables(i,j,k+1,ipress)
           dpdye(i,j,k) = (f2-f1)/(x2-x1) * &
                 10.0d0**alltables(i,j,k,ipress)
        enddo

        ! boundaries: one-sided derivative
        x1 = ye(1)
        f1 = alltables(i,j,1,ipress)
        x2 = ye(2)
        f2 = alltables(i,j,2,ipress)
        dpdye(i,1,k) = (f2-f1)/(x2-x1) * &
            10.0d0**alltables(1,j,k,ipress)

        x1 = ye(nye-1)
        f1 = alltables(i,j,nye-1,ipress)
        x2 = ye(nye)
        f2 = alltables(i,j,nye,ipress)
        dpdye(i,j,nye) = (f2-f1)/(x2-x1) * &
             10.0d0**alltables(i,1,k,ipress)

        ! dyLdye
        do k=2,nye-1
           x1 = ye(k-1)
           f1 = ye(k-1) + ynu(i,j,k-1)
           x2 = ye(k+1)
           f2 = ye(k+1) + ynu(i,j,k+1)
           dyLdye(i,j,k) = (f2-f1)/(x2-x1)
        enddo

        ! boundaries: one-sided derivative
        x1 = ye(1)
        f1 = ye(1) + ynu(i,j,1)
        x2 = ye(2)
        f2 = ye(2) + ynu(i,j,2)
        dyLdye(i,1,k) = (f2-f1)/(x2-x1)

        x1 = ye(nye-1)
        f1 = ye(nye-1) + ynu(i,j,nye-1)
        x2 = ye(nye)
        f2 = ye(nye) + ynu(i,j,nye)
        dyLdye(i,j,nye) = (f2-f1)/(x2-x1)

        ! dsdye
        do k=2,nye-1
           x1 = ye(k-1)
           f1 = log10(alltables(i,j,k-1,ientropy))
           x2 = ye(k+1)
           f2 = log10(alltables(i,j,k+1,ientropy))
           dsdye(i,j,k) = (f2-f1)/(x2-x1) * alltables(i,j,k,ientropy)
        enddo

        ! boundaries: one-sided derivative
        x1 = ye(1)
        f1 = log10(alltables(i,j,1,ientropy))
        x2 = ye(2)
        f2 = log10(alltables(i,j,2,ientropy))
        dsdye(i,1,k) = (f2-f1)/(x2-x1) * alltables(i,j,1,ientropy)

        x1 = ye(nye-1)
        f1 = log10(alltables(i,j,nye-1,ientropy))
        x2 = ye(nye)
        f2 = log10(alltables(i,j,nye,ientropy))
        dsdye(i,j,nye) = (f2-f1)/(x2-x1) * alltables(i,j,nye,ientropy)

     enddo
  enddo

  ! Derivatives at constant Ye and rho
  do k=1,nye
     do i=1,nrho
        ! dyLdT
        do j=2,ntemp-1
           x1 = logtemp(j-1)
           f1 = ye(k) + ynu(i,j-1,k)
           x2 = logtemp(j+1)
           f2 = ye(k) + ynu(i,j+1,k)
           dyLdT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)
        enddo

        ! boundaries: one-sided derivative
        x1 = logtemp(1)
        f1 = ye(k) + ynu(i,1,k)
        x2 = logtemp(2)
        f2 = ye(k) + ynu(i,2,k)
        dyLdT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1)

        x1 = logtemp(ntemp-1)
        f1 = ye(k) + ynu(i,ntemp-1,k)
        x2 = logtemp(ntemp)
        f2 = ye(k) + ynu(i,ntemp,k)
        dyLdT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(ntemp)

        ! dpdT
        do j=2,ntemp-1
           x1 = logtemp(j-1)
           f1 = alltables(i,j-1,k,ipress)
           x2 = logtemp(j+1)
           f2 = alltables(i,j+1,k,ipress)
           dpdT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j) * &
                10.0d0**alltables(i,j,k,ipress)
        enddo

        ! boundaries: one-sided derivative
        x1 = logtemp(1)
        f1 = alltables(i,1,k,ipress)
        x2 = logtemp(2)
        f2 = alltables(i,2,k,ipress)
        dpdT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
                10.0d0**alltables(i,1,k,ipress)

        x1 = logtemp(ntemp-1)
        f1 = alltables(i,ntemp-1,k,ipress)
        x2 = logtemp(ntemp)
        f2 = alltables(i,ntemp,k,ipress)
        dpdT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(ntemp) * &
                10.0d0**alltables(i,ntemp,k,ipress)

        ! dsdT
        do j=2,ntemp-1
           x1 = logtemp(j-1)
           f1 = log10(alltables(i,j-1,k,ientropy))
           x2 = logtemp(j+1)
           f2 = log10(alltables(i,j+1,k,ientropy))
           dsdT(i,j,k) = (f2-f1)/(x2-x1) &
               / 10.0d0**logtemp(j) * alltables(i,j,k,ientropy)
        enddo

        ! boundaries: one-sided derivative
        x1 = logtemp(1)
        f1 = log10(alltables(i,1,k,ientropy))
        x2 = logtemp(2)
        f2 = log10(alltables(i,2,k,ientropy))
        dsdT(i,1,k) = (f2-f1)/(x2-x1) &
            / 10.0d0**logtemp(1) * alltables(i,1,k,ientropy)

        x1 = logtemp(ntemp-1)
        f1 = log10(alltables(i,ntemp-1,k,ientropy))
        x2 = logtemp(ntemp)
        f2 = log10(alltables(i,ntemp,k,ientropy))
        dsdT(i,ntemp,k) = (f2-f1)/(x2-x1) &
            / 10.0d0**logtemp(ntemp) * alltables(i,ntemp,k,ientropy)

     enddo
  enddo

  ! Derivatives at constant Ye and rho
  do k=1,nye
     do j=1,ntemp
        ! dyLdrho
       do i=2,nrho-1
           x1 = logrho(i-1)
           f1 = ye(k) + ynu(i-1,j,k)
           x2 = logrho(i+1)
           f2 = ye(k) + ynu(i+1,j,k)
           dyLdrho(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(i)
        enddo

        ! boundaries: one-sided derivative
        x1 = logrho(1)
        f1 = ye(k) + ynu(1,j,k)
        x2 = logrho(2)
        f2 = ye(k) + ynu(2,j,k)
        dyLdrho(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1)

        x1 = logrho(nrho-1)
        f1 = ye(k) + ynu(nrho-1,j,k)
        x2 = logrho(nrho)
        f2 = ye(k) + ynu(nrho,j,k)
        dyLdrho(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho)


       ! dsdrho
       do i=2,nrho-1
           x1 = logrho(i-1)
           f1 = log10(alltables(i-1,j,k,ientropy))
           x2 = logrho(i+1)
           f2 = log10(alltables(i+1,j,k,ientropy))
           dsdrho(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(i) &
                * alltables(i,j,k,ientropy)
        enddo

        ! boundaries: one-sided derivative
        x1 = logrho(1)
        f1 = log10(alltables(1,j,k,ientropy))
        x2 = logrho(2)
        f2 = log10(alltables(2,j,k,ientropy))
        dsdrho(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
               * alltables(1,j,k,ientropy)

        x1 = logrho(nrho-1)
        f1 = log10(alltables(nrho-1,j,k,ientropy))
        x2 = logrho(nrho)
        f2 = log10(alltables(nrho,j,k,ientropy))
        dsdrho(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
                * alltables(nrho,j,k,ientropy)


       ! dpdrho
       do i=2,nrho-1
           x1 = logrho(i-1)
           f1 = alltables(i-1,j,k,ipress)
           x2 = logrho(i+1)
           f2 = alltables(i+1,j,k,ipress)
           dpdrho(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(i) &
                * alltables(i,j,k,ipress)
        enddo

        ! boundaries: one-sided derivative
        x1 = logrho(1)
        f1 = alltables(1,j,k,ipress)
        x2 = logrho(2)
        f2 = alltables(2,j,k,ipress)
        dpdrho(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
               * alltables(1,j,k,ipress)

        x1 = logrho(nrho-1)
        f1 = alltables(nrho-1,j,k,ipress)
        x2 = logrho(nrho)
        f2 = alltables(nrho,j,k,ipress)
        dpdrho(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
                * 10.0d0**alltables(nrho,j,k,ipress)

     enddo
   enddo

   ! -----------------------------------------
   dpdyL_s_rho = dpdye / dyLdye + dsdye / dyLdye * dpdye / dyLdye * &
        ( ( dyLdT - dyLdye / dpdye * dpdT) / &
        (dsdT - dsdye / dyLdye * dyLdT) )

   dpdye_s_rho = dpdye - dsdye * dpdT / dsdt

   ! -----------------------------------------
   dpds_rho_yL = - dpdye / dyLdye * &
       ( dyLdT - dyLdye / dpdye * dpdT ) &
     / ( dsdT - dsdye / dyLdye * dyLdT )

   dpds_rho_ye = dpdt / dsdt

   ! -----------------------------------------
   dpdrho_s_yL = alltables(:,:,:,ics2) - &
       ( dpdye - dpdT/dsdT * dsdye ) / &
       ( dyLdye - dyLdT/dsdT * dsdye ) * &
       (dyLdrho - dyLdT/dsdT * dsdrho )

   dpdrho_s_ye = alltables(:,:,:,ics2)

! NOW CALCULATE LOGARITHMIC DERIVATIVES LIKE IN PASCAL 2022
  allocate(dlnpdlns_rho_ye(nrho,ntemp,nye))
  allocate(dlnpdlnye_s_rho(nrho,ntemp,nye))
  allocate(dlnsdlnt(nrho,ntemp,nye))
  allocate(dlnpdlnt(nrho,ntemp,nye))
  allocate(dlnPdlnrho_s_ye(nrho,ntemp,nye))
  allocate(dlnPdlnrho(nrho,ntemp,nye))
  allocate(dlnsdlnrho(nrho,ntemp,nye))
  allocate(dlnsdlnye(nrho,ntemp,nye))
  allocate(dlnpdlnye(nrho,ntemp,nye))

  dlnpdlns_rho_ye(:,:,:) = 0.0d0
  dlnpdlnye_s_rho(:,:,:) = 0.0d0
  dlnsdlnt(:,:,:) = 0.0d0
  dlnpdlnt(:,:,:) = 0.0d0
  dlnPdlnrho_s_ye(:,:,:) = 0.0d0
  dlnPdlnrho(:,:,:) = 0.0d0
  dlnsdlnrho(:,:,:) = 0.0d0
  dlnsdlnye(:,:,:) = 0.0d0
  dlnpdlnye(:,:,:) = 0.0d0
  ! Derivatives at constant temperature
  do j=1,ntemp
     do i=1,nrho
        ! dpdye
        do k=2,nye-1
           x1 = log(ye(k-1))
           f1 = log(10.0**alltables(i,j,k-1,ipress))
           x2 = log(ye(k+1))
           f2 = log(10.0**alltables(i,j,k+1,ipress))
           dlnpdlnye(i,j,k) = (f2-f1)/(x2-x1)
        enddo

        ! boundaries: one-sided derivative
        x1 = log(ye(1))
        f1 = log(10.0**alltables(i,j,1,ipress))
        x2 = log(ye(2))
        f2 = log(10.0**alltables(i,j,2,ipress))
        dlnpdlnye(i,1,k) = (f2-f1)/(x2-x1)

        x1 = log(ye(nye-1))
        f1 = log(10.0**alltables(i,j,nye-1,ipress))
        x2 = log(ye(nye))
        f2 = log(10.0**alltables(i,j,nye,ipress))
        dlnpdlnye(i,j,nye) = (f2-f1)/(x2-x1)

        ! dsdye
        do k=2,nye-1
           x1 = log(ye(k-1))
           f1 = log(alltables(i,j,k-1,ientropy))
           x2 = log(ye(k+1))
           f2 = log(alltables(i,j,k+1,ientropy))
           dlnsdlnye(i,j,k) = (f2-f1)/(x2-x1)
        enddo

        ! boundaries: one-sided derivative
        x1 = log(ye(1))
        f1 = log(alltables(i,j,1,ientropy))
        x2 = log(ye(2))
        f2 = log(alltables(i,j,2,ientropy))
        dlnsdlnye(i,1,k) = (f2-f1)/(x2-x1)

        x1 = log(ye(nye-1))
        f1 = log(alltables(i,j,nye-1,ientropy))
        x2 = log(ye(nye))
        f2 = log(alltables(i,j,nye,ientropy))
        dlnsdlnye(i,j,nye) = (f2-f1)/(x2-x1)

     enddo
  enddo


  ! Derivatives at constant Ye and T
  do k=1,nye
     do j=1,ntemp
        ! dPdrho
        do i=2,nrho-1
           x1 = log(10.0**logrho(i-1))
           f1 = log(10.0**alltables(i-1,j,k,ipress))
           x2 = log(10.0**logrho(i+1))
           f2 = log(10.0**alltables(i+1,j,k,ipress))
           dlnPdlnrho(i,j,k) = (f2-f1)/(x2-x1)
        enddo

        ! boundaries: one-sided derivative
        x1 = log(10.0**logrho(1))
        f1 = log(10.0**alltables(1,j,k,ipress))
        x2 = log(10.0**logrho(2))
        f2 = log(10.0**alltables(2,j,k,ipress))
        dlnPdlnrho(1,j,k) = (f2-f1)/(x2-x1)

        x1 = log(10.0**logrho(nrho-1))
        f1 = log(10.0**alltables(nrho-1,j,k,ipress))
        x2 = log(10.0**logrho(nrho))
        f2 = log(10.0**alltables(nrho,j,k,ipress))
        dlnPdlnrho(nrho,j,k) = (f2-f1)/(x2-x1)

        ! dsdrho
        do i=2,nrho-1
           x1 = log(10.0**logrho(i-1))
           f1 = log(alltables(i-1,j,k,ientropy))
           x2 = log(10.0**logrho(i+1))
           f2 = log(alltables(i+1,j,k,ientropy))
           dlnsdlnrho(i,j,k) = (f2-f1)/(x2-x1)
        enddo

        ! boundaries: one-sided derivative
        x1 = log(10.0**logrho(1))
        f1 = log(alltables(1,j,k,ientropy))
        x2 = log(10.0**logrho(2))
        f2 = log(alltables(2,j,k,ientropy))
        dlnsdlnrho(1,j,k) = (f2-f1)/(x2-x1)

        x1 = log(10.0**logrho(nrho-1))
        f1 = log(alltables(nrho-1,j,k,ientropy))
        x2 = log(10.0**logrho(nrho))
        f2 = log(alltables(nrho,j,k,ientropy))
        dlnsdlnrho(nrho,j,k) = (f2-f1)/(x2-x1)

     enddo
  enddo

  ! Derivatives at constant Ye and rho
  do k=1,nye
     do i=1,nrho
        ! dpdT
        do j=2,ntemp-1
           x1 = log(10.0**logtemp(j-1))
           f1 = log(10.0**alltables(i,j-1,k,ipress))
           x2 = log(10.0**logtemp(j+1))
           f2 = log(10.0**alltables(i,j+1,k,ipress))
           dlnpdlnt(i,j,k) = (f2-f1)/(x2-x1)
        enddo

        ! boundaries: one-sided derivative
        x1 = log(10.0**logtemp(1))
        f1 = log(10.0**alltables(i,1,k,ipress))
        x2 = log(10.0**logtemp(2))
        f2 = log(10.0**alltables(i,2,k,ipress))
        dlnpdlnt(i,1,k) = (f2-f1)/(x2-x1)

        x1 = log(10.0**logtemp(ntemp-1))
        f1 = log(10.0**alltables(i,ntemp-1,k,ipress))
        x2 = log(10.0**logtemp(ntemp))
        f2 = log(10.0**alltables(i,ntemp,k,ipress))
        dlnpdlnt(i,ntemp,k) = (f2-f1)/(x2-x1)

        ! dsdT
        do j=2,ntemp-1
           x1 = log(10.0**logtemp(j-1))
           f1 = log(alltables(i,j-1,k,ientropy))
           x2 = log(10.0**logtemp(j+1))
           f2 = log(alltables(i,j+1,k,ientropy))
           dlnsdlnt(i,j,k) = (f2-f1)/(x2-x1)
        enddo

        ! boundaries: one-sided derivative
        x1 = log(10.0**logtemp(1))
        f1 = log(alltables(i,1,k,ientropy))
        x2 = log(10.0**logtemp(2))
        f2 = log(alltables(i,2,k,ientropy))
        dlnsdlnt(i,1,k) = (f2-f1)/(x2-x1)

        x1 = log(10.0**logtemp(ntemp-1))
        f1 = log(alltables(i,ntemp-1,k,ientropy))
        x2 = log(10.0**logtemp(ntemp))
        f2 = log(alltables(i,ntemp,k,ientropy))
        dlnsdlnt(i,ntemp,k) = (f2-f1)/(x2-x1)

     enddo
  enddo

  dlnpdlnye_s_rho = dlnpdlnye - dlnsdlnye * dlnpdlnt / dlnsdlnt

  dlnpdlns_rho_ye = dlnpdlnt / dlnsdlnt

  dlnpdlnrho_s_ye = dlnpdlnrho - dlnsdlnrho * dlnpdlnt / dlnsdlnt

end subroutine Calculate_additional_derivs

subroutine Interpolate_additional_derivs( xrho, xtemp, xye, &
        this_dpds_rho_ye, this_dpdye_s_rho, this_dpdrho_s_ye )

  use eosmodule

  real*8 xrho, xtemp, xye !input
  real*8 this_dpds_rho_ye, this_dpdye_s_rho, this_dpdrho_s_ye !output
  real*8 ffx(1,1)
  real*8 lr, lt, y

  lr = log10( xrho )
  lt = log10( xtemp )
  y = xye

  !dpds_rho_ye
  call intp3d_many(lr,lt,y,ffx,1,dpds_rho_ye, &
        nrho,ntemp,nye,1,logrho,logtemp,ye)

  this_dpds_rho_ye = ffx(1,1)

  ! dpdye_s_rho
  call intp3d_many(lr,lt,y,ffx,1,dpdye_s_rho, &
        nrho,ntemp,nye,1,logrho,logtemp,ye)

  this_dpdye_s_rho = ffx(1,1)

  ! dpdrho_s_ye
  call intp3d_many(lr,lt,y,ffx,1,dpdrho_s_ye, &
        nrho,ntemp,nye,1,logrho,logtemp,ye)

  this_dpdrho_s_ye = ffx(1,1)

end subroutine Interpolate_additional_derivs

subroutine Interpolate_additional_lnderivs( xrho, xtemp, xye, &
        this_dlnpdlns_rho_ye, this_dlnpdlnye_s_rho, this_dlnpdlnrho_s_ye )

  use eosmodule

  real*8 xrho, xtemp, xye !input
  real*8 this_dlnpdlns_rho_ye, this_dlnpdlnye_s_rho, this_dlnpdlnrho_s_ye !output
  real*8 ffx(1,1)
  real*8 lr, lt, y

  lr = log10( xrho )
  lt = log10( xtemp )
  y = xye

  !dpds_rho_ye
  call intp3d_many(lr,lt,y,ffx,1,dlnpdlns_rho_ye, &
        nrho,ntemp,nye,1,logrho,logtemp,ye)

  this_dlnpdlns_rho_ye = ffx(1,1)

  ! dpdye_s_rho
  call intp3d_many(lr,lt,y,ffx,1,dlnpdlnye_s_rho, &
        nrho,ntemp,nye,1,logrho,logtemp,ye)

  this_dlnpdlnye_s_rho = ffx(1,1)

  ! dpdrho_s_ye
  call intp3d_many(lr,lt,y,ffx,1,dlnpdlnrho_s_ye, &
        nrho,ntemp,nye,1,logrho,logtemp,ye)

  this_dlnpdlnrho_s_ye = ffx(1,1)

end subroutine Interpolate_additional_lnderivs