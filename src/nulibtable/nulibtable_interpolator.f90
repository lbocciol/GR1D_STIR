!this takes rho,temp,ye,species and return eas over energy range
subroutine nulibtable_single_species_range_energy(iSpecie, eas, eas_n1, eas_n2)
  
  use GR1D_module, only: rho_gf, rho, temp, ye, n1, ghosts1
  use nulibtable
  implicit none

  real*8 :: log_rho(n1), log_temp(n1) !log versions
  integer, intent(in) :: iSpecie
  integer, intent(in) :: eas_n1,eas_n2
  real*8, intent(out) :: eas(eas_n1,eas_n2,n1)
  real*8 :: xeas(eas_n1)
  integer :: startindex,endindex

  if(size(eas,1).ne.nulibtable_number_groups) then
     stop "nulibtable_single_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_easvariables) then
     stop "nulibtable_single_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif

  startindex = (iSpecie-1)*nulibtable_number_groups+1
  endindex = startindex + nulibtable_number_groups - 1
  
  log_rho = log10(rho(k)/rho_gf)
  log_temp = log10(temp(k))

  xeas = 0.0d0
  call intp3d_man1_vectorized(log_rho, log_temp, ye, n1, &
    nulibtable_emissivities(:,:,:,startindex:endindex), &
    nulibtable_nrho, nulibtable_ntemp, nulibtable_n1e, eas_n1, &
    nulibtable_logrho, nulibtable_logtemp, nulibtable_ye, xeas)
  eas(:,1,:) = 10.0d0**xeas(:,:) !(nvars, n1)

  xeas = 0.0d0
  call intp3d_man1_vectorized(xlrho, xltemp, xye, n1, &
    nulibtable_absopacity(:,:,:,startindex:endindex), &
    nulibtable_nrho, nulibtable_ntemp, nulibtable_n1e, eas_n1, &
    nulibtable_logrho, nulibtable_logtemp, nulibtable_ye, xeas)
  eas(:,2,:) = 10.0d0**xeas(:,:) !(nvars, n1)

  xeas = 0.0d0
  call intp3d_man1_vectorized(xlrho, xltemp, xye, n1, &
    nulibtable_scatopacity(:,:,:,startindex:endindex), &
    nulibtable_nrho, nulibtable_ntemp, nulibtable_n1e, eas_n1, &
    nulibtable_logrho, nulibtable_logtemp, nulibtable_ye, xeas)
  eas(:,3,:) = 10.0d0**xeas(:,:) !(nvars, n1)

end subroutine nulibtable_single_species_range_energy

!this takes rho,temp,ye and return eas over energy and species range
subroutine nulibtable_range_species_range_energy(eas,eas_n1,eas_n2,eas_n3)

  use GR1D_module, only: rho_gf, rho, temp, ye, n1, ghosts1
  use nulibtable
  implicit none

  real*8 :: log_rho(n1), log_temp(n1) !log versions
  integer, intent(in) :: eas_n1,eas_n2,eas_n3
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3,n1)
  integer :: ins,ing
  real*8 :: xeas(eas_n1*eas_n2)
  integer :: index

  if(size(eas,1).ne.nulibtable_number_species) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif  
  if(size(eas,3).ne.nulibtable_number_easvariables) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif

  xlrho = log10(xrho)
  xltemp = log10(xtemp)

  xeas = 0.0d0
  call intp3d_man1_vectorized(log_rho, log_temp, ye, n1, &
    nulibtable_emissivities, &
    nulibtable_nrho, nulibtable_ntemp, nulibtable_n1e, eas_n1, &
    nulibtable_logrho, nulibtable_logtemp, nulibtable_ye, xeas)

  do k=ghosts1,n1-ghosts1-1
    do ins=1,nulibtable_number_species
      do ing=1,nulibtable_number_groups
          index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
          eas(ins,ing,1,k) = 10.0d0**xeas(index,k)
      enddo
    enddo
  enddo

  xeas = 0.0d0
  call intp3d_man1_vectorized(log_rho, log_temp, ye, n1, &
    nulibtable_absopacity, &
    nulibtable_nrho, nulibtable_ntemp, nulibtable_n1e, eas_n1, &
    nulibtable_logrho, nulibtable_logtemp, nulibtable_ye, xeas)

  do k=ghosts1,n1-ghosts1-1
    do ins=1,nulibtable_number_species
      do ing=1,nulibtable_number_groups
          index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
          eas(ins,ing,2,k) = 10.0d0**xeas(index,k)
      enddo
    enddo
  enddo

  xeas = 0.0d0
  call intp3d_man1_vectorized(log_rho, log_temp, ye, n1, &
    nulibtable_scatopacity, &
    nulibtable_nrho, nulibtable_ntemp, nulibtable_n1e, eas_n1, &
    nulibtable_logrho, nulibtable_logtemp, nulibtable_ye, xeas)

  do k=ghosts1,n1-ghosts1-1
    do ins=1,nulibtable_number_species
      do ing=1,nulibtable_number_groups
          index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
          eas(ins,ing,3,k) = 10.0d0**xeas(index,k)
      enddo
    enddo
  enddo

end subroutine nulibtable_range_species_range_energy

!===========================================================
! Inelastic scattering: species range over both energies
!===========================================================
subroutine nulibtable_inelastic_range_species_range_energy2(xtemp, xeta, eas, eas_n1, eas_n2, eas_n3, eas_n4)

  use nulibtable
  implicit none

  ! Inputs
  real*8, intent(in) :: xtemp, xeta
  integer, intent(in) :: eas_n1, eas_n2, eas_n3, eas_n4

  ! Outputs
  real*8, intent(out) :: eas(eas_n1, eas_n2, eas_n3, eas_n4)

  ! Locals
  real*8 :: xltemp, xleta
  integer :: ins, ing_in, ing_out, index
  real*8 :: xeas(eas_n1*eas_n2*(eas_n2+1)/2)
  real*8 :: energy_conversion

  !---------------------------------------------------------
  ! Checks
  if (size(eas,1) .ne. nulibtable_number_species) stop "inelastic: wrong dim (1)"
  if (size(eas,2) .ne. nulibtable_number_groups)  stop "inelastic: wrong dim (2)"
  if (size(eas,3) .ne. nulibtable_number_groups)  stop "inelastic: wrong dim (3)"
  if (size(eas,4) .ne. 2)                         stop "inelastic: wrong dim (4)"

  ! Log transforms
  xltemp = log10(xtemp)
  xleta  = log10(xeta)

  if (xltemp < nulibtable_logItemp_min .or. xltemp > nulibtable_logItemp_max) stop "inelastic: temp out of bounds"
  if (xleta  < nulibtable_logIeta_min  .or. xleta  > nulibtable_logIeta_max)  stop "inelastic: eta out of bounds"

  energy_conversion = 1.60217733d-6 * 5.59424238d-55

  !---------------------------------------------------------
  ! Phi0
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_Itable_Phi0, &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*eas_n2*(eas_n2+1)/2, &
       nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_in=1,nulibtable_number_groups
        do ing_out=1,ing_in
           index = index + 1
           eas(ins,ing_in,ing_out,1) = 10.0d0**xeas(index)
        enddo
     enddo
     do ing_in=1,nulibtable_number_groups
        do ing_out=ing_in+1,nulibtable_number_groups
           eas(ins,ing_in,ing_out,1) = exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in)) / &
                (xtemp*energy_conversion)) * eas(ins,ing_out,ing_in,1)
        enddo
     enddo
  enddo

  !---------------------------------------------------------
  ! Phi1 (interpolated as ratio phi1/phi0)
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_Itable_Phi1, &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*eas_n2*(eas_n2+1)/2, &
       nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_in=1,nulibtable_number_groups
        do ing_out=1,ing_in
           index = index + 1
           eas(ins,ing_in,ing_out,2) = xeas(index) * eas(ins,ing_in,ing_out,1)
        enddo
     enddo
     do ing_in=1,nulibtable_number_groups
        do ing_out=ing_in+1,nulibtable_number_groups
           eas(ins,ing_in,ing_out,2) = exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in)) / &
                (xtemp*energy_conversion)) * eas(ins,ing_out,ing_in,2)
        enddo
     enddo
  enddo

end subroutine nulibtable_inelastic_range_species_range_energy2


!===========================================================
! Inelastic scattering: single species over both energies
!===========================================================
subroutine nulibtable_inelastic_single_species_range_energy2(xtemp, xeta, iSpecie, eas, eas_n1, eas_n2, eas_n3)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xeta
  integer, intent(in) :: iSpecie, eas_n1, eas_n2, eas_n3
  real*8, intent(out):: eas(eas_n1, eas_n2, eas_n3)

  real*8 :: xltemp, xleta
  integer :: ing_in, ing_out, index
  integer :: startindex, endindex
  real*8 :: xeas(eas_n1*(eas_n1+1)/2)
  real*8 :: energy_conversion

  ! Checks
  if (size(eas,1) .ne. nulibtable_number_groups) stop "inelastic single: wrong dim (1)"
  if (size(eas,2) .ne. nulibtable_number_groups) stop "inelastic single: wrong dim (2)"
  if (size(eas,3) .ne. 2)                       stop "inelastic single: wrong dim (3)"

  xltemp = log10(xtemp)
  xleta  = log10(xeta)

  if (xltemp < nulibtable_logItemp_min .or. xltemp > nulibtable_logItemp_max) stop "inelastic single: temp out of bounds"
  if (xleta  < nulibtable_logIeta_min  .or. xleta  > nulibtable_logIeta_max)  stop "inelastic single: eta out of bounds"

  energy_conversion = 1.60217733d-6 * 5.59424238d-55

  startindex = (iSpecie-1)*nulibtable_number_groups*(nulibtable_number_groups+1)/2 + 1
  endindex   = startindex + nulibtable_number_groups*(nulibtable_number_groups+1)/2 - 1

  ! Phi0
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_Itable_Phi0(:,:,startindex:endindex), &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*(eas_n1+1)/2, nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ing_in=1,nulibtable_number_groups
     do ing_out=1,ing_in
        index = index + 1
        eas(ing_in,ing_out,1) = 10.0d0**xeas(index)
     enddo
  enddo
  do ing_in=1,nulibtable_number_groups
     do ing_out=ing_in+1,nulibtable_number_groups
        eas(ing_in,ing_out,1) = exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in)) / &
             (xtemp*energy_conversion)) * eas(ing_out,ing_in,1)
     enddo
  enddo

  ! Phi1
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_Itable_Phi1(:,:,startindex:endindex), &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*(eas_n1+1)/2, nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ing_in=1,nulibtable_number_groups
     do ing_out=1,ing_in
        index = index + 1
        eas(ing_in,ing_out,2) = xeas(index) * eas(ing_in,ing_out,1)
     enddo
  enddo
  do ing_in=1,nulibtable_number_groups
     do ing_out=ing_in+1,nulibtable_number_groups
        eas(ing_in,ing_out,2) = exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in)) / &
             (xtemp*energy_conversion)) * eas(ing_out,ing_in,2)
     enddo
  enddo

end subroutine nulibtable_inelastic_single_species_range_energy2


!===========================================================
! e+e- annihilation: species range over both energies
!===========================================================
subroutine nulibtable_epannihil_range_species_range_energy2(xtemp, xeta, eas, eas_n1, eas_n2, eas_n3, eas_n4)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xeta
  integer, intent(in) :: eas_n1, eas_n2, eas_n3, eas_n4
  real*8, intent(out):: eas(eas_n1, eas_n2, eas_n3, eas_n4)

  real*8 :: xltemp, xleta
  integer :: ins, ing_this, ing_that, index
  real*8 :: xeas(eas_n1*eas_n2*eas_n2*2)

  ! Checks
  if (size(eas,1) .ne. nulibtable_number_species) stop "epannihil: wrong dim (1)"
  if (size(eas,2) .ne. nulibtable_number_groups)  stop "epannihil: wrong dim (2)"
  if (size(eas,3) .ne. nulibtable_number_groups)  stop "epannihil: wrong dim (3)"
  if (size(eas,4) .ne. 4)                         stop "epannihil: wrong dim (4)"

  xltemp = log10(xtemp)
  xleta  = log10(xeta)

  if (xltemp < nulibtable_logItemp_min .or. xltemp > nulibtable_logItemp_max) stop "epannihil: temp out of bounds"
  if (xleta  < nulibtable_logIeta_min  .or. xleta  > nulibtable_logIeta_max)  stop "epannihil: eta out of bounds"

  ! Phi0
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_epannihiltable_Phi0, &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*eas_n2*eas_n2*2, &
       nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_this=1,nulibtable_number_groups
        do ing_that=1,nulibtable_number_groups
           index = index + 1
           eas(ins,ing_this,ing_that,1) = 10.0d0**xeas(index)
           index = index + 1
           eas(ins,ing_this,ing_that,2) = 10.0d0**xeas(index)
        enddo
     enddo
  enddo

  ! Phi1
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_epannihiltable_Phi1, &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*eas_n2*eas_n2*2, &
       nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_this=1,nulibtable_number_groups
        do ing_that=1,nulibtable_number_groups
           index = index + 1
           eas(ins,ing_this,ing_that,3) = xeas(index) * eas(ins,ing_this,ing_that,1)
           index = index + 1
           eas(ins,ing_this,ing_that,4) = xeas(index) * eas(ins,ing_this,ing_that,2)
        enddo
     enddo
  enddo

end subroutine nulibtable_epannihil_range_species_range_energy2


!===========================================================
! e+e- annihilation: single species over both energies
!===========================================================
subroutine nulibtable_epannihil_single_species_range_energy2(xtemp, xeta, iSpecie, eas, eas_n1, eas_n2, eas_n3)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xeta
  integer, intent(in) :: iSpecie, eas_n1, eas_n2, eas_n3
  real*8, intent(out):: eas(eas_n1, eas_n2, eas_n3)

  real*8 :: xltemp, xleta
  integer :: ing_this, ing_that, index
  integer :: startindex, endindex
  real*8 :: xeas(eas_n1*eas_n2*2)

  ! Checks
  if (size(eas,1) .ne. nulibtable_number_groups) stop "epannihil single: wrong dim (1)"
  if (size(eas,2) .ne. nulibtable_number_groups) stop "epannihil single: wrong dim (2)"
  if (size(eas,3) .ne. 4)                       stop "epannihil single: wrong dim (3)"

  xltemp = log10(xtemp)
  xleta  = log10(xeta)

  if (xltemp < nulibtable_logItemp_min .or. xltemp > nulibtable_logItemp_max) stop "epannihil single: temp out of bounds"
  if (xleta  < nulibtable_logIeta_min  .or. xleta  > nulibtable_logIeta_max)  stop "epannihil single: eta out of bounds"

  startindex = (iSpecie-1)*nulibtable_number_groups*nulibtable_number_groups*2 + 1
  endindex   = startindex + nulibtable_number_groups*nulibtable_number_groups*2 - 1

  ! Phi0
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_epannihiltable_Phi0(:,:,startindex:endindex), &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*eas_n2*2, nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ing_this=1,nulibtable_number_groups
     do ing_that=1,nulibtable_number_groups
        index = index + 1
        eas(ing_this,ing_that,1) = 10.0d0**xeas(index)
        index = index + 1
        eas(ing_this,ing_that,2) = 10.0d0**xeas(index)
     enddo
  enddo

  ! Phi1
  call intp2d_man1_mod(xltemp, xleta, xeas, nulibtable_epannihiltable_Phi1(:,:,startindex:endindex), &
       nulibtable_nItemp, nulibtable_nIeta, eas_n1*eas_n2*2, nulibtable_logItemp, nulibtable_logIeta)

  index = 0
  do ing_this=1,nulibtable_number_groups
     do ing_that=1,nulibtable_number_groups
        index = index + 1
        eas(ing_this,ing_that,3) = xeas(index) * eas(ing_this,ing_that,1)
        index = index + 1
        eas(ing_this,ing_that,4) = xeas(index) * eas(ing_this,ing_that,2)
     enddo
  enddo

end subroutine nulibtable_epannihil_single_species_range_energy2
