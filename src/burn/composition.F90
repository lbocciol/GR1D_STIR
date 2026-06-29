! composition.F90
!
! The GR1D-level nuclear composition: the set of species GR1D tracks per zone in
! Yion(:,i), which is the reaction network's species set OPTIONALLY augmented with
! free nucleons (n, p).
!
! WHY THIS EXISTS
!   An alpha-chain network (e.g. aprox13) contains only N=Z nuclei, so it can only
!   represent Ye = 0.5 and its NSE charge constraint sum Z_i Y_i = Ye is not
!   independent of baryon number -- the NSE solve is singular for Ye /= 0.5.  Adding
!   free n and p as TRACKED-but-INERT species gives the composition the two extra
!   degrees of freedom needed to carry Ye /= 0.5: NSE becomes well-posed, and the
!   network simply never touches the appended nucleon slots (Step.F90 advances only
!   Yion(1:nspec_net,i), so n,p stay frozen).
!
!   A network that already contains n and p (e.g. aprox19, with weak rates) needs no
!   augmentation: set track_free_nucleons=.false. and this module is a pass-through
!   (nspec == nspec_net), so the rest of the code is unchanged and the nucleons burn
!   normally.  That is the whole point of routing the composition through here.
!
! LAYOUT
!   Slots 1..nspec_net  : the network species, in pynet's order (the network only
!                         ever sees these).
!   Slots nspec_net+1.. : appended free nucleons (n then p), present only when
!                         track_free_nucleons=.true.
!
! UNITS
!   aion, zion           : dimensionless mass and proton numbers
!   nuclei_binding_energy : nuclear binding energy B_i [MeV] (FLASH `bion`), same
!                          ordering; free nucleons (n,p) have B = 0.  Consumed by the
!                          NSE Saha solver (nse.F90).  NOTE: the composite-EOS energy
!                          zero-point is NOT computed here anymore -- it is matched
!                          numerically against nuc_eos in build_energy_offset_OttEOS (eos.F90).

module composition

  use pynet, only: pyn_nspec => nspec, pyn_aion => aion, pyn_zion => zion

  implicit none
  private

  ! Full composition (network species, plus appended n,p when tracked).
  integer, save :: nspec     = 0        ! full species count
  integer, save :: nspec_net = 0        ! network species count (= pynet nspec)
  real(8), allocatable, save :: aion(:), zion(:)
  real(8), allocatable, save :: nuclei_binding_energy(:), nuclei_mass_excess(:)

  public :: nspec, nspec_net, aion, zion, nuclei_binding_energy, nuclei_mass_excess
  public :: composition_init, composition_set_free_nucleons

  ! Define a structure to hold the data
  type :: isotope
      character(len=2) :: symbol
      integer          :: Z
      integer          :: A
      real(8)          :: mass_excess    ! in keV
      real(8)          :: binding_energy ! Total Binding Energy in keV
  end type isotope

  ! Array constructor for the 22 requested isotopes (sorted by Z, then A)
  ! Data based on AME2020 values
  type(isotope), parameter :: ISOTOPES(22) = [ &
      isotope('He',  2,  3,  14931.218,   7718.042), &
      isotope('He',  2,  4,   2424.915,  28295.663), &
      isotope('C ',  6, 12,      0.000,  92161.734), &
      isotope('O ',  8, 16,  -4737.002, 127619.314), &
      isotope('Ne', 10, 20,  -7041.931, 160644.821), &
      isotope('Mg', 12, 24, -13933.567, 198257.035), &
      isotope('Si', 14, 28, -21492.794, 236536.840), &
      isotope('S ', 16, 32, -26015.545, 271780.169), &
      isotope('S ', 16, 36, -30820.180, 308870.076), &
      isotope('Ar', 18, 36, -30231.541, 306716.743), &
      isotope('Ca', 20, 40, -34846.111, 342051.891), &
      isotope('Ti', 22, 44, -37548.800, 375475.158), &
      isotope('Ti', 22, 50, -48332.900, 434687.166), &
      isotope('Cr', 24, 48, -42817.500, 411464.436), &
      isotope('Fe', 26, 52, -48330.600, 447698.114), &
      isotope('Fe', 26, 54, -56253.700, 471763.850), &
      isotope('Fe', 26, 56, -60605.400, 492258.186), &
      isotope('Fe', 26, 58, -62153.500, 509948.922), &
      isotope('Fe', 26, 60, -61405.600, 525343.658), &
      isotope('Fe', 26, 62, -58872.400, 538953.094), &
      isotope('Ni', 28, 56, -53907.300, 483995.392), &
      isotope('Ni', 28, 62, -66746.100, 545262.100)  &
  ]

contains

  ! ---------------------------------------------------------------------------
  ! Build the full composition from the pynet network.  Called once at startup,
  ! BEFORE allocate_vars (which sizes Yion with nspec) and nse_init (which uses the
  ! full species arrays).
  !   track_free_nucleons = .true.  -> append free n,p (network must NOT already
  !                                    contain them)
  !                       = .false. -> full == network (use this for aprox19 etc.)
  subroutine composition_init(track_free_nucleons)
    logical, intent(in) :: track_free_nucleons
    integer :: k, nadd
    logical :: has_n, has_p

    nspec_net = pyn_nspec

    ! detect free nucleons (A==1) already present in the network
    has_n = .false.; has_p = .false.
    do k = 1, nspec_net
      if (nint(pyn_aion(k)) == 1 .and. nint(pyn_zion(k)) == 0) has_n = .true.
      if (nint(pyn_aion(k)) == 1 .and. nint(pyn_zion(k)) == 1) has_p = .true.
    end do

    if (track_free_nucleons) then
      if (has_n .or. has_p) then
        write(*,*) "composition_init: track_free_nucleons=.true. but the network ", &
                   "already contains free n and/or p -- set it .false. for this network."
        stop "Aborting!"
      end if
      nadd = 2
    else
      nadd = 0
    end if

    nspec = nspec_net + nadd
    allocate(aion(nspec), zion(nspec))
    allocate(nuclei_binding_energy(nspec), nuclei_mass_excess(nspec))

    aion(1:nspec_net) = pyn_aion
    zion(1:nspec_net) = pyn_zion

    if (nadd == 2) then
      ! appended free nucleons: neutron then proton.  A free nucleon is unbound,
      ! so its binding energy is 0 (it still contributes to Ye via its protons).
      aion(nspec_net+1) = 1.0d0; zion(nspec_net+1) = 0.0d0
      aion(nspec_net+2) = 1.0d0; zion(nspec_net+2) = 1.0d0
    end if

    call Get_BE_and_ME_given_species( &
          aion, zion, &
          nuclei_binding_energy, nuclei_mass_excess, nspec)

  end subroutine composition_init

  ! ---------------------------------------------------------------------------
  ! Fill the appended free-nucleon slots (n,p) of a composition vector so the FULL
  ! composition carries every baryon and reproduces the requested electron fraction:
  !
  !     sum_i Z_i Y_i = ye ,   sum_i A_i Y_i = 1 .
  !
  ! On entry Y(1:nspec_net) holds the network species (from the progenitor compo
  ! file, which typically sums to < 1 in mass because the rest is free nucleons);
  ! the network's charge is z_net = sum Z_i Y_i and its baryon (mass) fraction is
  ! f_net = sum A_i Y_i.  The free nucleons (A=1, so Y=X) take up the remainder:
  !     Y_p = ye - z_net ,   Y_n = (1 - f_net) - Y_p .
  ! This is what lets an all-N=Z network represent Ye/=0.5 and gives the EOS a
  ! physically complete composition (e.g. nucleon-dominated matter at high T,rho).
  !
  ! No-op when there are no appended nucleons (nspec == nspec_net, e.g. aprox19,
  ! whose progenitor file already supplies n,p).  Tiny negatives from rounding (or
  ! a network that already nearly fills the baryon budget at Ye/=0.5) are clamped.
  subroutine composition_set_free_nucleons(Y, ye)
    real(8), intent(inout) :: Y(nspec)   ! in: network slots filled; out: n,p set
    real(8), intent(in)    :: ye
    real(8) :: f_net, z_net, yp, yn

    if (nspec == nspec_net) return       ! network already includes n,p

    f_net = sum(aion(1:nspec_net) * Y(1:nspec_net))   ! network mass fraction
    z_net = sum(zion(1:nspec_net) * Y(1:nspec_net))   ! network charge

    yp = max(ye - z_net,          0.0d0)              ! free protons
    yn = max(1.0d0 - f_net - yp,  0.0d0)              ! remaining baryons -> neutrons

    Y(nspec_net+1) = yn                  ! neutron slot (see composition_init order)
    Y(nspec_net+2) = yp                  ! proton  slot
  end subroutine composition_set_free_nucleons

  subroutine Get_BE_and_ME_given_species(A, Z, BE, ME, nisos)
  
    integer, intent(in)  :: nisos
    real(8), intent(in)  :: A(nisos), Z(nisos)
    real(8), intent(out) :: BE(nisos), ME(nisos)

    integer :: iIso, i
    logical :: found

    ! 1. Initialize arrays to 0.0d0
    ! This ensures free nucleons (or unknown species) get a default Binding Energy of 0
    BE = 0.0d0
    ME = 0.0d0

    ! 2. Loop through the species requested by the solver
    do iIso = 1, nisos
        found = .false.

        ! 3. Search the ISOTOPES array for matching A and Z
        do i = 1, size(ISOTOPES)
            if (ISOTOPES(i)%A == INT(A(iIso)) .and. ISOTOPES(i)%Z == INT(Z(iIso))) then
                ! Extract values and convert from keV to MeV
                BE(iIso) = ISOTOPES(i)%binding_energy / 1000.0d0
                ME(iIso) = ISOTOPES(i)%mass_excess / 1000.0d0
                found = .true.
                exit
            end if
        end do

        ! 4. Handle exact Mass Excess for free nucleons (Optional but recommended for NSE)
        ! Even though BE is 0, free nucleons have a non-zero Mass Excess.
        if (.not. found) then
            if (Z(iIso) == 1 .and. A(iIso) == 1) then
                ME(iIso) = 7.288971d0  ! Proton Mass Excess in MeV
            else if (Z(iIso) == 0 .and. A(iIso) == 1) then
                ME(iIso) = 8.071318d0  ! Neutron Mass Excess in MeV
            end if
        end if
    end do

  end subroutine Get_BE_and_ME_given_species

end module composition


