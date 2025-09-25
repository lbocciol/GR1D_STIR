!--------------------------------------------------------------------
!  Fixed/optimized versions of intp2d_manyvectorized and
!  intp3d_many_vectorized
!--------------------------------------------------------------------
SUBROUTINE intp2d_manyvectorized(Y1, Y2, f, ft, nX1, nX2, nvars, X1t, X2t)
  implicit none

  integer, intent(in) :: nX1, nX2, nvars
  real*8, intent(in) :: ft(nX1,nX2,nvars)
  real*8, intent(in) :: Y1(:), Y2(:), X1t(nX1), X2t(nX2)
  real*8, intent(out) :: f(nvars, size(Y1))

  integer :: nY
  integer :: iY, iv
  real*8 :: dx, dy, dxi, dyi, dxyi
  real*8, allocatable :: delx(:), dely(:)
  integer, allocatable :: iY1(:), iY2(:)
  real*8, allocatable :: fh_all(:,:,:)  ! nY x nvars x 4

  nY = size(Y1)

  ! Determine spacing (assumes equidistant)
  dx = (X1t(nX1) - X1t(1)) / dble(nX1 - 1)
  dy = (X2t(nX2) - X2t(1)) / dble(nX2 - 1)
  dxi = 1.0d0 / dx
  dyi = 1.0d0 / dy
  dxyi = dxi * dyi

  allocate(delx(nY), dely(nY))
  allocate(iY1(nY), iY2(nY))
  allocate(fh_all(nY, nvars, 4))

  !---------------------------------------------------------
  ! Precompute indices and deltas for each query point
  ! (parallelizable)
  !---------------------------------------------------------
  !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(iY) &
  !$OMP& SHARED(Y1,Y2,X1t,X2t,iY1,iY2,delx,dely,dxi,dyi,nX1,nX2)
  do iY = 1, nY
     iY1(iY) = 2 + int( (Y1(iY) - X1t(1) - 1.0d-10) * dxi )
     iY2(iY) = 2 + int( (Y2(iY) - X2t(1) - 1.0d-10) * dyi )

     iY1(iY) = max(2, min(iY1(iY), nX1))
     iY2(iY) = max(2, min(iY2(iY), nX2))

     delx(iY) = X1t(iY1(iY)) - Y1(iY)
     dely(iY) = X2t(iY2(iY)) - Y2(iY)
  end do
  !$OMP END PARALLEL DO

  !---------------------------------------------------------
  ! Precompute surrounding table values fh_all for each point
  ! Collapse across (iY, iv) so threads can share work
  !---------------------------------------------------------
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(iY,iv) &
  !$OMP& SHARED(ft,fh_all,iY1,iY2)
  do iY = 1, nY
    do iv = 1, nvars
      fh_all(iY,iv,1) = ft(iY1(iY)  , iY2(iY)  , iv)
      fh_all(iY,iv,2) = ft(iY1(iY)-1, iY2(iY)  , iv)
      fh_all(iY,iv,3) = ft(iY1(iY)  , iY2(iY)-1, iv)
      fh_all(iY,iv,4) = ft(iY1(iY)-1, iY2(iY)-1, iv)
    end do
  end do
  !$OMP END PARALLEL DO

  !---------------------------------------------------------
  ! Compute interpolated values using precomputed fh_all and deltas
  ! Collapse across (iY, iv)
  !---------------------------------------------------------
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(iY,iv) &
  !$OMP& SHARED(f,fh_all,delx,dely,dxi,dyi,dxyi)
  do iY = 1, nY
    do iv = 1, nvars
      f(iv,iY) = fh_all(iY,iv,1) + &
                  dxi*(fh_all(iY,iv,2)-fh_all(iY,iv,1))*delx(iY) + &
                  dyi*(fh_all(iY,iv,3)-fh_all(iY,iv,1))*dely(iY) + &
                  dxyi*(fh_all(iY,iv,4)-fh_all(iY,iv,2)-fh_all(iY,iv,3)+fh_all(iY,iv,1))*delx(iY)*dely(iY)
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(delx, dely, iY1, iY2, fh_all)

END SUBROUTINE intp2d_manyvectorized

!--------------------------------------------------------------------
!  3D version — cleaned up to truly use precomputed index/delta arrays
!--------------------------------------------------------------------
SUBROUTINE intp3d_many_vectorized(Y1, Y2, Y3, nY, ft, nX1, nX2, nX3, nvars, X1t, X2t, X3t, f)
  implicit none

  integer, intent(in) :: nY, nX1, nX2, nX3, nvars
  real*8, intent(in) :: Y1(nY), Y2(nY), Y3(nY)
  real*8, intent(in) :: ft(nX1,nX2,nX3,nvars)
  real*8, intent(in) :: X1t(nX1), X2t(nX2), X3t(nX3)
  real*8, intent(out):: f(nvars,nY)

  integer :: iY, iv
  integer, allocatable :: iY1(:), iY2(:), iY3(:)
  real*8, allocatable :: delx(:), dely(:), delz(:)
  real*8 :: dx, dy, dz, dxi, dyi, dzi, dxyi, dxzi, dyzi, dxyzi
  real*8, allocatable :: fh_all(:,:,:)
  real*8, allocatable :: a(:,:,:)

  !---------------------------------------------------------
  ! Compute spacing and inverse spacing (equidistant table)
  !---------------------------------------------------------
  dx    = (X1t(nX1) - X1t(1)) / dble(nX1 - 1)
  dy    = (X2t(nX2) - X2t(1)) / dble(nX2 - 1)
  dz    = (X3t(nX3) - X3t(1)) / dble(nX3 - 1)

  dxi   = 1.0d0 / dx
  dyi   = 1.0d0 / dy
  dzi   = 1.0d0 / dz
  dxyi  = dxi * dyi
  dxzi  = dxi * dzi
  dyzi  = dyi * dzi
  dxyzi = dxi * dyi * dzi

  allocate(iY1(nY), iY2(nY), iY3(nY))
  allocate(delx(nY), dely(nY), delz(nY))
  allocate(fh_all(nY, nvars, 8))
  allocate(a(nY, nvars, 8))

  !---------------------------------------------------------
  ! Determine table indices and deltas for all points (parallel)
  !---------------------------------------------------------
  !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(iY) &
  !$OMP& SHARED(Y1,Y2,Y3,X1t,X2t,X3t,iY1,iY2,iY3,delx,dely,delz,dxi,dyi,dzi,nX1,nX2,nX3)
  do iY = 1, nY
     iY1(iY) = 2 + int( (Y1(iY) - X1t(1) - 1.0d-10) * dxi )
     iY2(iY) = 2 + int( (Y2(iY) - X2t(1) - 1.0d-10) * dyi )
     iY3(iY) = 2 + int( (Y3(iY) - X3t(1) - 1.0d-10) * dzi )

     iY1(iY) = max(2, min(iY1(iY), nX1))
     iY2(iY) = max(2, min(iY2(iY), nX2))
     iY3(iY) = max(2, min(iY3(iY), nX3))

     delx(iY) = X1t(iY1(iY)) - Y1(iY)
     dely(iY) = X2t(iY2(iY)) - Y2(iY)
     delz(iY) = X3t(iY3(iY)) - Y3(iY)
  end do
  !$OMP END PARALLEL DO

  !---------------------------------------------------------
  ! Precompute surrounding points fh_all for all variables
  ! Collapse across (iY, iv)
  !---------------------------------------------------------
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(iY,iv) &
  !$OMP& SHARED(ft,fh_all,iY1,iY2,iY3)
  do iY = 1, nY
    do iv = 1, nvars
      fh_all(iY,iv,1) = ft(iY1(iY)  , iY2(iY)  , iY3(iY)  , iv)
      fh_all(iY,iv,2) = ft(iY1(iY)-1, iY2(iY)  , iY3(iY)  , iv)
      fh_all(iY,iv,3) = ft(iY1(iY)  , iY2(iY)-1, iY3(iY)  , iv)
      fh_all(iY,iv,4) = ft(iY1(iY)  , iY2(iY)  , iY3(iY)-1, iv)
      fh_all(iY,iv,5) = ft(iY1(iY)-1, iY2(iY)-1, iY3(iY)  , iv)
      fh_all(iY,iv,6) = ft(iY1(iY)-1, iY2(iY)  , iY3(iY)-1, iv)
      fh_all(iY,iv,7) = ft(iY1(iY)  , iY2(iY)-1, iY3(iY)-1, iv)
      fh_all(iY,iv,8) = ft(iY1(iY)-1, iY2(iY)-1, iY3(iY)-1, iv)
    end do
  end do
  !$OMP END PARALLEL DO

  !---------------------------------------------------------
  ! Compute interpolation coefficients and evaluate
  ! Collapse across (iY, iv)
  !---------------------------------------------------------
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(iY,iv) &
  !$OMP& SHARED(a,fh_all,delx,dely,delz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi,f)
  do iY = 1, nY
    do iv = 1, nvars
      a(iY,iv,1) = fh_all(iY,iv,1)
      a(iY,iv,2) = dxi   * (fh_all(iY,iv,2) - fh_all(iY,iv,1))
      a(iY,iv,3) = dyi   * (fh_all(iY,iv,3) - fh_all(iY,iv,1))
      a(iY,iv,4) = dzi   * (fh_all(iY,iv,4) - fh_all(iY,iv,1))
      a(iY,iv,5) = dxyi  * (fh_all(iY,iv,5) - fh_all(iY,iv,2) - fh_all(iY,iv,3) + fh_all(iY,iv,1))
      a(iY,iv,6) = dxzi  * (fh_all(iY,iv,6) - fh_all(iY,iv,2) - fh_all(iY,iv,4) + fh_all(iY,iv,1))
      a(iY,iv,7) = dyzi  * (fh_all(iY,iv,7) - fh_all(iY,iv,3) - fh_all(iY,iv,4) + fh_all(iY,iv,1))
      a(iY,iv,8) = dxyzi * (fh_all(iY,iv,8) - fh_all(iY,iv,1) + fh_all(iY,iv,2) + fh_all(iY,iv,3) + fh_all(iY,iv,4) &
                 - fh_all(iY,iv,5) - fh_all(iY,iv,6) - fh_all(iY,iv,7))

      f(iv,iY) = a(iY,iv,1) + a(iY,iv,2) * delx(iY) + a(iY,iv,3) * dely(iY) + a(iY,iv,4) * delz(iY) &
                 + a(iY,iv,5) * delx(iY) * dely(iY) + a(iY,iv,6) * delx(iY) * delz(iY) &
                 + a(iY,iv,7) * dely(iY) * delz(iY) + a(iY,iv,8) * delx(iY) * dely(iY) * delz(iY)
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(iY1, iY2, iY3, delx, dely, delz, fh_all, a)

END SUBROUTINE intp3d_many_vectorized
