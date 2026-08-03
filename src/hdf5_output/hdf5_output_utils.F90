!-*-f90-*-
module hdf5_output_utils
  !! HDF5 Output Utilities Module
  !! Provides wrapper functions for writing GR1D_STIR output to HDF5 files
  !! Two files are created: xg.h5 (grid profiles, split into /hydro and /M1
  !! groups) and dat.h5 (scalar time series); with HAVE_BURN a third file
  !! burn.h5 holds all nuclear-burning output.
  !! Files are kept open within output scope for efficiency

  use hdf5
  implicit none

  private

  ! Public subroutines
  public :: hdf5_initialize
  public :: hdf5_create_output_files
  public :: hdf5_write_root_dataset_1d
  public :: hdf5_append_grid_data_1d
  public :: hdf5_append_xg_scalar
  public :: hdf5_append_scalar
  public :: hdf5_append_scalar_array
  public :: hdf5_finalize
  public :: hdf5_write_metadata
  public :: hdf5_output_init
  public :: hdf5_open_xg_file
  public :: hdf5_close_xg_file
  public :: hdf5_open_dat_file
  public :: hdf5_close_dat_file
#ifdef HAVE_BURN
  public :: hdf5_open_burn_file
  public :: hdf5_close_burn_file
  public :: hdf5_append_burn_scalar
  public :: hdf5_append_burn_1d
  public :: hdf5_append_burn_2d
  public :: hdf5_write_burn_root_dataset_1d
  public :: hdf5_write_burn_metadata
#endif

  ! Module variables
  character(len=1024), save :: xg_file_path = ""
  character(len=1024), save :: dat_file_path = ""

  ! File/group handle storage for persistent open files
  integer(HID_T), save :: xg_file_id = -1
  integer(HID_T), save :: xg_hydro_group_id = -1
  integer(HID_T), save :: xg_M1_group_id = -1
  integer(HID_T), save :: dat_file_id = -1
  integer(HID_T), save :: dat_scalars_group_id = -1
  logical, save :: xg_file_open = .false.
  logical, save :: dat_file_open = .false.

#ifdef HAVE_BURN
  character(len=1024), save :: burn_file_path = ""
  integer(HID_T), save :: burn_file_id = -1
  integer(HID_T), save :: burn_group_id = -1   ! the /fields group
  logical, save :: burn_file_open = .false.
#endif

  ! Compression settings (GZIP removed - using regular HDF5 files)
  ! integer, parameter :: GZIP_LEVEL = 4

contains

  subroutine hdf5_initialize()

    integer :: error

    ! THIS IS THE CRITICAL LINE
    call h5open_f(error)

    if (error /= 0) then
      write(*,*) "ERROR: Failed to initialize Fortran HDF5 library!"
      stop
    endif

  end subroutine hdf5_initialize

  subroutine hdf5_create_output_files(outdir)
    !! Create xg.h5 and dat.h5 files (and burn.h5 when compiled with HAVE_BURN)
    character(len=*), intent(in) :: outdir
    integer :: error
    integer(HID_T) :: file_id, fcpl_id, group_id

    ! Store file paths for later use
    xg_file_path = trim(adjustl(outdir))//"/xg.h5"
    dat_file_path = trim(adjustl(outdir))//"/dat.h5"

    ! Create file creation property list
    call h5pcreate_f(H5P_FILE_CREATE_F, fcpl_id, error)

    ! Create xg.h5 file
    call h5fcreate_f(trim(xg_file_path), H5F_ACC_TRUNC_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create xg.h5"
      stop
    endif

        ! Create /metadata group
    call h5gcreate_f(file_id, "/metadata", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /metadata group"
      stop
    endif
    call h5gclose_f(group_id, error)

    ! Create /hydro group (append-style time-series datasets, grid variables)
    call h5gcreate_f(file_id, "/hydro", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /hydro group"
      stop
    endif
    call h5gclose_f(group_id, error)

    ! Create /M1 group (neutrino variables and spectra)
    call h5gcreate_f(file_id, "/M1", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /M1 group"
      stop
    endif
    call h5gclose_f(group_id, error)

    call h5fclose_f(file_id, error)

    ! Create dat.h5 file with /scalars and /metadata groups
    call h5fcreate_f(trim(dat_file_path), H5F_ACC_TRUNC_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create dat.h5"
      stop
    endif

    ! Create /scalars group
    call h5gcreate_f(file_id, "/scalars", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /scalars group"
      stop
    endif

    ! Create /metadata group
    call h5gcreate_f(file_id, "/metadata", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /metadata group"
      stop
    endif

    call h5fclose_f(file_id, error)

#ifdef HAVE_BURN
    ! Create burn.h5 file with /fields and /metadata groups
    burn_file_path = trim(adjustl(outdir))//"/burn.h5"

    call h5fcreate_f(trim(burn_file_path), H5F_ACC_TRUNC_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create burn.h5"
      stop
    endif

    call h5gcreate_f(file_id, "/metadata", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /metadata group in burn.h5"
      stop
    endif
    call h5gclose_f(group_id, error)

    call h5gcreate_f(file_id, "/fields", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /fields group in burn.h5"
      stop
    endif
    call h5gclose_f(group_id, error)

    call h5fclose_f(file_id, error)
#endif

    call h5pclose_f(fcpl_id, error)

  end subroutine hdf5_create_output_files

  subroutine hdf5_create_or_open_group(file_id, group_name, group_id)
    !! Open a group under file_id, creating it first if it doesn't exist yet
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: group_name
    integer(HID_T), intent(out) :: group_id
    integer :: error
    logical :: group_exists

    call h5lexists_f(file_id, trim(group_name), group_exists, error)

    if (.not. group_exists) then
      call h5gcreate_f(file_id, trim(group_name), group_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to create ", trim(group_name), " group"
        return
      endif
    else
      call h5gopen_f(file_id, trim(group_name), group_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to open ", trim(group_name), " group"
        return
      endif
    endif

  end subroutine hdf5_create_or_open_group

  subroutine hdf5_open_xg_file()
    !! Open xg.h5 file and the /hydro and /M1 groups for persistent writing
    integer :: error

    ! Open file if not already open
    if (.not. xg_file_open) then
      call h5fopen_f(trim(xg_file_path), H5F_ACC_RDWR_F, xg_file_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to open xg.h5 for writing"
        return
      endif
      xg_file_open = .true.
    endif

    call hdf5_create_or_open_group(xg_file_id, "/hydro", xg_hydro_group_id)
    call hdf5_create_or_open_group(xg_file_id, "/M1", xg_M1_group_id)

  end subroutine hdf5_open_xg_file

  subroutine hdf5_close_xg_file()
    !! Close xg.h5 file and groups
    integer :: error

    if (xg_file_open) then
      if (xg_hydro_group_id >= 0) then
        call h5gclose_f(xg_hydro_group_id, error)
        xg_hydro_group_id = -1
      endif
      if (xg_M1_group_id >= 0) then
        call h5gclose_f(xg_M1_group_id, error)
        xg_M1_group_id = -1
      endif
      call h5fclose_f(xg_file_id, error)
      xg_file_id = -1
      xg_file_open = .false.
    endif

  end subroutine hdf5_close_xg_file

#ifdef HAVE_BURN
  subroutine hdf5_open_burn_file()
    !! Open burn.h5 file and the /fields group for persistent writing
    integer :: error

    ! Open file if not already open
    if (.not. burn_file_open) then
      call h5fopen_f(trim(burn_file_path), H5F_ACC_RDWR_F, burn_file_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to open burn.h5 for writing"
        return
      endif
      burn_file_open = .true.
    endif

    call hdf5_create_or_open_group(burn_file_id, "/fields", burn_group_id)

  end subroutine hdf5_open_burn_file

  subroutine hdf5_close_burn_file()
    !! Close burn.h5 file and group
    integer :: error

    if (burn_file_open) then
      if (burn_group_id >= 0) then
        call h5gclose_f(burn_group_id, error)
        burn_group_id = -1
      endif
      call h5fclose_f(burn_file_id, error)
      burn_file_id = -1
      burn_file_open = .false.
    endif

  end subroutine hdf5_close_burn_file
#endif

  subroutine hdf5_open_dat_file()
    !! Open dat.h5 file and get /scalars group for persistent writing
    integer :: error

    ! Open file if not already open
    if (.not. dat_file_open) then
      call h5fopen_f(trim(dat_file_path), H5F_ACC_RDWR_F, dat_file_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to open dat.h5 for writing"
        return
      endif

      ! Open /scalars group
      call h5gopen_f(dat_file_id, "/scalars", dat_scalars_group_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to open /scalars group"
        call h5fclose_f(dat_file_id, error)
        return
      endif

      dat_file_open = .true.
    endif

  end subroutine hdf5_open_dat_file

  subroutine hdf5_close_dat_file()
    !! Close dat.h5 file and scalars group
    integer :: error

    if (dat_file_open) then
      if (dat_scalars_group_id >= 0) then
        call h5gclose_f(dat_scalars_group_id, error)
        dat_scalars_group_id = -1
      endif
      call h5fclose_f(dat_file_id, error)
      dat_file_id = -1
      dat_file_open = .false.
    endif

  end subroutine hdf5_close_dat_file

  subroutine hdf5_write_root_dataset_1d_at(file_path, varname, data, n1)
    !! Write 1D grid variable directly to the root of the given HDF5 file
    !! Opens file, creates dataset, writes data, closes file
    character(len=*), intent(in) :: file_path
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1
    real(kind=8), intent(in) :: data(n1)

    integer :: error
    integer(HID_T) :: file_id, dset_id, dspace_id, dcpl_id
    integer(HSIZE_T) :: dims(1)

    ! Open file
    call h5fopen_f(trim(file_path), H5F_ACC_RDWR_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to open ", trim(file_path), " for writing"
      return
    endif

    dims(1) = n1

    ! Create dataspace
    call h5screate_simple_f(1, dims, dspace_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create dataspace for ", trim(varname)
      call h5fclose_f(file_id, error)
      return
    endif

    ! Create dataset creation property list (no compression)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)

    ! Create dataset directly at the root level using file_id
    call h5dcreate_f(file_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, error, dcpl_id)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create dataset for ", trim(varname)
      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)
      call h5fclose_f(file_id, error)
      return
    endif

    ! Write data
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to write data for ", trim(varname)
    endif

    ! Close dataset and resources
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
    call h5pclose_f(dcpl_id, error)
    call h5fclose_f(file_id, error)

  end subroutine hdf5_write_root_dataset_1d_at

  subroutine hdf5_write_root_dataset_1d(varname, data, n1)
    !! Write 1D grid variable directly to the root of xg.h5
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1
    real(kind=8), intent(in) :: data(n1)

    call hdf5_write_root_dataset_1d_at(xg_file_path, varname, data, n1)

  end subroutine hdf5_write_root_dataset_1d

#ifdef HAVE_BURN
  subroutine hdf5_write_burn_root_dataset_1d(varname, data, n1)
    !! Write 1D variable directly to the root of burn.h5
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1
    real(kind=8), intent(in) :: data(n1)

    call hdf5_write_root_dataset_1d_at(burn_file_path, varname, data, n1)

  end subroutine hdf5_write_burn_root_dataset_1d
#endif

  subroutine hdf5_append_grid_data_1d(group, varname, data, n1)
    !! Append one row (one dump) of a grid variable to /hydro or /M1 in xg.h5
    !! Assumes xg_file_id and the group ids are already open via hdf5_open_xg_file()
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1
    real(kind=8), intent(in) :: data(n1)
    integer(HID_T) :: loc_id

    ! Validate that file is open
    if (.not. xg_file_open) then
      write(*,*) "ERROR: xg.h5 file not open. Call hdf5_open_xg_file() first"
      return
    endif

    select case (trim(group))
    case ("hydro")
      loc_id = xg_hydro_group_id
    case ("M1")
      loc_id = xg_M1_group_id
    case default
      write(*,*) "ERROR: unknown xg.h5 group '", trim(group), "' for ", trim(varname)
      stop "hdf5_append_grid_data_1d: group must be 'hydro' or 'M1'"
    end select

    if (loc_id < 0) then
      write(*,*) "ERROR: xg.h5 group '", trim(group), "' not open. Call hdf5_open_xg_file() first"
      return
    endif

    call hdf5_append_row(loc_id, varname, data, n1, 1_HSIZE_T)

  end subroutine hdf5_append_grid_data_1d

  subroutine hdf5_append_xg_scalar(varname, value)
    !! Append a scalar (e.g. time) to the root of xg.h5 - one shared time
    !! axis for both /hydro and /M1
    !! Assumes xg_file_id is already open via hdf5_open_xg_file()
    character(len=*), intent(in) :: varname
    real(kind=8), intent(in) :: value

    ! Validate that file is open
    if (.not. xg_file_open) then
      write(*,*) "ERROR: xg.h5 file not open. Call hdf5_open_xg_file() first"
      return
    endif

    call hdf5_append_scalar_at(xg_file_id, varname, value, 1024_HSIZE_T)

  end subroutine hdf5_append_xg_scalar

#ifdef HAVE_BURN
  subroutine hdf5_append_burn_scalar(varname, value)
    !! Append a scalar (e.g. time) to /fields in burn.h5
    !! Assumes burn_file_id and burn_group_id are already open via hdf5_open_burn_file()
    character(len=*), intent(in) :: varname
    real(kind=8), intent(in) :: value

    if (.not. burn_file_open .or. burn_group_id < 0) then
      write(*,*) "ERROR: burn.h5 file/group not open. Call hdf5_open_burn_file() first"
      return
    endif

    call hdf5_append_scalar_at(burn_group_id, varname, value, 1024_HSIZE_T)

  end subroutine hdf5_append_burn_scalar

  subroutine hdf5_append_burn_1d(varname, data, n)
    !! Append one row (one dump) of a grid variable to /fields in burn.h5
    !! Assumes burn_file_id and burn_group_id are already open via hdf5_open_burn_file()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n
    real(kind=8), intent(in) :: data(n)

    if (.not. burn_file_open .or. burn_group_id < 0) then
      write(*,*) "ERROR: burn.h5 file/group not open. Call hdf5_open_burn_file() first"
      return
    endif

    call hdf5_append_row(burn_group_id, varname, data, n, 1_HSIZE_T)

  end subroutine hdf5_append_burn_1d

  subroutine hdf5_append_burn_2d(varname, values, na, nb)
    !! Append one (na,nb) slab (one dump) to /fields in burn.h5 (e.g. Yion)
    !! Assumes burn_file_id and burn_group_id are already open via hdf5_open_burn_file()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: na, nb
    real(kind=8), intent(in) :: values(na,nb)

    if (.not. burn_file_open .or. burn_group_id < 0) then
      write(*,*) "ERROR: burn.h5 file/group not open. Call hdf5_open_burn_file() first"
      return
    endif

    call hdf5_append_row_2d(burn_group_id, varname, values, na, nb, 1_HSIZE_T)

  end subroutine hdf5_append_burn_2d
#endif

  subroutine hdf5_append_scalar(varname, value)
    !! Append scalar value to time series in dat.h5 using persistent file handles
    !! Assumes dat_file_id and dat_scalars_group_id are already open via hdf5_open_dat_file()
    character(len=*), intent(in) :: varname
    real(kind=8), intent(in) :: value

    ! Validate that file is open
    if (.not. dat_file_open .or. dat_scalars_group_id < 0) then
      write(*,*) "ERROR: dat.h5 file/group not open. Call hdf5_open_dat_file() first"
      return
    endif

    call hdf5_append_scalar_at(dat_scalars_group_id, varname, value, 100_HSIZE_T)

  end subroutine hdf5_append_scalar

  subroutine hdf5_append_scalar_at(loc_id, varname, value, chunk_len)
    !! Append scalar value to an extendible 1D dataset under loc_id,
    !! creating the dataset (chunked, unlimited) on first call
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: varname
    real(kind=8), intent(in) :: value
    integer(HSIZE_T), intent(in) :: chunk_len
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(1), maxdims(1), offset(1), count(1)
    real(kind=8) :: data_buf(1)
    logical :: exists

    ! Check if dataset exists
    call h5lexists_f(loc_id, trim(varname), exists, error)

    if (.not. exists) then
      ! Create new resizable dataset
      dims(1) = 1
      maxdims(1) = H5S_UNLIMITED_F

      call h5screate_simple_f(1, dims, dspace_id, error, maxdims)

      ! Create dataset creation property list with chunking (required for resizable datasets)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
      call h5pset_chunk_f(dcpl_id, 1, [chunk_len], error)

      call h5dcreate_f(loc_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, error, dcpl_id)

      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)

      ! Write first value
      data_buf(1) = value
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_buf, [1_HSIZE_T], error)

      call h5dclose_f(dset_id, error)

    else
      ! Append to existing dataset
      call h5dopen_f(loc_id, trim(varname), dset_id, error)

      ! Get current size
      call h5dget_space_f(dset_id, dspace_id, error)
      call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)

      ! Resize dataset
      dims(1) = dims(1) + 1
      call h5dset_extent_f(dset_id, dims, error)

      ! Get updated dataspace
      call h5dget_space_f(dset_id, dspace_id, error)

      ! Create memory space for single value
      call h5screate_simple_f(1, [1_HSIZE_T], memspace_id, error)

      ! Set offset to end of dataset
      offset(1) = dims(1) - 1
      count(1) = 1
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error)

      ! Write value
      data_buf(1) = value
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_buf, [1_HSIZE_T], error, memspace_id, dspace_id)

      call h5sclose_f(memspace_id, error)
      call h5sclose_f(dspace_id, error)
      call h5dclose_f(dset_id, error)
    endif

  end subroutine hdf5_append_scalar_at

  subroutine hdf5_append_scalar_array(varname, values, n)
    !! Append array of scalars to 2D time series in dat.h5 using persistent file handles
    !! Assumes dat_file_id and dat_scalars_group_id are already open via hdf5_open_dat_file()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n
    real(kind=8), intent(in) :: values(n)

    ! Validate that file is open
    if (.not. dat_file_open .or. dat_scalars_group_id < 0) then
      write(*,*) "ERROR: dat.h5 file/group not open. Call hdf5_open_dat_file() first"
      return
    endif

    call hdf5_append_row(dat_scalars_group_id, varname, values, n, 100_HSIZE_T)

  end subroutine hdf5_append_scalar_array

  subroutine hdf5_append_row(loc_id, varname, values, n, chunk_rows)
    !! Append one row to an extendible 2D dataset (ntime_unlimited, n) under loc_id,
    !! creating the dataset (chunked, unlimited along time) on first call
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n
    real(kind=8), intent(in) :: values(n)
    integer(HSIZE_T), intent(in) :: chunk_rows
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(2), maxdims(2), offset(2), count(2)
    logical :: exists

    ! Check if dataset exists
    call h5lexists_f(loc_id, trim(varname), exists, error)

    if (.not. exists) then
      ! Create new resizable 2D dataset
      dims(1) = 1
      dims(2) = n
      maxdims(1) = H5S_UNLIMITED_F
      maxdims(2) = n

      call h5screate_simple_f(2, dims, dspace_id, error, maxdims)

      ! Create dataset creation property list with chunking (required for resizable datasets)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
      call h5pset_chunk_f(dcpl_id, 2, [chunk_rows, int(n, HSIZE_T)], error)

      call h5dcreate_f(loc_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, error, dcpl_id)

      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)

      ! Write first row
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, reshape(values, [1, n]), [1_HSIZE_T, int(n, HSIZE_T)], error)

      call h5dclose_f(dset_id, error)

    else
      ! Append to existing dataset
      call h5dopen_f(loc_id, trim(varname), dset_id, error)

      ! Get current size
      call h5dget_space_f(dset_id, dspace_id, error)
      call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)

      ! Resize dataset
      dims(1) = dims(1) + 1
      call h5dset_extent_f(dset_id, dims, error)

      ! Get updated dataspace
      call h5dget_space_f(dset_id, dspace_id, error)

      ! Create memory space for single row
      call h5screate_simple_f(2, [1_HSIZE_T, int(n, HSIZE_T)], memspace_id, error)

      ! Set offset to end of dataset
      offset(1) = dims(1) - 1
      offset(2) = 0
      count(1) = 1
      count(2) = n
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error)

      ! Write row
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, reshape(values, [1, n]), [1_HSIZE_T, int(n, HSIZE_T)], &
                      error, memspace_id, dspace_id)

      call h5sclose_f(memspace_id, error)
      call h5sclose_f(dspace_id, error)
      call h5dclose_f(dset_id, error)
    endif

  end subroutine hdf5_append_row

  subroutine hdf5_append_row_2d(loc_id, varname, values, na, nb, chunk_rows)
    !! Append one (na,nb) slab to an extendible 3D dataset (ntime_unlimited, na, nb)
    !! under loc_id, creating the dataset (chunked, unlimited along time) on first call
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: varname
    integer, intent(in) :: na, nb
    real(kind=8), intent(in) :: values(na,nb)
    integer(HSIZE_T), intent(in) :: chunk_rows
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(3), maxdims(3), offset(3), count(3)
    logical :: exists

    ! Check if dataset exists
    call h5lexists_f(loc_id, trim(varname), exists, error)

    if (.not. exists) then
      ! Create new resizable 3D dataset
      dims(1) = 1
      dims(2) = na
      dims(3) = nb
      maxdims(1) = H5S_UNLIMITED_F
      maxdims(2) = na
      maxdims(3) = nb

      call h5screate_simple_f(3, dims, dspace_id, error, maxdims)

      ! Create dataset creation property list with chunking (required for resizable datasets)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
      call h5pset_chunk_f(dcpl_id, 3, [chunk_rows, int(na, HSIZE_T), int(nb, HSIZE_T)], error)

      call h5dcreate_f(loc_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, error, dcpl_id)

      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)

      ! Write first slab
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, reshape(values, [1, na, nb]), &
                      [1_HSIZE_T, int(na, HSIZE_T), int(nb, HSIZE_T)], error)

      call h5dclose_f(dset_id, error)

    else
      ! Append to existing dataset
      call h5dopen_f(loc_id, trim(varname), dset_id, error)

      ! Get current size
      call h5dget_space_f(dset_id, dspace_id, error)
      call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)

      ! Resize dataset
      dims(1) = dims(1) + 1
      call h5dset_extent_f(dset_id, dims, error)

      ! Get updated dataspace
      call h5dget_space_f(dset_id, dspace_id, error)

      ! Create memory space for single slab
      call h5screate_simple_f(3, [1_HSIZE_T, int(na, HSIZE_T), int(nb, HSIZE_T)], memspace_id, error)

      ! Set offset to end of dataset
      offset(1) = dims(1) - 1
      offset(2) = 0
      offset(3) = 0
      count(1) = 1
      count(2) = na
      count(3) = nb
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error)

      ! Write slab
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, reshape(values, [1, na, nb]), &
                      [1_HSIZE_T, int(na, HSIZE_T), int(nb, HSIZE_T)], &
                      error, memspace_id, dspace_id)

      call h5sclose_f(memspace_id, error)
      call h5sclose_f(dspace_id, error)
      call h5dclose_f(dset_id, error)
    endif

  end subroutine hdf5_append_row_2d

  subroutine hdf5_write_metadata(file_path, n1, ghosts1, eoskey, do_M1, do_rotation, &
                                 do_turbulence, GR, &
                                 number_groups, number_species)

    !! Write simulation metadata to dat.h5 and xg.h5

    character(len=256), intent(in) :: file_path
    integer, intent(in) :: n1, ghosts1, eoskey, number_groups, number_species
    logical, intent(in) :: do_M1, do_rotation, do_turbulence, GR
    integer :: error
    integer(HID_T) :: file_id, metadata_group_id, aid, asid
    integer :: int_val

    ! Open file
    call h5fopen_f(trim(file_path), H5F_ACC_RDWR_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to open ", file_path ," for metadata writing"
      return
    endif

    ! Open /metadata group
    call h5gopen_f(file_id, "/metadata", metadata_group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to open /metadata group"
      call h5fclose_f(file_id, error)
      return
    endif

    ! Create scalar attribute space
    call h5screate_f(H5S_SCALAR_F, asid, error)

    ! Write integer attributes
    call h5acreate_f(metadata_group_id, "n1", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, n1, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5acreate_f(metadata_group_id, "eoskey", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, eoskey, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5acreate_f(metadata_group_id, "number_groups", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, number_groups, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5acreate_f(metadata_group_id, "number_species", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, number_species, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    ! Write logical attributes (as integers: 0=false, 1=true)
    int_val = merge(1, 0, do_M1)
    call h5acreate_f(metadata_group_id, "do_M1", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, int_val, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    int_val = merge(1, 0, do_rotation)
    call h5acreate_f(metadata_group_id, "do_rotation", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, int_val, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    int_val = merge(1, 0, do_turbulence)
    call h5acreate_f(metadata_group_id, "do_turbulence", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, int_val, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    int_val = merge(1, 0, GR)
    call h5acreate_f(metadata_group_id, "GR", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, int_val, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5sclose_f(asid, error)
    call h5gclose_f(metadata_group_id, error)
    call h5fclose_f(file_id, error)

  end subroutine hdf5_write_metadata

#ifdef HAVE_BURN
  subroutine hdf5_write_burn_metadata(n1, ghosts1, nspec, nspec_net, &
                                      track_free_nucleons, T_eos_high, T_eos_low)

    !! Write nuclear-burning metadata to burn.h5 /metadata

    integer, intent(in) :: n1, ghosts1, nspec, nspec_net
    logical, intent(in) :: track_free_nucleons
    real(kind=8), intent(in) :: T_eos_high, T_eos_low
    integer :: error
    integer(HID_T) :: file_id, metadata_group_id, aid, asid
    integer :: int_val

    ! Open file
    call h5fopen_f(trim(burn_file_path), H5F_ACC_RDWR_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to open burn.h5 for metadata writing"
      return
    endif

    ! Open /metadata group
    call h5gopen_f(file_id, "/metadata", metadata_group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to open /metadata group in burn.h5"
      call h5fclose_f(file_id, error)
      return
    endif

    ! Create scalar attribute space
    call h5screate_f(H5S_SCALAR_F, asid, error)

    ! Write integer attributes
    call h5acreate_f(metadata_group_id, "n1", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, n1, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5acreate_f(metadata_group_id, "ghosts1", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, ghosts1, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5acreate_f(metadata_group_id, "nspec", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, nspec, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5acreate_f(metadata_group_id, "nspec_net", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, nspec_net, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    ! Logical attribute (as integer: 0=false, 1=true)
    int_val = merge(1, 0, track_free_nucleons)
    call h5acreate_f(metadata_group_id, "track_free_nucleons", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, int_val, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    ! Write double attributes (regime thresholds, Kelvin)
    call h5acreate_f(metadata_group_id, "T_eos_high", H5T_NATIVE_DOUBLE, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_DOUBLE, T_eos_high, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5acreate_f(metadata_group_id, "T_eos_low", H5T_NATIVE_DOUBLE, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_DOUBLE, T_eos_low, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)

    call h5sclose_f(asid, error)
    call h5gclose_f(metadata_group_id, error)
    call h5fclose_f(file_id, error)

  end subroutine hdf5_write_burn_metadata
#endif

  subroutine hdf5_finalize()

    !! Close HDF5 library
    integer :: error
    call h5close_f(error)

    if (error /= 0) then
      write(*,*) "ERROR: Failed to initialize Fortran HDF5 library!"
      stop
    endif

  end subroutine hdf5_finalize

  subroutine hdf5_output_init()

   use GR1D_module
#ifdef HAVE_BURN
   use composition, only: nspec, nspec_net
#endif
    ! Initialize HDF5 library
    call hdf5_initialize()

    ! Create output files
    call hdf5_create_output_files(outdir)

    ! Write metadata to dat.h5
    call hdf5_write_metadata(dat_file_path, n1, ghosts1, eoskey, do_M1, &
                             do_rotation, do_turbulence, GR, &
                             number_groups, number_species)
    call hdf5_write_metadata(xg_file_path, n1, ghosts1, eoskey, do_M1, &
                             do_rotation, do_turbulence, GR, &
                             number_groups, number_species)

#ifdef HAVE_BURN
    call hdf5_write_burn_metadata(n1, ghosts1, nspec, nspec_net, &
                                  track_free_nucleons, T_eos_high, T_eos_low)
#endif

    write(*,*) "HDF5 output initialized successfully"
    write(*,*) "  Grid output: ", trim(adjustl(xg_file_path))
    write(*,*) "  Scalar output: ", trim(adjustl(dat_file_path))
#ifdef HAVE_BURN
    write(*,*) "  Burn output: ", trim(adjustl(burn_file_path))
#endif

    call hdf5_finalize()

  end subroutine hdf5_output_init

end module hdf5_output_utils
