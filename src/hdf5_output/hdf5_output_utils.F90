!-*-f90-*-
module hdf5_output_utils
  !! HDF5 Output Utilities Module
  !! Provides wrapper functions for writing GR1D_STIR output to HDF5 files
  !! Two files are created: xg.h5 (grid profiles) and dat.h5 (scalar time series)
  !! Files are kept open within output scope for efficiency
  
  use hdf5
  implicit none
  
  private
  
  ! Public subroutines
  public :: hdf5_initialize
  public :: hdf5_create_output_files
  public :: hdf5_write_root_dataset_1d
  public :: hdf5_write_grid_data_1d
  public :: hdf5_write_grid_data_2d
  public :: hdf5_append_scalar
  public :: hdf5_append_scalar_array
  public :: hdf5_finalize
  public :: hdf5_write_metadata
  public :: hdf5_increment_output_counter
  public :: hdf5_output_init
  public :: hdf5_open_xg_for_output
  public :: hdf5_close_xg_file
  public :: hdf5_open_dat_for_scalars
  public :: hdf5_close_dat_file
  
  ! Module variables
  character(len=1024), save :: xg_file_path = ""
  character(len=1024), save :: dat_file_path = ""
  integer, save :: output_counter = 0
  
  ! File/group handle storage for persistent open files
  integer(HID_T), save :: xg_file_id = -1
  integer(HID_T), save :: xg_group_id = -1
  integer(HID_T), save :: dat_file_id = -1
  integer(HID_T), save :: dat_scalars_group_id = -1
  logical, save :: xg_file_open = .false.
  logical, save :: dat_file_open = .false.
  
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
    !! Create xg.h5 and dat.h5 files
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
    call h5pclose_f(fcpl_id, error)
    
  end subroutine hdf5_create_output_files

  subroutine hdf5_open_xg_for_output()
    !! Open xg.h5 file and create/open output_#### group for persistent writing
    integer :: error
    character(len=256) :: group_name
    logical :: group_exists
    integer(HID_T) :: aid, asid
    
    ! Open file if not already open
    if (.not. xg_file_open) then
      call h5fopen_f(trim(xg_file_path), H5F_ACC_RDWR_F, xg_file_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to open xg.h5 for writing"
        return
      endif
      xg_file_open = .true.
    endif
    
    ! Create group name for this output timestep
    write(group_name, '(A,I4.4)') "/output_", output_counter
    
    ! Check if group exists
    call h5lexists_f(xg_file_id, trim(group_name), group_exists, error)
    
    if (.not. group_exists) then
      ! Create new group
      call h5gcreate_f(xg_file_id, trim(group_name), xg_group_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to create group ", trim(group_name)
        return
      endif
      
      ! Write attributes
      call h5screate_f(H5S_SCALAR_F, asid, error)
      
      ! Write time attribute
      call h5acreate_f(xg_group_id, "time", H5T_NATIVE_DOUBLE, asid, aid, error)
      call h5awrite_f(aid, H5T_NATIVE_DOUBLE, [0.0d0], [1_HSIZE_T], error)
      call h5aclose_f(aid, error)
      
      ! Write timestep attribute
      call h5acreate_f(xg_group_id, "timestep", H5T_NATIVE_INTEGER, asid, aid, error)
      call h5awrite_f(aid, H5T_NATIVE_INTEGER, [0], [1_HSIZE_T], error)
      call h5aclose_f(aid, error)
      
      ! Write output_number attribute
      call h5acreate_f(xg_group_id, "output_number", H5T_NATIVE_INTEGER, asid, aid, error)
      call h5awrite_f(aid, H5T_NATIVE_INTEGER, output_counter, [1_HSIZE_T], error)
      call h5aclose_f(aid, error)
      
      call h5sclose_f(asid, error)
    else
      ! Open existing group
      call h5gopen_f(xg_file_id, trim(group_name), xg_group_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to open group ", trim(group_name)
        return
      endif
    endif
    
  end subroutine hdf5_open_xg_for_output

  subroutine hdf5_close_xg_file()
    !! Close xg.h5 file and group
    integer :: error
    
    if (xg_file_open) then
      if (xg_group_id >= 0) then
        call h5gclose_f(xg_group_id, error)
        xg_group_id = -1
      endif
      call h5fclose_f(xg_file_id, error)
      xg_file_id = -1
      xg_file_open = .false.
    endif
    
  end subroutine hdf5_close_xg_file

  subroutine hdf5_open_dat_for_scalars()
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
    
  end subroutine hdf5_open_dat_for_scalars

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

  subroutine hdf5_write_root_dataset_1d(varname, data, n1)
    !! Write 1D grid variable directly to the root of the HDF5 file
    !! Opens file, creates dataset, writes data, closes file
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1
    real(kind=8), intent(in) :: data(n1)
    
    integer :: error
    integer(HID_T) :: file_id, dset_id, dspace_id, dcpl_id
    integer(HSIZE_T) :: dims(1)
    
    ! Open file
    call h5fopen_f(trim(xg_file_path), H5F_ACC_RDWR_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to open xg.h5 for writing"
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
    
  end subroutine hdf5_write_root_dataset_1d

  subroutine hdf5_write_grid_data_1d(varname, data, n1)
    !! Write 1D grid variable to xg.h5 using persistent file handles
    !! Assumes xg_file_id and xg_group_id are already open via hdf5_open_xg_for_output()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1
    real(kind=8), intent(in) :: data(n1)
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id
    integer(HSIZE_T) :: dims(1)
    
    ! Validate that file is open
    if (.not. xg_file_open .or. xg_group_id < 0) then
      write(*,*) "ERROR: xg.h5 file/group not open. Call hdf5_open_xg_for_output() first"
      return
    endif
    
    dims(1) = n1
    
    ! Create dataspace
    call h5screate_simple_f(1, dims, dspace_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create dataspace for ", trim(varname)
      return
    endif
    
    ! Create dataset creation property list (no compression)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
    
    ! Create dataset in the group
    call h5dcreate_f(xg_group_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, error, dcpl_id)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create dataset for ", trim(varname)
      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)
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
    
  end subroutine hdf5_write_grid_data_1d

  subroutine hdf5_write_grid_data_2d(varname, data, n1, n2)
    !! Write 2D grid variable (e.g., spectra) to xg.h5 using persistent file handles
    !! Assumes xg_file_id and xg_group_id are already open via hdf5_open_xg_for_output()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1, n2
    real(kind=8), intent(in) :: data(n1, n2)
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id
    integer(HSIZE_T) :: dims(2)
    
    ! Validate that file is open
    if (.not. xg_file_open .or. xg_group_id < 0) then
      write(*,*) "ERROR: xg.h5 file/group not open. Call hdf5_open_xg_for_output() first"
      return
    endif
    
    dims(1) = n1
    dims(2) = n2
    
    ! Create dataspace
    call h5screate_simple_f(2, dims, dspace_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create dataspace for ", trim(varname)
      return
    endif
    
    ! Create dataset creation property list (no compression)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)

    ! Create dataset in the group
    call h5dcreate_f(xg_group_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, error, dcpl_id)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create dataset for ", trim(varname)
      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)
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
    
  end subroutine hdf5_write_grid_data_2d

  subroutine hdf5_append_scalar(varname, value)
    !! Append scalar value to time series in dat.h5 using persistent file handles
    !! Assumes dat_file_id and dat_scalars_group_id are already open via hdf5_open_dat_for_scalars()
    character(len=*), intent(in) :: varname
    real(kind=8), intent(in) :: value
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(1), maxdims(1), offset(1), count(1)
    real(kind=8) :: data_buf(1)
    logical :: exists
    
    ! Validate that file is open
    if (.not. dat_file_open .or. dat_scalars_group_id < 0) then
      write(*,*) "ERROR: dat.h5 file/group not open. Call hdf5_open_dat_for_scalars() first"
      return
    endif
    
    ! Check if dataset exists
    call h5lexists_f(dat_scalars_group_id, trim(varname), exists, error)
    
    if (.not. exists) then
      ! Create new resizable dataset
      dims(1) = 1
      maxdims(1) = H5S_UNLIMITED_F
      
      call h5screate_simple_f(1, dims, dspace_id, error, maxdims)
      
      ! Create dataset creation property list with chunking (required for resizable datasets)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
      call h5pset_chunk_f(dcpl_id, 1, [100_HSIZE_T], error)
      
      call h5dcreate_f(dat_scalars_group_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, error, dcpl_id)
      
      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)
      
      ! Write first value
      data_buf(1) = value
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_buf, [1_HSIZE_T], error)
      
      call h5dclose_f(dset_id, error)
      
    else
      ! Append to existing dataset
      call h5dopen_f(dat_scalars_group_id, trim(varname), dset_id, error)
      
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
    
  end subroutine hdf5_append_scalar

  subroutine hdf5_append_scalar_array(varname, values, n)
    !! Append array of scalars to 2D time series in dat.h5 using persistent file handles
    !! Assumes dat_file_id and dat_scalars_group_id are already open via hdf5_open_dat_for_scalars()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n
    real(kind=8), intent(in) :: values(n)
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(2), maxdims(2), offset(2), count(2)
    logical :: exists
    
    ! Validate that file is open
    if (.not. dat_file_open .or. dat_scalars_group_id < 0) then
      write(*,*) "ERROR: dat.h5 file/group not open. Call hdf5_open_dat_for_scalars() first"
      return
    endif
    
    ! Check if dataset exists
    call h5lexists_f(dat_scalars_group_id, trim(varname), exists, error)
    
    if (.not. exists) then
      ! Create new resizable 2D dataset
      dims(1) = 1
      dims(2) = n
      maxdims(1) = H5S_UNLIMITED_F
      maxdims(2) = n
      
      call h5screate_simple_f(2, dims, dspace_id, error, maxdims)
      
      ! Create dataset creation property list with chunking (required for resizable datasets)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
      call h5pset_chunk_f(dcpl_id, 2, [100_HSIZE_T, int(n, HSIZE_T)], error)
      
      call h5dcreate_f(dat_scalars_group_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, error, dcpl_id)
      
      call h5sclose_f(dspace_id, error)
      call h5pclose_f(dcpl_id, error)
      
      ! Write first row
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, reshape(values, [1, n]), [1_HSIZE_T, int(n, HSIZE_T)], error)
      
      call h5dclose_f(dset_id, error)
      
    else
      ! Append to existing dataset
      call h5dopen_f(dat_scalars_group_id, trim(varname), dset_id, error)
      
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
    
  end subroutine hdf5_append_scalar_array

  subroutine hdf5_write_metadata(n1, eoskey, do_M1, do_rotation, do_turbulence, GR, &
                                  number_groups, number_species)
    !! Write simulation metadata to dat.h5
    integer, intent(in) :: n1, eoskey, number_groups, number_species
    logical, intent(in) :: do_M1, do_rotation, do_turbulence, GR
    integer :: error
    integer(HID_T) :: file_id, metadata_group_id, aid, asid
    integer :: int_val
    
    ! Open file
    call h5fopen_f(trim(dat_file_path), H5F_ACC_RDWR_F, file_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to open dat.h5 for metadata writing"
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

  subroutine hdf5_increment_output_counter()

    !! Increment output counter (called after each output)
    output_counter = output_counter + 1

  end subroutine hdf5_increment_output_counter

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
    ! Initialize HDF5 library
    call hdf5_initialize()
    
    ! Create output files
    call hdf5_create_output_files(outdir)
    
    ! Write metadata to dat.h5
    call hdf5_write_metadata(n1, eoskey, do_M1, do_rotation, do_turbulence, GR, &
                             number_groups, number_species)
    
    write(*,*) "HDF5 output initialized successfully"
    write(*,*) "  Grid output: ", trim(adjustl(outdir))//"/xg.h5"
    write(*,*) "  Scalar output: ", trim(adjustl(outdir))//"/dat.h5"
    
    call hdf5_finalize()

  end subroutine hdf5_output_init

end module hdf5_output_utils
