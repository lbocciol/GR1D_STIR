!-*-f90-*-
module hdf5_output_utils
  !! HDF5 Output Utilities Module
  !! Provides wrapper functions for writing GR1D_STIR output to HDF5 files
  !! Two files are created: xg.h5 (grid profiles) and dat.h5 (scalar time series)
  !! Files are kept open within output scope for efficiency
  !!
  !! xg.h5 layout 2 (/metadata attribute xg_layout = 2): every grid variable is one dataset at
  !! the file root that grows by one output at a time, shape (n, nT) as seen from python/C
  !! (n = n1 or number_groups), plus /time of shape (nT). One output per chunk, so each output
  !! only writes new chunks. /time is written last in every output: readers take nT = len(time)
  !! and ignore any extra column left by a crash in the middle of an output.
  !! (Layout 1, before: one group /output_NNNNN per output with one dataset per variable.)
  
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
  public :: hdf5_open_xg_file
  public :: hdf5_close_xg_file
  public :: hdf5_open_dat_file
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

  ! Layout of xg.h5 written by this version (stored as /metadata attribute xg_layout)
  integer, parameter :: XG_LAYOUT = 2
  
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
    
        ! Create /metadata group
    call h5gcreate_f(file_id, "/metadata", group_id, error)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to create /metadata group"
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

  subroutine hdf5_open_xg_file()
    !! Open xg.h5 file for persistent writing of one output (time series datasets at the root)
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
    
  end subroutine hdf5_open_xg_file

  subroutine hdf5_close_xg_file()
    !! Close xg.h5 file and group
    integer :: error
    
    if (xg_file_open) then
      call h5fclose_f(xg_file_id, error)
      xg_file_id = -1
      xg_file_open = .false.
    endif
    
  end subroutine hdf5_close_xg_file

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
    !! Append one output of a grid variable to its time series /<varname> in xg.h5
    !! (shape (n1, nT) from python/C, 1D (nT) if n1 = 1, e.g. time), creating it at the first
    !! output. Assumes xg.h5 is open via hdf5_open_xg_file()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1
    real(kind=8), intent(in) :: data(n1)
    integer :: error, rank
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(2), maxdims(2), offset(2), count(2), chunk(2)
    logical :: exists
    
    ! Validate that file is open
    if (.not. xg_file_open) then
      write(*,*) "ERROR: xg.h5 file not open. Call hdf5_open_xg_file() first"
      return
    endif
    
    ! Fortran dimension order: (time, n1); a single value per output is stored as 1D (time)
    rank = 2
    if (n1 == 1) rank = 1
    count(1) = 1
    count(2) = n1
    
    call h5lexists_f(xg_file_id, trim(varname), exists, error)
    
    if (.not. exists) then
      ! Create resizable dataset with one output, chunked one output at a time
      dims = count
      maxdims(1) = H5S_UNLIMITED_F
      maxdims(2) = n1
      chunk = count
      call h5screate_simple_f(rank, dims(1:rank), dspace_id, error, maxdims(1:rank))
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
      call h5pset_chunk_f(dcpl_id, rank, chunk(1:rank), error)
      call h5dcreate_f(xg_file_id, trim(varname), H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, error, dcpl_id)
      call h5pclose_f(dcpl_id, error)
      if (error /= 0) then
        write(*,*) "ERROR: Failed to create dataset for ", trim(varname)
        call h5sclose_f(dspace_id, error)
        return
      endif
      call h5sclose_f(dspace_id, error)
      offset = 0
    else
      ! Extend the existing dataset by one output
      call h5dopen_f(xg_file_id, trim(varname), dset_id, error)
      call h5dget_space_f(dset_id, dspace_id, error)
      call h5sget_simple_extent_dims_f(dspace_id, dims(1:rank), maxdims(1:rank), error)
      call h5sclose_f(dspace_id, error)
      offset = 0
      offset(1) = dims(1)
      dims(1) = dims(1) + 1
      call h5dset_extent_f(dset_id, dims(1:rank), error)
    endif
    
    ! Write the new output into the last position along time
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset(1:rank), count(1:rank), error)
    call h5screate_simple_f(rank, count(1:rank), memspace_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count(1:rank), error, memspace_id, dspace_id)
    if (error /= 0) then
      write(*,*) "ERROR: Failed to write data for ", trim(varname)
    endif
    
    call h5sclose_f(memspace_id, error)
    call h5sclose_f(dspace_id, error)
    call h5dclose_f(dset_id, error)
    
  end subroutine hdf5_write_grid_data_1d

  subroutine hdf5_write_grid_data_2d(varname, data, n1, n2)
    !! Write 2D grid variable (e.g., spectra) to xg.h5 using persistent file handles
    !! Assumes xg_file_id and xg_group_id are already open via hdf5_open_xg_file()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n1, n2
    real(kind=8), intent(in) :: data(n1, n2)
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id
    integer(HSIZE_T) :: dims(2)
    
    ! Validate that file is open
    if (.not. xg_file_open .or. xg_group_id < 0) then
      write(*,*) "ERROR: xg.h5 file/group not open. Call hdf5_open_xg_file() first"
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
    !! Assumes dat_file_id and dat_scalars_group_id are already open via hdf5_open_dat_file()
    character(len=*), intent(in) :: varname
    real(kind=8), intent(in) :: value
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(1), maxdims(1), offset(1), count(1)
    real(kind=8) :: data_buf(1)
    logical :: exists
    
    ! Validate that file is open
    if (.not. dat_file_open .or. dat_scalars_group_id < 0) then
      write(*,*) "ERROR: dat.h5 file/group not open. Call hdf5_open_dat_file() first"
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
    !! Assumes dat_file_id and dat_scalars_group_id are already open via hdf5_open_dat_file()
    character(len=*), intent(in) :: varname
    integer, intent(in) :: n
    real(kind=8), intent(in) :: values(n)
    integer :: error
    integer(HID_T) :: dset_id, dspace_id, dcpl_id, memspace_id
    integer(HSIZE_T) :: dims(2), maxdims(2), offset(2), count(2)
    logical :: exists
    
    ! Validate that file is open
    if (.not. dat_file_open .or. dat_scalars_group_id < 0) then
      write(*,*) "ERROR: dat.h5 file/group not open. Call hdf5_open_dat_file() first"
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

  subroutine hdf5_write_xg_layout()
    !! Store the xg.h5 layout version as /metadata attribute xg_layout (xg.h5 only)
    integer :: error
    integer(HID_T) :: file_id, metadata_group_id, aid, asid

    call h5fopen_f(trim(xg_file_path), H5F_ACC_RDWR_F, file_id, error)
    call h5gopen_f(file_id, "/metadata", metadata_group_id, error)
    call h5screate_f(H5S_SCALAR_F, asid, error)
    call h5acreate_f(metadata_group_id, "xg_layout", H5T_NATIVE_INTEGER, asid, aid, error)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, XG_LAYOUT, [1_HSIZE_T], error)
    call h5aclose_f(aid, error)
    call h5sclose_f(asid, error)
    call h5gclose_f(metadata_group_id, error)
    call h5fclose_f(file_id, error)

  end subroutine hdf5_write_xg_layout

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
    call hdf5_write_metadata(dat_file_path, n1, ghosts1, eoskey, do_M1, &
                             do_rotation, do_turbulence, GR, &
                             number_groups, number_species)
    call hdf5_write_metadata(xg_file_path, n1, ghosts1, eoskey, do_M1, &
                             do_rotation, do_turbulence, GR, &
                             number_groups, number_species)
    call hdf5_write_xg_layout()

    write(*,*) "HDF5 output initialized successfully"
    write(*,*) "  Grid output: ", trim(adjustl(dat_file_path))
    write(*,*) "  Scalar output: ", trim(adjustl(xg_file_path))
    
    call hdf5_finalize()

  end subroutine hdf5_output_init

end module hdf5_output_utils
