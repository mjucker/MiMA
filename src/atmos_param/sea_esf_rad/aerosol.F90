                    module aerosol_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code to initialize/allocate aerosol climatology
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes prescribed aerosol climatology from input file,
!  allocates necessary memory space and interpolate the aerosol climatology
!  to the model specification. Afterwards the memory space is deallocated,
!  the aerosol climatology information is freed.
! </DESCRIPTION>
!  shared modules:

use time_manager_mod,  only: time_type, time_manager_init, operator(+),&
                             set_date, operator(-), print_date, &
                             set_time, days_in_month, get_date, &
                             operator(>), operator(/=)
use diag_manager_mod,  only: diag_manager_init, get_base_time
use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, WARNING, close_file
use interpolator_mod,  only: interpolate_type, interpolator_init, &
                             interpolator, interpolator_end, &
                             CONSTANT, INTERP_WEIGHTED_P
use mpp_io_mod,        only: mpp_open, mpp_close, MPP_RDONLY,   &
                             MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI,  &
                             MPP_SINGLE, mpp_io_init
use constants_mod,     only: constants_init, RADIAN

!  shared radiation package modules:

use rad_utilities_mod, only: aerosol_type, rad_utilities_init

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    aersosol_mod defines the aerosol distribution that is to be seen
!    by the radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id: aerosol.F90,v 12.0 2005/04/14 15:44:04 fms Exp $'
character(len=128) :: tagname = '$Name: lima $'


!-----------------------------------------------------------------------
!------  interfaces -------

public          &
        aerosol_init, aerosol_driver, aerosol_end, &
        aerosol_dealloc

!private         &


!---------------------------------------------------------------------
!------  namelist ------

 character(len=32)      ::      &
         aerosol_data_source = 'climatology'
                                   ! source of aerosol data, either
                                   ! 'climatology' file (default) or
                                   ! single column 'input' file or
                                   ! calculate a column for location and
                                   ! time specified ('calculate_column')
integer, parameter     ::      &
        MAX_DATA_FIELDS = 100     ! maximum number of aerosol species
integer, parameter     ::      &
        MAX_AEROSOL_FAMILIES = 9 ! maximum number of aerosol families
character(len=64)      ::      &
        data_names(MAX_DATA_FIELDS) = '  ' 
                                  ! names of active aerosol species 
character(len=64)      ::      &
        filename = '  '           ! name of netcdf file containing 
                                  ! aerosol species to be activated
character(len=64)      ::      &
      family_names(MAX_AEROSOL_FAMILIES) = '  ' 
                                  ! names of active aerosol families 
logical, dimension(MAX_DATA_FIELDS) :: in_family1 = .false.
                                  ! aerosol n is in family 1 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family2 = .false.
                                  ! aerosol n is in family 2 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family3 = .false.
                                  ! aerosol n is in family 3 ?
logical,dimension(MAX_DATA_FIELDS) :: in_family4 = .false.
                                  ! aerosol n is in family 4 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family5 = .false.
                                  ! aerosol n is in family 5 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family6 = .false.
                                  ! aerosol n is in family 6 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family7 = .false.
                                  ! aerosol n is in family 7 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family8 = .false.
                                  ! aerosol n is in family 8 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family9 = .false.
                                  ! aerosol n is in family 9 ?
logical         :: use_aerosol_timeseries = .false.
                                  ! use a timeseries providing inter-
                                  ! annual aerosol variation ?
integer, dimension(6,MAX_DATA_FIELDS) :: aerosol_dataset_entry  =  1
                      ! time in aerosol data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
logical, dimension(MAX_AEROSOL_FAMILIES) ::   &
                                    volc_in_fam_col_opt_depth = .false.
                      ! is the volcanic contribution to column optical
                      ! depth to be included for this family in the
                      ! netcdf output fields ?
real,dimension(2)   :: lonb_col = (/-999., -999./)
                      ! longitudes defining the region to use for column
                      ! data calculation
real,dimension(2)   :: latb_col = (/-999., -999./)
                      ! latitudes defining the region to use for column
                      ! data calculation
integer, dimension(6)  :: time_col = (/0,0,0,0,0,0/)
                      ! time to use for column data calculation


namelist /aerosol_nml/                            &
                           aerosol_data_source,   &
                           lonb_col, latb_col, time_col, &
                           data_names, filename,  &
                           family_names,   &
                           use_aerosol_timeseries, &
                           aerosol_dataset_entry,  &               
                           in_family1, in_family2, in_family3, &
                           in_family4, in_family5, in_family6, &
                           in_family7, in_family8, in_family9, &
                           volc_in_fam_col_opt_depth
                           

!---------------------------------------------------------------------
!---- public data ----


!---------------------------------------------------------------------
!---- private data ----


!-------------------------------------------------------------------
!    specified_aerosol contains the column input aerosol concentration
!    ratio (kg/m**2).  used when aerosol_data_source = 'input'.
!-------------------------------------------------------------------
real, dimension (:), allocatable     ::  specified_aerosol
 
!---------------------------------------------------------------------
!   the following is an interpolate_type variable containing the
!   information about the aerosol species.
!---------------------------------------------------------------------
type(interpolate_type), save  :: Aerosol_interp

!--------------------------------------------------------------------
!    miscellaneous variables
!--------------------------------------------------------------------
logical  :: do_column_aerosol = .false.      ! using single column aero-
                                             ! sol data ?
integer  :: nfields=0                        ! number of active aerosol 
                                             ! species
integer  :: nfamilies=0                      ! number of active aerosol 
                                             ! families
logical  :: module_is_initialized = .false.  ! module has been 
                                             ! initialized  ?

type(time_type) :: Model_init_time  ! initial calendar time for model  
                                    ! [ time_type ]
type(time_type), dimension(:), allocatable ::   &
                   Aerosol_offset   ! difference between model initial
                                    ! time and aerosol timeseries app-
                                    ! lied at model initial time
                                    ! [ time_type ]
type(time_type), dimension(:), allocatable ::   &
                   Aerosol_entry    ! time in aerosol timeseries which
                                    ! is mapped to model initial time
                                    ! [ time_type ]
type(time_type)  ::   &
              Aerosol_column_time   ! time for which aerosol data is
                                    ! extracted from aerosol timeseries
                                    ! in 'calculate_columns' case
                                    ! [ time_type ]
logical, dimension(:), allocatable    ::     &
                   negative_offset 
                                    ! the model initial time is later 
                                    ! than the aerosol_dataset_entry 
                                    ! time  ?
integer , dimension(:), allocatable :: data_out_of_bounds, vert_interp
logical, dimension(:), allocatable :: using_fixed_year_data 
                                    ! are we using a fixed year
                                    !  of data from a timeseries file ?

#include "netcdf.inc"
 
!---------------------------------------------------------------------
!---------------------------------------------------------------------


                           contains


 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------
! <SUBROUTINE NAME="aerosol_init">
!  <OVERVIEW>
!   Subroutine to initialize/interpolate prescribed aerosol climatology
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize/interpolate prescribed aerosol climatology
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_init(lonb, latb, aerosol_names)
!  </TEMPLATE>
!  <IN NAME="lonb" TYPE="real">
!   Array of model longitudes on cell boundaries in [radians]
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   Array of model latitudes on cell boundaries in [radians]
!  </IN>
!  <IN NAME="aerosol_names" TYPE="character">
!   names of the activated aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_init (lonb, latb, aerosol_names,   &
                         aerosol_family_names)

!-----------------------------------------------------------------------
!    aerosol_init is the constructor for aerosol_mod.
!-----------------------------------------------------------------------

real, dimension(:),              intent(in)  :: lonb,latb
character(len=64), dimension(:), pointer     :: aerosol_names
character(len=64), dimension(:), pointer     :: aerosol_family_names

!----------------------------------------------------------------------
!  intent(in) variables:
!
!       lonb           array of model longitudes on cell boundaries 
!                      [ radians ]
!       latb           array of model latitudes at cell boundaries 
!                      [ radians ]
!
!   pointer variables:
!
!       aerosol_names  names of the activated aerosol species
!       aerosol_family_names  names of the activated aerosol families
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:
      
      integer   ::   unit, ierr, io       
      integer   ::   n

!---------------------------------------------------------------------
!    local variables:
!
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!         n          do-loop index
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call mpp_io_init
      call fms_init
      call diag_manager_init
      call rad_utilities_init
      call time_manager_init
      call constants_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosol_nml, iostat=io,  &
               end=10)
        ierr = check_nml_error(io,'aerosol_nml')
        end do
10      call close_file (unit)   
      endif                      
                                  
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                        write (stdlog(), nml=aerosol_nml)

!---------------------------------------------------------------------
!    when running standalone code on other than FMS level structure,
!    aerosol_data_source must be 'input'.
!---------------------------------------------------------------------
      if (trim(aerosol_data_source)== 'input') then
        do_column_aerosol = .true.
        call obtain_input_file_data
        nfields = 1
        nfamilies = 0
        allocate (aerosol_names(nfields))
        aerosol_names (1) = 'total_aerosol'
        allocate (aerosol_family_names(nfamilies))
      else

!-----------------------------------------------------------------------
!    count number of activated aerosol species. allocate a pointer 
!    array to return the names of the activated species to the calling 
!    routine.
!-----------------------------------------------------------------------
        do n=1,MAX_DATA_FIELDS
          if (data_names(n) /= ' '  ) then
            nfields = n
          else
            exit
          endif
        end do

        allocate (aerosol_names(nfields))
        aerosol_names (:) = data_names(1:nfields)

!----------------------------------------------------------------------
!    allocate module variables.
!----------------------------------------------------------------------
        allocate (Aerosol_offset(nfields), Aerosol_entry(nfields), &
                negative_offset(nfields), using_fixed_year_data(nfields)) 
        Aerosol_offset = set_time (0,0)
        Aerosol_entry = set_time (0,0)
        negative_offset = .false.
        using_fixed_year_data = .false.

!----------------------------------------------------------------------
!    define the model base time  (defined in diag_table) 
!----------------------------------------------------------------------
        Model_init_time = get_base_time()

!----------------------------------------------------------------------
!    process the aerosol timeseries file.  
!----------------------------------------------------------------------
        do n=1,nfields           
          if (use_aerosol_timeseries) then
            using_fixed_year_data(n) = .false.
  
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when a aerosol timeseries
!    file is to be used.
!---------------------------------------------------------------------
            if (aerosol_dataset_entry(1,n) == 1 .and. &
                aerosol_dataset_entry(2,n) == 1 .and. &
                aerosol_dataset_entry(3,n) == 1 .and. &
                aerosol_dataset_entry(4,n) == 1 .and. &
                aerosol_dataset_entry(5,n) == 1 .and. &
                aerosol_dataset_entry(6,n) == 1 ) then
              call error_mesg ('aerosol_mod', &
               'must set aerosol_dataset_entry when using  &
                                   &aerosol timeseries', FATAL)
            endif

!----------------------------------------------------------------------
!    define the offset from model base time to aerosol_dataset_entry
!    as a time_type variable.
!----------------------------------------------------------------------
            Aerosol_entry(n) = set_date (aerosol_dataset_entry(1,n), &
                                         aerosol_dataset_entry(2,n), &
                                         aerosol_dataset_entry(3,n), &
                                         aerosol_dataset_entry(4,n), &
                                         aerosol_dataset_entry(5,n), &
                                         aerosol_dataset_entry(6,n))
            call error_mesg ( 'aerosol_mod', &
             'PROCESSING AEROSOL TIMESERIES FOR ' //&
                                      & trim(aerosol_names(n)), NOTE)
            call print_date (Aerosol_entry(n) ,   &
              str= ' Data from aerosol timeseries at time: ')
            call print_date (Model_init_time , str=' This data is &
                                   &mapped to model time:')
            Aerosol_offset(n) = Aerosol_entry(n) - Model_init_time

            if (Model_init_time > Aerosol_entry(n)) then
              negative_offset(n) = .true.
            else
              negative_offset(n) = .false.
            endif

!---------------------------------------------------------------------
!    if a single entry from the timeseries is to be used, define 
!    Aerosol_entry as feb 01 of the desired year, given by the first
!    element of aerosol_dataset_entry.
!---------------------------------------------------------------------
          else
            if (aerosol_dataset_entry(1,n) == 1 .and. &
                aerosol_dataset_entry(2,n) == 1 .and. &
                aerosol_dataset_entry(3,n) == 1 .and. &
                aerosol_dataset_entry(4,n) == 1 .and. &
                aerosol_dataset_entry(5,n) == 1 .and. &
                aerosol_dataset_entry(6,n) == 1 ) then
              using_fixed_year_data(n) = .false.
              if (mpp_pe() == mpp_root_pe() ) then
               print *, 'Aerosol data for ', trim(aerosol_names(n)),   &
                  ' obtained from single year climatology file '
              endif
            else
              using_fixed_year_data(n) = .true.
              Aerosol_entry(n) = set_date (aerosol_dataset_entry(1,n), &
                                           2, 1, 0, 0, 0)
              call error_mesg ('aerosol_mod', &
                  'Aerosol data is defined from a single annual cycle &
                 &for ' // trim(aerosol_names(n)) //   &
                  &' - no interannual variation', NOTE)
              if (mpp_pe() == mpp_root_pe() ) then
                print *, 'Aerosol data for ', trim(aerosol_names(n)),  &
                  ' obtained from aerosol timeseries &
                   &for year:', aerosol_dataset_entry(1,n)
              endif
            endif
          endif
        end do

!-----------------------------------------------------------------------
!    count number of activated aerosol families. allocate a pointer 
!    array to return the names of the activated species to the calling 
!    routine.
!-----------------------------------------------------------------------
        do n=1,MAX_AEROSOL_FAMILIES
          if (family_names(n) /= ' '  ) then
            nfamilies = n
          else
            exit
          endif
        end do
        allocate (aerosol_family_names(nfamilies))
        aerosol_family_names (:) = family_names(1:nfamilies)

!-----------------------------------------------------------------------
!    initialize the interpolation module if any aerosol species have
!    been activated.
!-----------------------------------------------------------------------
        allocate (data_out_of_bounds(nfields))
        allocate (vert_interp       (nfields))
        data_out_of_bounds = CONSTANT
        vert_interp = INTERP_WEIGHTED_P
        if (trim(aerosol_data_source) == 'calculate_column') then
          do n=1,2
            if (lonb_col(n) < 0. .or. lonb_col(n) > 360.) then
              call error_mesg ('aerosol_mod', &
                   ' invalid value for lonb_col', FATAL)
            endif
            if (latb_col(n) < -90. .or. latb_col(n) > 90.) then
              call error_mesg ('aerosol_mod', &
               ' invalid value for latb_col', FATAL)
            endif
          end do
          if (time_col(1) == 0) then
            call error_mesg ('aerosol_mod', &
                   'invalid time specified for time_col', FATAL)
          endif
          call interpolator_init (Aerosol_interp, filename,  &
                                  lonb_col/RADIAN,  &
                                  latb_col/RADIAN,&
                                  data_names(:nfields),   &
                                  data_out_of_bounds=  &
                                                data_out_of_bounds, &
                                  vert_interp=vert_interp )
          Aerosol_column_time = set_date (time_col(1), time_col(2), &
                                          time_col(3), time_col(4), &
                                          time_col(5), time_col(6))
          call print_date (Aerosol_column_time,   &
             str= ' Aerosol data used is from aerosol timeseries at time: ')
          if (mpp_pe() == mpp_root_pe() ) then
            print *, 'Aerosol data is averaged over latitudes',  &
                    latb_col(1), ' to', latb_col(2), ' and longitudes',&
                    lonb_col(1), ' to', lonb_col(2)
          endif
        else
          if (nfields > 0) then
            call interpolator_init (Aerosol_interp, filename, lonb, &
                                    latb, data_names(:nfields),   &
                                    data_out_of_bounds=   &
                                                  data_out_of_bounds, &
                                    vert_interp=vert_interp )
          endif
        endif
      endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------



end subroutine aerosol_init



!######################################################################
! <SUBROUTINE NAME="aerosol_driver">
!  <OVERVIEW>
!   Interpolate aerosol verical profile based on prescribed aerosol
!   climatology input and model set up.
!  </OVERVIEW>
!  <TEMPLATE>
!   call aerosol_driver (is, js, model_time, p_half, Aerosol)
!  </TEMPLATE>
!  <INOUT NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatology input
!  </INOUT>
!  <IN NAME="model_time" TYPE="time_type">
!   The internal model simulation time, i.e. Jan. 1 1982
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   The array of model layer pressure values
!  </IN>
!  <IN NAME="is" TYPE="integer">
!   The longitude index of model physics window domain
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   The latitude index of model physics window domain
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_driver (is, js, model_time, p_half, Aerosol)

!-----------------------------------------------------------------------
!    aerosol_driver returns the names and concentrations of activated 
!    aerosol species at model grid points at the model_time to the 
!    calling routine in aerosol_type variable Aerosol. 
!-----------------------------------------------------------------------

integer,                  intent(in)     :: is,js
type(time_type),          intent(in)     :: model_time
real, dimension(:,:,:),   intent(in)     :: p_half
type(aerosol_type),       intent(inout)  :: Aerosol

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       is, js           starting subdomain i,j indices of data in 
!                        the physics_window being integrated
!       model_time       time for which aerosol data is desired
!                        [ time_type ]
!       p_half           model pressure at interface levels 
!                        [ Pa ]
!      
!   intent(inout) variables:
!
!       Aerosol    aerosol_type variable. the following components will
!                  be returned from this routine:
!                   aerosol      concentration of each active aerosol 
!                                species at each model grid point
!                                [ kg / m**2 ]
!                   aerosol_names 
!                                names assigned to each active species
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension(1,1, size(p_half,3)-1, nfields) :: aerosol_data
      real, dimension(1,1, size(p_half,3))            :: p_half_col
      type(time_type) :: Aerosol_time          ! time for which data is
                                               ! obtained from aerosol
                                               ! timeseries
      logical         :: make_separate_calls   ! aerosol interpolation
                                               ! to be done one at a 
                                               ! time
      integer         :: n, k, j, i, na            ! do-loop index

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

      if (do_column_aerosol) then
 
!---------------------------------------------------------------------
!    allocate the components of the aerosol_type variable. here all
!    aerosol is consolidated into a single variable.
!----------------------------------------------------------------------
        allocate (Aerosol%aerosol_names (nfields))
        Aerosol%aerosol_names(1) = 'total_aerosol'
        allocate (Aerosol%family_members(nfields, nfamilies))
        allocate (Aerosol%aerosol(size(p_half,1), size(p_half,2), &
                                  size(p_half,3) - 1, nfields))
        do k=1, size(Aerosol%aerosol,3)
          Aerosol%aerosol(:,:,k,1) = specified_aerosol(k)
        end do
      else
      
!--------------------------------------------------------------------
!    allocate and define an array to hold the activated aerosol names.
!---------------------------------------------------------------------
      allocate (Aerosol%aerosol_names (nfields))
      Aerosol%aerosol_names = data_names(:nfields) 

!--------------------------------------------------------------------
!    allocate and define an array which defines the members of the
!    requested aerosol families. the additional field is the volcanic 
!    aerosol.
!---------------------------------------------------------------------
      allocate (Aerosol%family_members(nfields+1, nfamilies))
      if (nfamilies > 0) then
        do n=1,nfamilies
          do na = 1, nfields
          select case(n)
            case (1)
             Aerosol%family_members(na,1) = in_family1(na)
            case (2)
             Aerosol%family_members(na,2) = in_family2(na)
            case (3)
             Aerosol%family_members(na,3) = in_family3(na)
            case (4)
             Aerosol%family_members(na,4) = in_family4(na)
            case (5)
             Aerosol%family_members(na,5) = in_family5(na)
            case (6)
             Aerosol%family_members(na,6) = in_family6(na)
            case (7)
             Aerosol%family_members(na,7) = in_family7(na)
            case (8)
             Aerosol%family_members(na,8) = in_family8(na)
            case (9)
             Aerosol%family_members(na,9) = in_family9(na)
            case DEFAULT
          end select
        end do
        if (volc_in_fam_col_opt_depth(n)) then
          Aerosol%family_members(nfields+1,n) = .true.
        else
          Aerosol%family_members(nfields+1,n) = .false.
        endif
        end do
      endif
      

!---------------------------------------------------------------------
!    allocate an array to hold the aerosol amounts for each species at
!    each grid point. 
!----------------------------------------------------------------------
      allocate (Aerosol%aerosol(size(p_half,1), size(p_half,2), &
                                size(p_half,3) - 1, nfields))
      
!----------------------------------------------------------------------
!    determine if separate calls to interpolator  must be made for 
!    each aerosol species, or if all variables in the file may be
!    interpolated together.  reasons for separate calls include differ-
!    ent data times desired for different aerosols, different vertical
!    interpolation procedures and different treatment of undefined
!    data.
!----------------------------------------------------------------------
      make_separate_calls = .false.
      do n=2,nfields
        if (using_fixed_year_data(n) .and.   &
            (.not. using_fixed_year_data(n-1) ) ) then
          make_separate_calls = .true.
          exit
        endif
        if (Aerosol_entry(n) /= Aerosol_entry(n-1)) then
          make_separate_calls = .true.
          exit
        endif
        if (data_out_of_bounds(n) /= data_out_of_bounds(n-1)) then
          make_separate_calls = .true.
          exit
        endif
        if (vert_interp       (n) /= vert_interp       (n-1)) then
          make_separate_calls = .true.
          exit
        endif
      end do

!--------------------------------------------------------------------
!    if separate calls are required for each aerosol species, loop over
!    the individual species.
!--------------------------------------------------------------------
      if (make_separate_calls) then
        do n=1,nfields

!--------------------------------------------------------------------
!    if the data timeseries is to be used for species n, define the
!    time for which data is desired, and then call interpolator to
!    obtain the data.
!--------------------------------------------------------------------
          if (use_aerosol_timeseries) then

!----------------------------------------------------------------------
!    if separate calls are required, define the Aerosol_time for 
!    aerosol n and call interpolator. store the aerosol amount in 
!    Aerosol%aerosol.
!----------------------------------------------------------------------
            if (negative_offset(n)) then
              Aerosol_time = model_time - Aerosol_offset(n)
            else
              Aerosol_time = model_time + Aerosol_offset(n)
            endif
            call interpolator (Aerosol_interp, Aerosol_time, p_half, &
                                 Aerosol%aerosol(:,:,:,n),    &
                                 Aerosol%aerosol_names(n), is, js)

!--------------------------------------------------------------------
!    if the data timeseries is not to be used for species n, define the
!    time for which data is desired, and then call interpolator to
!    obtain the data.
!--------------------------------------------------------------------
          else  ! (use_aerosol_timeseries)

!---------------------------------------------------------------------
!    if a fixed year has not been specified, obtain data relevant for
!    the current model year.
!---------------------------------------------------------------------
            if ( .not. using_fixed_year_data(n)) then
              call interpolator (Aerosol_interp, model_time, p_half, &
                                 Aerosol%aerosol(:,:,:,n),    &
                                 Aerosol%aerosol_names(n), is, js)

!----------------------------------------------------------------------
!    if a fixed year has been specified, call set_aerosol_time to define
!    the Aerosol_time to be used for aerosol n. call interpolator to 
!    obtain the aerosol values and store the aerosol amounts in 
!    Aerosol%aerosol.
!----------------------------------------------------------------------
            else 
              call set_aerosol_time (model_time, Aerosol_entry(n), &
                                     Aerosol_time)
              call interpolator (Aerosol_interp, Aerosol_time, p_half, &
                                 Aerosol%aerosol(:,:,:,n),    &
                                 Aerosol%aerosol_names(n), is, js)
            endif  
          endif ! (use_aerosol_timeseries)
        end do  !(nfields)

!----------------------------------------------------------------------
!    if separate calls are not required, use the first aerosol char-
!    acteristics to define Aerosol_time and make a single call to 
!    interpolator. store the aerosol amounts in Aerosol%aerosol.
!----------------------------------------------------------------------
      else ! (make_separate_calls)

!--------------------------------------------------------------------
!    if the data timeseries is to be used for species n, define the
!    time for which data is desired, and then call interpolator to
!    obtain the data.
!--------------------------------------------------------------------
        if (use_aerosol_timeseries) then
          if (negative_offset(1)) then
              Aerosol_time = model_time - Aerosol_offset(1)
            else
              Aerosol_time = model_time + Aerosol_offset(1)
            endif

!--------------------------------------------------------------------
!    if 'calculate_column' is being used, obtain the aerosol values for
!    each column, one at a time, using the pressure profile for that
!    column. this allows each column to see the same aerosol fields,
!    but distributed appropriately for its pressure structure.
!--------------------------------------------------------------------
            if (trim(aerosol_data_source) == 'calculate_column') then
              do j=1, size(p_half,2)
                do i=1, size(p_half,1)
                  p_half_col(1,1,:) = p_half(i,j,:)
                  call interpolator (Aerosol_interp,&
                                     Aerosol_column_time,  &
                                     p_half_col, aerosol_data, &
!                                    Aerosol%aerosol_names(1), is, js)
                                     Aerosol%aerosol_names(1), 1, 1  )
                  Aerosol%aerosol(i,j,:,:) = aerosol_data(1,1,:,:)
                end do
              end do
            else

              call interpolator (Aerosol_interp, Aerosol_time, p_half, &
                                 Aerosol%aerosol,    &
                                 Aerosol%aerosol_names(1), is, js)
            endif 

!--------------------------------------------------------------------
!    if the data timeseries is not to be used for species n, define the
!    time for which data is desired, and then call interpolator to
!    obtain the data.
!--------------------------------------------------------------------
        else 

!---------------------------------------------------------------------
!    if a fixed year has not been specified, obtain data relevant for
!    the current model year.
!---------------------------------------------------------------------
          if (.not. using_fixed_year_data(1)) then
            call interpolator (Aerosol_interp, model_time, p_half, &
                               Aerosol%aerosol,    &
                               Aerosol%aerosol_names(1), is, js)

!----------------------------------------------------------------------
!    if a fixed year has been specified, call set_aerosol_time to define
!    the Aerosol_time. call interpolator to obtain the aerosol values 
!    and store the aerosol amounts in Aerosol%aerosol.
!----------------------------------------------------------------------
          else
            call set_aerosol_time (model_time, Aerosol_entry(1), &
                                   Aerosol_time)
            call interpolator (Aerosol_interp, Aerosol_time, p_half, &
                               Aerosol%aerosol,    &
                               Aerosol%aerosol_names(1), is, js)
          endif ! (using_fixed_year)
        endif  ! (use_aerosol_timeseries)
      endif ! (make_separate_calls   )
    endif  ! (do_column_aerosol)

!-------------------------------------------------------------------- 

end subroutine aerosol_driver



!#####################################################################
! <SUBROUTINE NAME="aerosol_end">
!  <OVERVIEW>
!   aerosol_end is the destructor for aerosol_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   aerosol_end is the destructor for aerosol_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosol_end

!----------------------------------------------------------------------
!    aerosol_end is the destructor for aerosol_mod.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!---------------------------------------------------------------------
!    call interpolator_end to release the interpolate_type variable 
!    used in this module.
!---------------------------------------------------------------------
      if (.not. do_column_aerosol) then
        if (nfields > 0) then
          call interpolator_end (Aerosol_interp)
        endif
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------



end subroutine aerosol_end


!#####################################################################

subroutine set_aerosol_time (Model_time, Entry, Aerosol_time)

type(time_type), intent(in)   :: Model_time, Entry
type(time_type), intent(out)  :: Aerosol_time

      integer :: mo_yr, yr, mo, dy, hr, mn, sc, dum, dayspmn

      call get_date (Model_time, mo_yr, mo, dy, hr, mn, sc)
      call get_date (Entry, yr, dum,dum,dum,dum,dum)
      if (mo ==2 .and. dy == 29) then
        dayspmn = days_in_month(Entry)
        if (dayspmn /= 29) then
          Aerosol_time = set_date (yr, mo, dy-1, hr, mn, sc)
        else
          Aerosol_time = set_date (yr, mo, dy, hr, mn, sc)
        endif
      else
        Aerosol_time = set_date (yr, mo, dy, hr, mn, sc)
      endif

!--------------------------------------------------------------------


end subroutine set_aerosol_time


!#####################################################################
! <SUBROUTINE NAME="aerosol_dealloc">
!  <OVERVIEW>
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_dealloc
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosol_dealloc (Aerosol)

!---------------------------------------------------------------------
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!---------------------------------------------------------------------

type(aerosol_type), intent(inout) :: Aerosol

!---------------------------------------------------------------------
!  intent(inout) variables:
! 
!      Aerosol       aerosol_type variable containing information on
!                    the activated aerosol species
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!---------------------------------------------------------------------
!    deallocate the components of the aerosol_type variable.
!---------------------------------------------------------------------
      deallocate (Aerosol%aerosol)
      deallocate (Aerosol%aerosol_names)
      deallocate (Aerosol%family_members)
 
!----------------------------------------------------------------------


end subroutine aerosol_dealloc


!#####################################################################

! <SUBROUTINE NAME="obtain_input_file_data">
!  <OVERVIEW>
!   obtain_input_file_data reads an input file containing a single
!    column aerosol profile. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   obtain_input_file_data reads an input file containing a single
!    column aerosol profile. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  obtain_input_file_data
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine obtain_input_file_data 

!---------------------------------------------------------------------
!    obtain_input_file_data reads an input file containing a single
!    column aerosol profile.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer   :: iounit    ! unit to read file on
      integer   :: kmax_file ! number of levels of data in file
      integer   :: k         ! do-loop index
      character*31, dimension(200) :: dimnam
      integer(kind=4), dimension(200) :: dimsiz
      integer(kind=4)                 :: ncid, rcode, nvars, ndims, &
                                         ngatts, recdim
      integer    :: i, j
      integer, PARAMETER :: MAXDIMS = 10
      integer(kind=4), dimension(MAXDIMS) :: start, count, vdims
      integer(kind=4)                     :: ivarid, ntp, nvdim, nvs, &
                                             ndsize
      character*31   dummy
      


!-------------------------------------------------------------------
!    determine if the input data input file exists in ascii format. if 
!    so, read the number of data records in the file.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/id1aero') ) then
        iounit = open_namelist_file ('INPUT/id1aero')
        read (iounit,FMT = '(i4)') kmax_file

!-------------------------------------------------------------------
!    allocate space for the input data. read the data set. close the 
!    file upon completion.
!---------------------------------------------------------------------
         allocate (specified_aerosol(kmax_file) )
         read (iounit,FMT = '(5e18.10)')   &
                          (specified_aerosol(k),k=1,kmax_file)
         call close_file (iounit)

!-------------------------------------------------------------------
!    determine if a netcdf input data file exists. if so, read the 
!    number of data records in the file.
!---------------------------------------------------------------------
      else if (file_exist ( 'INPUT/id1aero.nc') ) then
        ncid = ncopn ('INPUT/id1aero.nc', 0, rcode)
        call ncinq (ncid, ndims, nvars, ngatts, recdim, rcode)
        do i=1,ndims
          call ncdinq (ncid, i, dimnam(i), dimsiz(i), rcode)
          if (dimnam(i) == 'lev') then
            kmax_file = dimsiz(i)
          endif
        end do
             
!-------------------------------------------------------------------
!    allocate space for the input data. read the data set. close the 
!    file upon completion.
!---------------------------------------------------------------------
        allocate (specified_aerosol(kmax_file) )
        ivarid = ncvid(ncid, 'aerosol', rcode)
        call ncvinq (ncid, ivarid, dummy, ntp, nvdim, vdims, nvs, rcode)
        do j=1,nvdim
          call ncdinq (ncid, vdims(j), dummy, ndsize, rcode)
          start(j) = 1
          count(j) = ndsize
        end do
       call ncvgt (ncid, ivarid, start, count, specified_aerosol, rcode)

         call ncclos (ncid, rcode)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
       else
         call error_mesg ( 'aerosol_mod', &
              'desired aerosol input file is not present',FATAL)
       endif

!----------------------------------------------------------------------


end subroutine obtain_input_file_data 


!###################################################################### 



                  end module aerosol_mod 



!=======================================================================



#ifdef test_aerosol

program main

use aerosol_mod
use mpp_mod
use mpp_io_mod
use mpp_domains_mod
use time_manager_mod
use diag_manager_mod
use rad_utilities_mod




implicit none

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer, parameter :: NLON=20, NLAT=10,NLEV=8
integer, parameter :: MAX_AERSOL_NAMES = 100
real :: latb(NLAT+1),lonb(NLON+1),pi,phalf(NLON,NLAT,NLEV+1)
integer :: i,nspecies
type(time_type) :: model_time
character(len=64), dimension(MAX_AEROSOL_NAMES) :: names
type(aerosol_type)  :: Aerosol

pi = 4.*atan(1.)

call mpp_init
call mpp_io_init
call mpp_domains_init
call diag_manager_init
call set_calendar_type(JULIAN)

do i = 1,NLAT+1
   latb(i) = -90. + 180.*REAL(i-1)/REAL(NLAT)
end do
do i = 1,NLON+1
   lonb(i) = -180. + 360.*REAL(i-1)/REAL(NLON)
end do

latb(:) = latb(:) * pi/180.
lonb(:) = lonb(:) * pi/180.

do i = 1,NLEV+1
   phalf(:,:,i) = 101325. * REAL(i-1) / REAL(NLEV)
end do

model_time = set_date(1980,1,1,0,0,0)

call aerosol_init (lonb, latb, names)

call aerosol_driver (1,1,model_time, phalf, Aerosol)

call aerosol_dealloc (Aerosol)

call aerosol_end

call mpp_exit

end program main

#endif
