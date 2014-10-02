                 module cloud_spec_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
! </OVERVIEW>
! <DESCRIPTION>
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present, or from a prescribed formula based on
!    prescribed water paths for high, middle and low clouds.
! </DESCRIPTION>

!   shared modules:

use time_manager_mod,         only: time_type, time_manager_init, &
                                    set_time, operator (+)
use fms_mod,                  only: open_namelist_file, mpp_pe, &
                                    mpp_root_pe, stdlog,  fms_init, &
                                    mpp_clock_id, mpp_clock_begin, &
                                   mpp_clock_end, CLOCK_ROUTINE, &
                                   MPP_CLOCK_SYNC, &
                                    write_version_number, file_exist, & 
                                    check_nml_error, error_mesg,   &
                                    FATAL, NOTE, WARNING, close_file
use tracer_manager_mod,       only:         &
!                                   tracer_manager_init,  &
                                    get_tracer_index
use field_manager_mod,        only:       &
                                    field_manager_init, &
                                    MODEL_ATMOS
use data_override_mod,        only: data_override

! shared radiation package modules:

use rad_utilities_mod,        only: rad_utilities_init, &
                                    cld_specification_type, &
                                    atmos_input_type, &
                                    surface_type, &
                                    Rad_control, &
                                    microphysics_type,  &         
                                    Cldrad_control
use esfsw_parameters_mod,     only: esfsw_parameters_init, Solar_spect

! interface modules to various cloud parameterizations:

use strat_clouds_W_mod,       only: strat_clouds_W_init,   &
                                    strat_clouds_amt, strat_clouds_W_end
use diag_clouds_W_mod,        only: diag_clouds_W_init,   &
                                    diag_clouds_amt, &
                                    diag_clouds_W_end
use zetac_clouds_W_mod,       only: zetac_clouds_W_init,   &
                                    zetac_clouds_amt, &
                                    zetac_clouds_W_end
use specified_clouds_W_mod,   only: specified_clouds_W_init, &
                                    specified_clouds_amt, &
                                    specified_clouds_W_end
use rh_based_clouds_mod,      only: rh_based_clouds_init,  &
                                    rh_clouds_amt, &
                                    rh_based_clouds_end
use donner_deep_clouds_W_mod, only: donner_deep_clouds_W_init, &
                                    donner_deep_clouds_amt, &
                                    donner_deep_clouds_W_end
use mgrp_prscr_clds_mod,      only: mgrp_prscr_clds_init, &
                                    prscr_clds_amt,  & 
                                    mgrp_prscr_clds_end 
use standalone_clouds_mod,    only: standalone_clouds_init, &
                                    standalone_clouds_amt, &
                                    standalone_clouds_end
                                 
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present, or from a prescribed formula based on
!    prescribed water paths for high, middle and low clouds.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloud_spec.f90,v 12.0 2005/04/14 15:44:28 fms Exp $'
character(len=128)  :: tagname =  '$Name: lima $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloud_spec_init, cloud_spec,    &
         cloud_spec_dealloc, cloud_spec_end

private    &

!  called from cloud_spec:
         initialize_cldamts, microphys_presc_conc,  &
         combine_cloud_properties


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)  ::      &
              cloud_type_form = '     ' ! cloud parameterization being 
                                        ! used; either 'strat', 'rh', 
                                        ! 'deep',  'stratdeep', 'zonal',
                                        ! 'obs', 'prescribed', 'diag', 
                                        ! 'none', 'specified', 'zetac'
                                        ! 'specified_strat'
                                        ! or 'not_sea_esf'       
real :: wtr_cld_reff=10.                ! assumed cloud drop efective
                                        ! radius [ microns ]  
real :: ice_cld_reff=50.                ! assumed ice cloud effective
                                        ! size [ microns ]
real :: rain_reff=250.                  ! assumed rain drop effective
                                        ! radius [ microns ]
character(len=16) :: overlap_type = 'random'    
                                        ! cloud overlap assumption; 
                                        ! allowable values are 'random'
                                        ! or 'max-random'  
logical :: doing_data_override=.false.

namelist /cloud_spec_nml / cloud_type_form, wtr_cld_reff,   &
                           ice_cld_reff, rain_reff, overlap_type, &
                           doing_data_override

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!    assumed water paths.
!--------------------------------------------------------------------
real   ::  lwpath_hi  = 6.313929   ! assumed water path for high clouds
                                   ! [ grams / m**2 ]
real   ::  lwpath_mid = 18.94179   ! assumed water path for middle 
                                   ! clouds [ grams / m**2 ]
real   ::  lwpath_low = 75.76714   ! assumed water path for low clouds
                                   ! [ grams / m**2 ]

!---------------------------------------------------------------------
!    logical  flags.

logical :: module_is_initialized = .false.   ! module initialized ?

!---------------------------------------------------------------------
!    time-step related constants.

integer :: num_pts       !  number of grid columns processed so far that
                         !  have cloud data present (used to identify
                         !  module coldstart condition)
integer :: tot_pts       !  total number of grid columns in the 
                         !  processor's domain

!---------------------------------------------------------------------
!     indices for cloud tracers

integer :: nql           ! tracer index for liquid water
integer :: nqi           ! tracer index for ice water
integer :: nqa           ! tracer index for cloud area

!---------------------------------------------------------------------
!     miscellaneous variables:

integer :: num_slingo_bands  ! number of radiative bands over which 
                             ! cloud optical depth is calculated in the
                             ! gordon diag_cloud parameterization
integer :: id, jd, kmax

!----------------------------------------------------------------------
!----------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="cloud_spec_init">
!  <OVERVIEW>
!   Contructor of cloud_spec_package module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Contructor of cloud_spec_package module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec_init ( pref, lonb, latb, axes, Time)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   reference pressure levels containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   the longitude array of the model grid point
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   the latitude array of the model grid point
!  </IN>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
! 
subroutine cloud_spec_init (pref, lonb, latb, axes, Time)

!---------------------------------------------------------------------
!    cloud_spec_init is the constructor for cloud_spec_mod.
!---------------------------------------------------------------------

real, dimension(:,:),     intent(in)   ::  pref        
real, dimension(:),       intent(in)   ::  lonb, latb
integer, dimension(4),    intent(in)   ::  axes
type(time_type),          intent(in)   ::  Time

!-------------------------------------------------------------------
!    intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      array of model longitudes on cell boundaries 
!                 [ radians ]
!       latb      array of model latitudes at cell boundaries [radians]
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
 
      integer   ::   unit, ierr, io
      integer   ::   ndum

!--------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!      ndum     dummy argument needed for call to field_manager_init
!
!--------------------------------------------------------------------
 
!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call rad_utilities_init
      call field_manager_init (ndum)
      call esfsw_parameters_init
!  not yet compliant:
!     call tracer_manager_init  ! not public
 
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=cloud_spec_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'cloud_spec_nml')
        enddo
10      call close_file (unit)
      endif

!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
           write (stdlog(), nml=cloud_spec_nml)

      id = size(lonb(:)) - 1
      jd = size(latb(:)) - 1
      kmax = size(pref,1) - 1

!--------------------------------------------------------------------
!    verify a valid type of cloud overlap. set logical variables
!    based on the namelist value.
!--------------------------------------------------------------------
      if (trim(overlap_type) == 'random') then
        Cldrad_control%do_random_overlap = .true.
      else if (trim(overlap_type) == 'max-random') then
        Cldrad_control%do_max_random_overlap = .true.
      else
        call error_mesg ('cloud_spec_mod',  &
         ' invalid specification of overlap_type', FATAL)
      endif

!-------------------------------------------------------------------
!    set the variables indicating that the above control variables have
!    been set.
!--------------------------------------------------------------------
      Cldrad_control%do_random_overlap_iz = .true.
      Cldrad_control%do_max_random_overlap_iz = .true.

!--------------------------------------------------------------------
!    if the sea-esf radiation package is not being used, then 
!    cloud_type_form will have been set to 'not_sea_esf'. in such a
!    case, the clouds will be specified internally within the 
!    radiation_driver_mod, so simply return.
!--------------------------------------------------------------------
      if (trim(cloud_type_form) == 'not_sea_esf')  return
        
!-------------------------------------------------------------------
!    verify that the nml variable cloud_type_form specifies a valid
!    cloud parameterization. set the appropriate logical control
!    variable(s) to .true.. call the constructor modules for the
!    specific cloud scheme(s) requested.
!-------------------------------------------------------------------
      if (trim(cloud_type_form) == 'strat')  then

!-------------------------------------------------------------------
!    cloud fractions, heights are predicted by the model based on klein 
!    parameterization. strat is an acceptable option both for standalone
!    and gcm applications.
!-------------------------------------------------------------------
        Cldrad_control%do_strat_clouds = .true.
        call strat_clouds_W_init(latb, lonb)

      else if (trim(cloud_type_form) == 'specified_strat')  then
        Cldrad_control%do_specified_strat_clouds = .true.
        Cldrad_control%do_strat_clouds = .true.
        call strat_clouds_W_init(latb, lonb)
        call standalone_clouds_init (pref, lonb, latb)

!-------------------------------------------------------------------
!    cloud fractions, heights are diagnosed based on model relative 
!    humidity.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'rh')   then
            Cldrad_control%do_rh_clouds = .true.
            call rh_based_clouds_init 

!-------------------------------------------------------------------
!    cloud fractions, heights are predicted by the donner deep cloud 
!    (cell cloud, anvil cloud) scheme.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'deep')  then
            Cldrad_control%do_donner_deep_clouds = .true.
            call donner_deep_clouds_W_init (pref, lonb, latb,   &
                                            axes, Time)

!-------------------------------------------------------------------
!    cloud fractions, heights are a combination of the donner
!    deep cloud (cell cloud, anvil cloud) and klein large-scale cloud
!    parameterizations.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'stratdeep')  then
            Cldrad_control%do_strat_clouds = .true.
            Cldrad_control%do_donner_deep_clouds = .true.
            call strat_clouds_W_init(latb, lonb)
            call donner_deep_clouds_W_init (pref, lonb, latb,   &
                                            axes, Time)

!-------------------------------------------------------------------
!    cloud fractions, heights are prescribed as zonally uniform using
!    the original fms specification.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'zonal')  then
            Cldrad_control%do_zonal_clouds = .true.
            call specified_clouds_W_init (lonb, latb)

!-------------------------------------------------------------------
!    cloud fractions, heights are based on observed data set.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'obs')  then
            Cldrad_control%do_obs_clouds = .true.
            call specified_clouds_W_init (lonb, latb)

!-------------------------------------------------------------------
!    cloud fractions, heights are prescribed as zonally invariant, using
!    the formulation from skyhi, with the ability to use a prescribed
!    microphysics.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'prescribed')  then
            Cldrad_control%do_mgroup_prescribed = .true.
            call mgrp_prscr_clds_init (pref, latb)

!-------------------------------------------------------------------
!    model is run with gordon diagnostic clouds.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'diag')  then
            Cldrad_control%do_diag_clouds = .true.
            call diag_clouds_W_init (num_slingo_bands)

!-------------------------------------------------------------------
!    model is run with zetac clouds.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'zetac')  then
            Cldrad_control%do_zetac_clouds = .true.
            call zetac_clouds_W_init 
 
!-------------------------------------------------------------------
!    model is run without clouds.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'none')  then
            Cldrad_control%do_no_clouds = .true.

!-------------------------------------------------------------------
!    model is run with specified clouds and cloud properties.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'specified')  then
            Cldrad_control%do_specified_clouds = .true.
            call standalone_clouds_init (pref, lonb, latb)

!-------------------------------------------------------------------
!    failure message if none of the above options was chosen.
!-------------------------------------------------------------------
          else
            call error_mesg ('cloud_spec_mod',  &
              'invalid cloud_type_form specified', FATAL)
      endif  ! (strat)

!--------------------------------------------------------------------
!    define the dimensions of the model subdomain assigned to the 
!    processor.
!--------------------------------------------------------------------
      tot_pts = (size(latb(:))-1)*(size(lonb(:))-1)

!--------------------------------------------------------------------
!    determine if the current run is cold-starting this module. if a 
!    restart file is present, then this is not a coldstart. in that case
!    set num_pts to tot_pts so that if cloud data is not available an 
!    error message can be generated. if this is a coldstart, cloud data
!    will not be available until num_pts equals or exceeds tot_pts, so
!    continue processing without issuing an error message. 
!--------------------------------------------------------------------
      if (file_exist ('INPUT/tracer_cld_amt.res') .or.  &
          file_exist ('INPUT/strat_cloud.res') ) then
        num_pts = tot_pts
      else
        num_pts = 0
      endif

!---------------------------------------------------------------------
!    obtain the tracer indices for the strat_cloud variables when
!    running gcm.
!---------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds .or.  &
            Cldrad_control%do_zetac_clouds) then
          nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
          nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
          nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
          if (mpp_pe() == mpp_root_pe()) &
            write (stdlog(),'(a,3i4)') 'Stratiform cloud tracer ind&
                &ices: nql,nqi,nqa =',nql,nqi,nqa
          if (min(nql,nqi,nqa) <= 0)   &
             call error_mesg ('cloud_spec_mod', &
             'stratiform cloud tracer(s) not found', FATAL)
          if (nql == nqi .or. nqa == nqi .or. nql == nqa)   &
              call error_mesg ('cloud_spec_mod',  &
            'tracers indices cannot be the same (i.e., nql=nqi=nqa).', &
                                                              FATAL)
        endif

!---------------------------------------------------------------------
!    define the variables indicating that the cloud parameterization
!    control variables have been defined.
!---------------------------------------------------------------------
      Cldrad_control%do_rh_clouds_iz = .true.
      Cldrad_control%do_strat_clouds_iz = .true.
      Cldrad_control%do_zonal_clouds_iz = .true.
      Cldrad_control%do_mgroup_prescribed_iz = .true.
      Cldrad_control%do_obs_clouds_iz = .true.
      Cldrad_control%do_no_clouds_iz = .true.
      Cldrad_control%do_diag_clouds_iz = .true.
      Cldrad_control%do_specified_clouds_iz = .true.
      Cldrad_control%do_specified_strat_clouds_iz = .true.
      Cldrad_control%do_donner_deep_clouds_iz = .true.
      Cldrad_control%do_zetac_clouds_iz = .true.
      Cldrad_control%do_stochastic_clouds_iz = .true.
 
!---------------------------------------------------------------------
!    mark the module initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------



end subroutine cloud_spec_init



!######################################################################
! <SUBROUTINE NAME="cloud_spec">
!  <OVERVIEW>
!    cloud_radiative_properties defines the cloud radiative properties 
!    appropriate for the radiation options that are active.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_radiative_properties defines the cloud radiative properties 
!    appropriate for the radiation options that are active.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec (is, ie, js, je, lat, z_half, z_full, Rad_time,
!                       Atmos_input, &
!                       Surface, Cld_spec, Lsc_microphys,  &
!                       Meso_microphys, Cell_microphys, cloud_water_in, &
!                       cloud_ice_in, cloud_area_in, r, kbot, mask)
!  </TEMPLATE>
!  <IN NAME="is,ie,js,je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!   time at which radiation calculation is to apply
!  </IN>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </INOUT>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </INOUT>
!  <INOUT NAME="Surface" TYPE="Surface">
!   Surface boundary condition to radiation package
!  </INOUT>
!  <IN NAME="cloud_water_in" TYPE="real">
!   OPTIONAL: cloud water mixing ratio  present when running 
!    standalone columns or sa_gcm
!  </IN>
!  <IN NAME="cloud_ice_in" TYPE="real">
!   OPTIONAL: cloud ice mixing ratio  present when running 
!    standalone columns or sa_gcm
!  </IN>
!  <IN NAME="cloud_area_in" TYPE="real">
!   OPTIONAL: fractional cloud area, present when running 
!                        standalone columns or sa_gcm
!  </IN>
!  <IN NAME="r" TYPE="real">
!   OPTIONAL: model tracer fields on the current time step
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine cloud_spec (is, ie, js, je, lat, z_half, z_full, Rad_time, &
                       Atmos_input, Surface, Cld_spec, Lsc_microphys, &
                       Meso_microphys, Cell_microphys, cloud_water_in, &
                       cloud_ice_in, cloud_area_in, r, kbot, mask)

!----------------------------------------------------------------------
!    cloud_spec specifies the cloud field seen by the radiation package.
!----------------------------------------------------------------------


!----------------------------------------------------------------------
integer,                      intent(in)             :: is, ie, js, je
real, dimension(:,:),         intent(in)             :: lat
real, dimension(:,:,:),       intent(in)             :: z_half, z_full
type(time_type),              intent(in)             :: Rad_time
type(atmos_input_type),       intent(inout)          :: Atmos_input
type(surface_type),           intent(inout)          :: Surface       
type(cld_specification_type), intent(inout)          :: Cld_spec    
type(microphysics_type),      intent(inout)          :: Lsc_microphys, &
                                                        Meso_microphys,&
                                                        Cell_microphys
real, dimension(:,:,:),       intent(in),   optional :: cloud_water_in,&
                                                        cloud_ice_in, &
                                                        cloud_area_in
real, dimension(:,:,:,:),     intent(in),   optional :: r
integer, dimension(:,:),      intent(in),   optional :: kbot
real, dimension(:,:,:),       intent(in),   optional :: mask
!-------------------------------------------------------------------
 
!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      lat               latitude of model points  [ radians ]
!      z_half            height asl at half levels [ m ]
!      z_full            height asl at full levels [ m ]
!      Rad_time          time at which radiation calculation is to apply
!                        [ time_type (days, seconds) ] 
!
!   intent(inout) variables:
!
!      Atmos_input       atmospheric input fields on model grid,
!                        [ atmos_input_type ] 
!      Surface           variables defining the surface albedo and land
!                        fraction
!                        [ surface_type ]
!      Cld_spec          variables on the model grid which define all or
!                        some of the following, dependent on the 
!                        specific cloud parameterization: cloud optical 
!                        paths, particle sizes, cloud fractions, cloud 
!                        thickness, number of clouds in a column, 
!                        and /or cloud type (high/mid/low, ice/liq or 
!                        random/max overlap)
!                        [ cld_specification_type ]
!      Lsc_microphys     variables describing the microphysical proper-
!                        ties of the large-scale clouds
!                        [ microphysics_type ]
!      Meso_microphys    variables describing the microphysical proper-
!                        ties of the meso-scale clouds
!                        [ microphysics_type ]
!      Cell_microphys    variables describing the microphysical proper-
!                        ties of the convective cell-scale clouds
!                        [ microphysics_type ]
!
!   intent(in), optional variables:
!
!      cloud_water_in    cloud water mixing ratio (or specific humidity 
!                        ????), present when running standalone columns
!                        or sa_gcm
!                        [ non-dimensional ]
!      cloud_ice_in      cloud ice mixing ratio (or specific humidity 
!                         ????), present when running standalone columns
!                        or sa_gcm
!                        [ non-dimensional ]
!      cloud_area_in     fractional cloud area, present when running 
!                        standalone columns or sa_gcm
!                        [ non-dimensional ]
!      r                 model tracer fields on the current time step
!      kbot              present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!      mask              present when running eta vertical coordinate,
!                        mask to remove points below ground
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: ix, jx, kx
      integer   :: ierr
      logical   :: override
      type(time_type) :: Data_time
      real, dimension (id, jd, kmax) :: ql_proc, qi_proc, qa_proc

!---------------------------------------------------------------------
!   local variables:
!
!        ix      number of grid points in x direction (on processor)
!        jx      number of grid points in y direction (on processor)
!        kx      number of model layers
!        ierr
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    check for the presence of optional input arguments.
!---------------------------------------------------------------------
      ierr = 1
      if (present(cloud_water_in) .and.   &
          present(cloud_ice_in) .and. &
          present(cloud_area_in)) then
        ierr = 0
      else
        if (present (r)) then
          ierr = 0
        else
          call error_mesg ('cloud_spec_mod', &
              'must input either r or cloud_water_in and  &
              &cloud_ice_in when using predicted cloud microphysics', &
                                                               FATAL)
        endif
      endif

!----------------------------------------------------------------------
!    define model dimensions.
!----------------------------------------------------------------------
      ix = size(Atmos_input%deltaz,1)
      jx = size(Atmos_input%deltaz,2)
      kx = size(Atmos_input%deltaz,3)

!----------------------------------------------------------------------
!    call initialize_cldamts to allocate and initialize the arrays
!    contained in the structures used to specify the cloud amounts, 
!    types and locations and the microphysical parameters.
!----------------------------------------------------------------------
      call initialize_cldamts (ix, jx, kx, Lsc_microphys,   &
                               Meso_microphys, Cell_microphys, Cld_spec)

!---------------------------------------------------------------------
!    define the cloud_water, cloud_ice and cloud_area components of 
!    Cld_spec.
!---------------------------------------------------------------------
      if (present (cloud_ice_in) .and. &
          present (cloud_water_in) ) then
        Cld_spec%cloud_ice   = cloud_ice_in
        Cld_spec%cloud_water = cloud_water_in
      endif
      if (present (cloud_area_in)  )then
        Cld_spec%cloud_area = cloud_area_in
      endif

!----------------------------------------------------------------------
!    if a cloud scheme is activated (in contrast to running without any
!    clouds), call the appropriate subroutine to define the cloud
!    location, type, amount or whatever other arrays the particular 
!    parameterization uses to specify its clouds. if the model is being
!    run with do_no_clouds = .true., exit from this routine, leaving
!    the cloud specification variables as they were initialized (to a
!    condition of no clouds).
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!---------------------------------------------------------------------
!    when running in standalone columns mode, call standalone_clouds_amt
!    to obtain the cloud specification variables. 
!---------------------------------------------------------------------
        if (Cldrad_control%do_specified_clouds .or. &
            Cldrad_control%do_specified_strat_clouds )   then

          call standalone_clouds_amt (is, ie, js, je, lat,     &
                                      Atmos_input%press, Cld_spec)

!---------------------------------------------------------------------
!    if the rh diagnostic cloud scheme is active, call rh_clouds_amt
!    to define the needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_rh_clouds) then
          call rh_clouds_amt (is, ie, js, je, Atmos_input%press, lat,  &
                              Cld_spec)

!---------------------------------------------------------------------
!    if either zonal clouds or obs clouds is active, call 
!    specified_clouds_amt to obtain the needed cloud specification
!    variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_zonal_clouds .or. &
                 Cldrad_control%do_obs_clouds) then
          call specified_clouds_amt (is, ie, js, je, Rad_time, lat,    &
                                     Atmos_input%pflux, Cld_spec)

!---------------------------------------------------------------------
!    if mgrp_prscr_clds is active, call prscr_clds_amt to obtain the 
!    needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_mgroup_prescribed) then
          call prscr_clds_amt (is, ie, js, je, Cld_spec)

!----------------------------------------------------------------------
!    if gordon diagnostic clouds are active, call diag_clouds_amt to 
!    obtain the needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_diag_clouds) then
          call diag_clouds_amt (is, ie, js, je, lat, Atmos_input%pflux,&
                                Atmos_input%press, Rad_time, Cld_spec, &
                                Lsc_microphys)

!----------------------------------------------------------------------
!    if zetac clouds are active, call zetac_clouds_amt to 
!    obtain the needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_zetac_clouds) then
          if (present (r)) then
            Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
            Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
            Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
!           ierr = 0
!         else
!           call error_mesg ('cloud_spec_mod', &
!             ' must pass tracer array r when using zetac clouds', &
!                                                             FATAL)
          endif
          call zetac_clouds_amt (is, ie, js, je, z_half, z_full, &
                                 Surface%land, Atmos_input%phalf, &
                                 Atmos_input%deltaz, Cld_spec, &
                                 Lsc_microphys)

        endif ! (do_rh_clouds)
!--------------------------------------------------------------------
!    if klein prognostic clouds are active, call strat_clouds_amt to 
!    obtain the needed cloud specification variables.
!--------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds) then

!---------------------------------------------------------------------
!    if the gcm is being executed, call strat_cloud_avg to obtain the
!    appropriate (either instantaneous or time-averaged) values of
!    cloud water, cloud ice and cloud fraction. if the sa_gcm or the
!    standalone columns mode is being executed with the strat cloud
!    option, then values for the cloud water, cloud ice and when needed
!    cloud area have been input as optional arguments to this sub-
!    routine.
!---------------------------------------------------------------------
            if (present (r)) then
              Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
              Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
              Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
            endif

!---------------------------------------------------------------------
!    if the cloud input data is to be overriden, define the time slice
!    of data which is to be used. allocate storage for the cloud data.
!---------------------------------------------------------------------
          if (doing_data_override) then
            Data_time = Rad_time +    &
                          set_time (Rad_control%rad_time_step, 0)
 
!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    water data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current 
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qlnew', ql_proc,   &
                                  Data_time, override=override)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'ql => cloud_water not overridden successfully', FATAL)
            else
              Cld_spec%cloud_water(:,:,:) = ql_proc(is:ie,js:je,:)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    ice data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current 
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qinew', qi_proc,   &
                                Data_time, override=override)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'qi => cloud_ice   not overridden successfully', FATAL)
            else
              Cld_spec%cloud_ice(:,:,:) = qi_proc(is:ie,js:je,:)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    fraction data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qanew', qa_proc,   &
                                 Data_time, override=override)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
               'qa => cloud_area not overridden successfully', FATAL)
            else
              Cld_spec%cloud_area(:,:,:) = qa_proc(is:ie,js:je,:)
            endif
            ierr = 0
          endif ! (doing_override)

!---------------------------------------------------------------------
!    if values for the cloud variables have been successfully obtained,
!    call strat_clouds_amt to define the appropriate cloud specification
!    variables.
!---------------------------------------------------------------------
          if (ierr == 0) then
            call strat_clouds_amt (is, ie, js, je, Rad_time, &
                                   Atmos_input%pflux,  &
                                   Atmos_input%press,   &
                                   Atmos_input%cloudtemp, Surface%land,&
                                   Cld_spec, Lsc_microphys)

!----------------------------------------------------------------------
!    if ierr is non-zero, then cloud data was not successfully obtained.
!    if this is not the coldstart step, write an error message and 
!    stop execution.
!----------------------------------------------------------------------
          else 
            if (num_pts >= tot_pts) then
              call error_mesg ('cloud_spec_mod',  &
                     'no strat cloud data available; ierr /= 0', FATAL)

!----------------------------------------------------------------------
!    if this is the coldstart step, retain the input values corres-
!    ponding to no clouds, increment the points counter, and continue. 
!----------------------------------------------------------------------
            else
              num_pts = num_pts + size(Atmos_input%press,1)*   &
                                  size(Atmos_input%press,2)
            endif
          endif
        endif ! (do_strat_clouds)

!--------------------------------------------------------------------
!    since donner_deep_clouds may be active along with strat clouds, 
!    the associated properties are determined outside of the above loop.
!    these properties are placed in Cell_microphys and Meso_microphys.
!----------------------------------------------------------------------
        if (Cldrad_control%do_donner_deep_clouds) then
          call donner_deep_clouds_amt (is, ie, js, je,  &
                                       Cell_microphys, Meso_microphys)
        endif

!---------------------------------------------------------------------
!    obtain the microphysical properties (sizes and concentrations) if
!    a prescribed microphysics scheme is active. 
!---------------------------------------------------------------------
        if (Cldrad_control%do_presc_cld_microphys) then
          call microphys_presc_conc (is, ie, js, je,   &
                                     Atmos_input%clouddeltaz,   &
                                     Atmos_input%cloudtemp, &
                                     Cld_spec, Lsc_microphys)
        endif

!---------------------------------------------------------------------
!    call combine_cloud_properties to combine (if necessary) the cloud 
!    properties from multiple cloud types (large-scale, meso, cell) into
!    a single set for use by the radiation package. this is only needed
!    when microphysically-based properties are present, and when
!    either strat clouds and / or donner deep clouds is activated.
!---------------------------------------------------------------------
        if ( .not. Cldrad_control%do_specified_strat_clouds ) then
          if (Cldrad_control%do_sw_micro  .or.    &
              Cldrad_control%do_lw_micro) then
            if (Cldrad_control%do_strat_clouds .or.    &
                Cldrad_control%do_donner_deep_clouds) then
              call combine_cloud_properties (Lsc_microphys,    &
                                             Meso_microphys,   &
                                             Cell_microphys,   &
                                             Cld_spec)
            endif
          endif
        endif
      endif  !  (.not. do_no_clouds)

!--------------------------------------------------------------------
!    if microphysics is active and strat_clouds is not, define the water
!    paths (in units of kg / m**2).  if strat_clouds is active, these 
!    values will have already been defined. when microphysics is active,
!    define the effective sizes for the liquid and ice particles.
!--------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro .or.    &
          Cldrad_control%do_sw_micro)  then
        if (.not. Cldrad_control%do_strat_clouds .and.   &
            .not. Cldrad_control%do_zetac_clouds ) then
          Cld_spec%lwp = 1.0E-03*Lsc_microphys%conc_drop(:,:,:)* &
                         Atmos_input%clouddeltaz(:,:,:)
          Cld_spec%iwp = 1.0E-03*Lsc_microphys%conc_ice(:,:,:)*  &
                         Atmos_input%clouddeltaz(:,:,:)
        endif
        Cld_spec%reff_liq_micro = Lsc_microphys%size_drop
        Cld_spec%reff_ice_micro = Lsc_microphys%size_ice
      endif

!---------------------------------------------------------------------


end subroutine cloud_spec    



!######################################################################
! <SUBROUTINE NAME="cloud_spec_dealloc">
!  <OVERVIEW>
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec and the microphysics_type
!    structures Lsc_microphys, Meso_microphys and Cell_microphys. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec and the microphysics_type
!    structures Lsc_microphys, Meso_microphys and Cell_microphys. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec_dealloc (Cld_spec, Lsc_microphys, Meso_microphys,&
!                               Cell_microphys)
!  </TEMPLATE>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </INOUT>
! </SUBROUTINE>
! 
subroutine cloud_spec_dealloc (Cld_spec, Lsc_microphys, Meso_microphys,&
                               Cell_microphys)

!---------------------------------------------------------------------
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec and the microphysics_type
!    structures Lsc_microphys, Meso_microphys and Cell_microphys. 
!----------------------------------------------------------------------

type(cld_specification_type), intent(inout) :: Cld_spec
type(microphysics_type),      intent(inout) :: Lsc_microphys,   &
                                               Meso_microphys, &
                                               Cell_microphys

integer :: ier

!----------------------------------------------------------------------
!    deallocate the array elements of Cld_spec.
!----------------------------------------------------------------------
      deallocate (Cld_spec%camtsw         )  
      deallocate (Cld_spec%cmxolw         )
      deallocate (Cld_spec%crndlw         )
      deallocate (Cld_spec%ncldsw         )
      deallocate (Cld_spec%nmxolw         )
      deallocate (Cld_spec%nrndlw         )
      if (Cldrad_control%do_stochastic_clouds) then
        deallocate (Cld_spec%camtsw_band    )
        deallocate (Cld_spec%ncldsw_band    )
        deallocate (Cld_spec%cld_thickness_sw_band  )
        deallocate (Cld_spec%lwp_sw_band            )
        deallocate (Cld_spec%iwp_sw_band            )
        deallocate (Cld_spec%reff_liq_sw_band       )
        deallocate (Cld_spec%reff_ice_sw_band       )
      endif
      if (Cldrad_control%do_stochastic_clouds) then
        deallocate (Cld_spec%crndlw_band    )
        deallocate (Cld_spec%nrndlw_band    )
        deallocate (Cld_spec%cld_thickness_lw_band  )
        deallocate (Cld_spec%lwp_lw_band            )
        deallocate (Cld_spec%iwp_lw_band            )
        deallocate (Cld_spec%reff_liq_lw_band       )
        deallocate (Cld_spec%reff_ice_lw_band       )
      endif
      deallocate (Cld_spec%tau            )
      deallocate (Cld_spec%lwp            )
      deallocate (Cld_spec%iwp            )
      deallocate (Cld_spec%reff_liq       )
      deallocate (Cld_spec%reff_ice       )
      deallocate (Cld_spec%reff_liq_micro )
      deallocate (Cld_spec%reff_ice_micro )
      deallocate (Cld_spec%liq_frac       )
      deallocate (Cld_spec%cld_thickness  )
      deallocate (Cld_spec%hi_cloud       )
      deallocate (Cld_spec%mid_cloud      )
      deallocate (Cld_spec%low_cloud      )
      deallocate (Cld_spec%ice_cloud      )
      deallocate (Cld_spec%cloud_water    )
      deallocate (Cld_spec%cloud_ice      )
      deallocate (Cld_spec%cloud_area     )

!--------------------------------------------------------------------
!    deallocate the elements of Lsc_microphys.
!---------------------------------------------------------------------
      deallocate (Lsc_microphys%conc_drop   )
      deallocate (Lsc_microphys%conc_ice    )
      deallocate (Lsc_microphys%conc_rain   )
      deallocate (Lsc_microphys%conc_snow   )
      deallocate (Lsc_microphys%size_drop   )
      deallocate (Lsc_microphys%size_ice    )
      deallocate (Lsc_microphys%size_rain   )
      deallocate (Lsc_microphys%size_snow   )
      deallocate (Lsc_microphys%cldamt      )
      if (Cldrad_control%do_stochastic_clouds) then
        deallocate (Lsc_microphys%lw_stoch_conc_drop   , stat=ier)
        deallocate (Lsc_microphys%lw_stoch_conc_ice    , stat=ier)
        deallocate (Lsc_microphys%lw_stoch_size_drop   , stat=ier)
        deallocate (Lsc_microphys%lw_stoch_size_ice    , stat=ier)
        deallocate (Lsc_microphys%lw_stoch_cldamt      , stat=ier)

        deallocate (Lsc_microphys%sw_stoch_conc_drop   , stat=ier)
        deallocate (Lsc_microphys%sw_stoch_conc_ice    , stat=ier)
        deallocate (Lsc_microphys%sw_stoch_size_drop   , stat=ier)
        deallocate (Lsc_microphys%sw_stoch_size_ice    , stat=ier)
        deallocate (Lsc_microphys%sw_stoch_cldamt      , stat=ier)

        deallocate (Lsc_microphys%stoch_conc_drop   )
        deallocate (Lsc_microphys%stoch_conc_ice    )
        deallocate (Lsc_microphys%stoch_size_drop   )
        deallocate (Lsc_microphys%stoch_size_ice    )
        deallocate (Lsc_microphys%stoch_cldamt      )
      endif

!--------------------------------------------------------------------
!    deallocate the elements of Cell_microphys.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
        deallocate (Cell_microphys%conc_drop   )
        deallocate (Cell_microphys%conc_ice    )
        deallocate (Cell_microphys%conc_rain   )
        deallocate (Cell_microphys%conc_snow   )
        deallocate (Cell_microphys%size_drop   )
        deallocate (Cell_microphys%size_ice    )
        deallocate (Cell_microphys%size_rain   )
        deallocate (Cell_microphys%size_snow   )
        deallocate (Cell_microphys%cldamt      )

!--------------------------------------------------------------------
!    deallocate the elements of Meso_microphys.
!---------------------------------------------------------------------
        deallocate (Meso_microphys%conc_drop   )
        deallocate (Meso_microphys%conc_ice    )
        deallocate (Meso_microphys%conc_rain   )
        deallocate (Meso_microphys%conc_snow   )
        deallocate (Meso_microphys%size_drop   )
        deallocate (Meso_microphys%size_ice    )
        deallocate (Meso_microphys%size_rain   )
        deallocate (Meso_microphys%size_snow   )
        deallocate (Meso_microphys%cldamt      )
      endif

!---------------------------------------------------------------------


end subroutine cloud_spec_dealloc 



!#####################################################################

subroutine cloud_spec_end

!---------------------------------------------------------------------
!    cloud_spec_end is the destructor for cloud_spec_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    close the modules that were initialized by this module.
!--------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!-------------------------------------------------------------------
!    mgroup prescribed clouds.
!-------------------------------------------------------------------
        if (Cldrad_control%do_mgroup_prescribed) then
          call mgrp_prscr_clds_end 

!-------------------------------------------------------------------
!    rh-based diagnostic clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_rh_clouds) then
          call rh_based_clouds_end 

!-------------------------------------------------------------------
!    zonal or observed clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_zonal_clouds .or.  &
                 Cldrad_control%do_obs_clouds)  then
          call specified_clouds_W_end             

!-------------------------------------------------------------------
!    klein predicted clouds. if this option is active, donner_deep 
!    clouds may also be active. additionally, this may also have been
!    activated when running in standalone columns mode, in which case
!    standalone_clouds_end must be called.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_strat_clouds) then
          call strat_clouds_W_end
          if (Cldrad_control%do_donner_deep_clouds) then
            call donner_deep_clouds_W_end
          endif
          if (Cldrad_control%do_specified_strat_clouds .or. &
              Cldrad_control%do_specified_clouds ) then 
            call standalone_clouds_end
          endif

!-------------------------------------------------------------------
!    gordon diagnostic clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_diag_clouds) then
          call diag_clouds_W_end

!-------------------------------------------------------------------
!    zetac clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_zetac_clouds) then
          call zetac_clouds_W_end

!-------------------------------------------------------------------
!    donner deep convection clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_donner_deep_clouds) then
          call donner_deep_clouds_W_end

!-------------------------------------------------------------------
!    standalone specified clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_specified_clouds) then
          call standalone_clouds_end
        endif
      endif  ! (not do_no_clouds)

!--------------------------------------------------------------------
!    mark the module as no longer initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine cloud_spec_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!#####################################################################
! <SUBROUTINE NAME="initialize_cldamts">
!  <OVERVIEW>
!    initialize_cldamts allocates and initializes the array components 
!    of the structures used to specify the model cloud and microphysics
!    fields.
!  </OVERVIEW>
!  <DESCRIPTION>
!    initialize_cldamts allocates and initializes the array components 
!    of the structures used to specify the model cloud and microphysics
!    fields.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call initialize_cldamts (ix, jx, kx, Lsc_microphys,    &
!                               Meso_microphys, Cell_microphys, Cld_spec)
!  </TEMPLATE>
!  <IN NAME="ix, jx, kx" TYPE="integer">
!       ix             size of i dimension of physics window
!       jx             size of j dimension of physics window
!       kx             size of k dimension of physics window
!  </IN>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </INOUT>
! </SUBROUTINE>
! 
subroutine initialize_cldamts (ix, jx, kx, Lsc_microphys,    &
                               Meso_microphys, Cell_microphys, Cld_spec)

!---------------------------------------------------------------------
!    initialize_cldamts allocates and initializes the array components 
!    of the structures used to specify the model cloud and microphysics
!    fields.
!---------------------------------------------------------------------

integer,                      intent(in)     :: ix, jx, kx
type(microphysics_type),      intent(inout)  :: Lsc_microphys,   &
                                                Meso_microphys, &
                                                Cell_microphys
type(cld_specification_type), intent(inout)  :: Cld_spec

!----------------------------------------------------------------------
!    intent(in) variables:
! 
!       ix             size of i dimension of physics window
!       jx             size of j dimension of physics window
!       kx             size of k dimension of physics window
!
!    intent(inout) variables:
!
!       Lsc_microphys  microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!       Meso_microphys microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!       Cell_microphys microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!
!            the following elements are components of these  
!            microphysics_type variables which are allocated and
!            initialized here:
!
!            %conc_ice   ice particle concentration [ g / m**3 ]
!            %conc_drop  cloud droplet concentration [ g / m**3 ]
!            %conc_rain  rain drop concentration [ g / m**3 ]
!            %conc_snow  snow concentration [ g / m**3 ]
!            %size_ice   effective ice crystal diameter [ microns ]
!            %size_drop  effective cloud drop diameter [ microns ]
!            %size_rain  effective rain drop diameter [ microns ]
!            %size_snow  effective snow flake diameter [ microns ]
!            %cldamt     total cloud fraction (crystal + droplet)
!                        [ dimensionless ]
!            %lw_stoch_conc_ice
!                        ice particle concentration as a function of
!                        lw parameterization band [ g / m**3 ]
!            %lw_stoch_conc_drop
!                        cloud droplet concentration as a function of
!                        lw parameterization band [ g / m**3 ]
!            %lw_stoch_size_ice
!                        effective ice crystal diameter as a function
!                        of lw parameterization band [ microns ]
!            %lw_stoch_size_drop
!                        effective cloud drop diameter as a function of
!                        lw parameterization band [ microns ]
!            %lw_stoch_cldamt
!                        total cloud fraction (crystal + droplet) as a
!                        function of lw parameterization band
!                        [ dimensionless ]
!            %sw_stoch_conc_ice
!                        ice particle concentration as a function of
!                        sw parameterization band [ g / m**3 ]
!            %sw_stoch_conc_drop
!                        cloud droplet concentration as a function of
!                        sw parameterization band [ g / m**3 ]
!            %sw_stoch_size_ice
!                        effective ice crystal diameter as a function
!                        of sw parameterization band [ microns ]
!            %sw_stoch_size_drop
!                        effective cloud drop diameter as a function of
!                        sw parameterization band [ microns ]
!            %sw_stoch_cldamt
!                        total cloud fraction (crystal + droplet) as a
!                        function of sw parameterization band
!                        [ dimensionless ]
!
!       Cld_spec       variables on the model grid which define all or
!                      some of the following, dependent on the specific
!                      cloud parameterization: cloud optical paths, 
!                      particle sizes, cloud fractions, cloud thickness,
!                      number of clouds in a column, and /or cloud type 
!                      (high/mid/low, ice/liq or random/max overlap)
!                      [ cld_specification_type ]
!
!            the following elements are components of this  
!            cld_specification_type variable which is allocated and
!            initialized here:
!
!            %cmxolw         amount of maximally overlapped longwave 
!                            clouds [ dimensionless ]
!            %crndlw         amount of randomly overlapped longwave 
!                            clouds [ dimensionless ]
!            %nmxolw         number of maximally overlapped longwave 
!                            clouds in each grid column.
!            %nrndlw         number of maximally overlapped longwave 
!                            clouds in each grid column.
!            %camtsw         shortwave cloud amount. the sum of the max-
!                            imally overlapped and randomly overlapped 
!                            longwave cloud amounts. [ dimensionless ]
!            %ncldsw         number of shortwave clouds in each grid 
!                            column.
!            %camtsw_band    shortwave cloud amount. the sum of the max-
!                            imally overlapped and randomly overlapped
!                            longwave cloud amounts, differing with sw
!                            parameterization band. [ dimensionless ]
!            %crndlw_band    amount of randomly overlapped longwave
!                            clouds, differing with lw parameterization
!                            band [ dimensionless ]
!            %hi_cloud       logical mask for high clouds 
!            %mid_cloud      logical mask for middle clouds
!            %low_cloud      logical mask for low clouds
!            %ice_cloud      logical mask for ice clouds
!            %iwp            ice water path  [ kg / m**2 ]
!            %lwp            liquid water path [ kg / m**2 ]
!            %reff_liq       effective cloud drop radius  used with
!                            bulk cloud physics scheme [ microns ]
!            %reff_ice       effective ice crystal radius used with
!                            bulk cloud physics scheme [ microns ]
!            %reff_liq_micro effective cloud drop radius used with 
!                            microphysically based scheme [ microns ]
!            %reff_ice_micro effective ice crystal radius used with
!                            microphysically based scheme [ microns ]
!            %tau            extinction optical path  [ dimensionless ]
!            %liq_frac       fraction of cloud in a box which is liquid
!                            [ dimensionless ]
!            %cld_thickness  number of model layers contained in cloud  
!            %cloud_water    liquid cloud content [ kg liq / kg air ]
!            %cloud_ice      ice cloud content [ kg ice / kg air ]
!            %cloud_area     saturated volume fraction [ dimensionless ]
!
!---------------------------------------------------------------------

      integer  :: k, n

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters of the 
!    large-scale cloud scheme, including any precipitation fields and
!    the total cloud fraction. concentrations and fractions are init-
!    ialized to 0.0, and the effective sizes are set to small numbers to
!    avoid potential divides by zero.
!---------------------------------------------------------------------
      allocate (Lsc_microphys%conc_drop  (ix, jx, kx) )
      allocate (Lsc_microphys%conc_ice   (ix, jx, kx) )
      allocate (Lsc_microphys%conc_rain  (ix, jx, kx) )
      allocate (Lsc_microphys%conc_snow  (ix, jx, kx) )
      allocate (Lsc_microphys%size_drop  (ix, jx, kx) )
      allocate (Lsc_microphys%size_ice   (ix, jx, kx) )
      allocate (Lsc_microphys%size_rain  (ix, jx, kx) )
      allocate (Lsc_microphys%size_snow  (ix, jx, kx) )
      allocate (Lsc_microphys%cldamt     (ix, jx, kx) )
      Lsc_microphys%conc_drop(:,:,:) = 0.
      Lsc_microphys%conc_ice(:,:,:)  = 0.
      Lsc_microphys%conc_rain(:,:,:) = 0.
      Lsc_microphys%conc_snow(:,:,:) = 0.
      Lsc_microphys%size_drop(:,:,:) = 1.0e-20 
      Lsc_microphys%size_ice(:,:,:)  = 1.0e-20 
      Lsc_microphys%size_rain(:,:,:) = 1.0e-20        
      Lsc_microphys%size_snow(:,:,:) = 1.0e-20 
      Lsc_microphys%cldamt(:,:,:)    = 0.0
      if (Cldrad_control%do_stochastic_clouds) then
        allocate (Lsc_microphys%stoch_conc_drop  &
                                  (ix, jx, kx, Cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_conc_ice &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_size_drop  &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_size_ice   &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_cldamt  &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
  
        Lsc_microphys%lw_stoch_conc_drop => Lsc_microphys%stoch_conc_drop(:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_conc_drop => Lsc_microphys%stoch_conc_drop(:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_conc_ice  => Lsc_microphys%stoch_conc_ice (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_conc_ice  => Lsc_microphys%stoch_conc_ice (:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_size_drop => Lsc_microphys%stoch_size_drop(:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_size_drop => Lsc_microphys%stoch_size_drop(:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_size_ice  => Lsc_microphys%stoch_size_ice (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_size_ice  => Lsc_microphys%stoch_size_ice (:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_cldamt    => Lsc_microphys%stoch_cldamt   (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_cldamt    => Lsc_microphys%stoch_cldamt   (:, :, :, Cldrad_control%nlwcldb+1:)

       do n=1,Cldrad_control%nlwcldb + Solar_spect%nbands
        Lsc_microphys%stoch_conc_drop(:,:,:,n) = 0.
        Lsc_microphys%stoch_conc_ice(:,:,:,n)  = 0.
        Lsc_microphys%stoch_size_drop(:,:,:,n) = 1.0e-20
        Lsc_microphys%stoch_size_ice(:,:,:,n)  = 1.0e-20
        Lsc_microphys%stoch_cldamt(:,:,:,n)    = 0.0
       end do
      endif

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters of the 
!    meso-scale cloud scheme, including any precipitation fields and
!    the total cloud fraction. concentrations and fractions are init-
!    ialized to 0.0, and the effective sizes are set to small numbers to
!    avoid potential divides by zero.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
      allocate (Meso_microphys%conc_drop  (ix, jx, kx) )
      allocate (Meso_microphys%conc_ice   (ix, jx, kx) )
      allocate (Meso_microphys%conc_rain  (ix, jx, kx) )
      allocate (Meso_microphys%conc_snow  (ix, jx, kx) )
      allocate (Meso_microphys%size_drop  (ix, jx, kx) )
      allocate (Meso_microphys%size_ice   (ix, jx, kx) )
      allocate (Meso_microphys%size_rain  (ix, jx, kx) )
      allocate (Meso_microphys%size_snow  (ix, jx, kx) )
      allocate (Meso_microphys%cldamt     (ix, jx, kx) )
      Meso_microphys%conc_drop = 0.
      Meso_microphys%conc_ice  = 0.
      Meso_microphys%conc_rain = 0.
      Meso_microphys%conc_snow = 0.
      Meso_microphys%size_drop = 1.0e-20 
      Meso_microphys%size_ice  = 1.0e-20  
      Meso_microphys%size_rain = 1.0e-20                            
      Meso_microphys%size_snow = 1.0e-20 
      Meso_microphys%cldamt    = 0.0                               
      endif

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters of the 
!    cell-scale cloud scheme, including any precipitation fields and
!    the total cloud fraction. concentrations and fractions are init-
!    ialized to 0.0, and the effective sizes are set to small numbers to
!    avoid potential divides by zero.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
      allocate (Cell_microphys%conc_drop  (ix, jx, kx) )
      allocate (Cell_microphys%conc_ice   (ix, jx, kx) )
      allocate (Cell_microphys%conc_rain  (ix, jx, kx) )
      allocate (Cell_microphys%conc_snow  (ix, jx, kx) )
      allocate (Cell_microphys%size_drop  (ix, jx, kx) )
      allocate (Cell_microphys%size_ice   (ix, jx, kx) )
      allocate (Cell_microphys%size_rain  (ix, jx, kx) )
      allocate (Cell_microphys%size_snow  (ix, jx, kx) )
      allocate (Cell_microphys%cldamt     (ix, jx, kx) )
      Cell_microphys%conc_drop = 0.
      Cell_microphys%conc_ice  = 0.
      Cell_microphys%conc_rain = 0.
      Cell_microphys%conc_snow = 0.
      Cell_microphys%size_drop = 1.0e-20 
      Cell_microphys%size_ice  = 1.0e-20  
      Cell_microphys%size_rain = 1.0e-20                          
      Cell_microphys%size_snow = 1.0e-20 
      Cell_microphys%cldamt     = 0.
      endif

!---------------------------------------------------------------------
!    allocate arrays to hold the cloud fractions seen by the shortwave
!    and the random and maximum overlap fractions seen by the longwave
!    radiation, and then the number of each of these types of cloud in
!    each column. initialize the cloud fractions and number of clouds
!    to zero.
!---------------------------------------------------------------------
      allocate ( Cld_spec%camtsw (ix, jx, kx ) )
      allocate ( Cld_spec%cmxolw (ix, jx, kx ) )
      allocate ( Cld_spec%crndlw (ix, jx, kx ) )
      allocate ( Cld_spec%ncldsw (ix, jx     ) )
      allocate ( Cld_spec%nmxolw (ix, jx     ) )
      allocate ( Cld_spec%nrndlw (ix, jx     ) )
      Cld_spec%cmxolw(:,:,:) = 0.0E+00
      Cld_spec%crndlw(:,:,:) = 0.0E+00
      Cld_spec%camtsw(:,:,:) = 0.0E+00
      Cld_spec%nmxolw (:,:)  = 0
      Cld_spec%nrndlw (:,:)  = 0
      Cld_spec%ncldsw (:,:)  = 0
      if (Cldrad_control%do_stochastic_clouds) then
        allocate ( Cld_spec%camtsw_band    &
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate ( Cld_spec%ncldsw_band    &
                                 (ix, jx, Solar_spect%nbands) )
        allocate (Cld_spec%cld_thickness_sw_band & 
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%iwp_sw_band    &
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%lwp_sw_band   &
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%reff_liq_sw_band   &  
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%reff_ice_sw_band   &   
                                 (ix, jx, kx, Solar_spect%nbands) )
        do n=1,Solar_spect%nbands
        Cld_spec%camtsw_band(:,:,:,n) = 0.0E+00
        Cld_spec%ncldsw_band(:,:,n) = 0
        Cld_spec%cld_thickness_sw_band(:,:,:,n) = 0              
        Cld_spec%lwp_sw_band(:,:,:,n)   = 0.0
        Cld_spec%iwp_sw_band(:,:,:,n)   = 0.0
        Cld_spec%reff_liq_sw_band(:,:,:,n)      = 10.0
        Cld_spec%reff_ice_sw_band(:,:,:,n)      = 30.0
        end do
        allocate ( Cld_spec%crndlw_band    &
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate ( Cld_spec%nrndlw_band    &
                                 (ix, jx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%cld_thickness_lw_band & 
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%iwp_lw_band    &
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%lwp_lw_band   &
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%reff_liq_lw_band   &  
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%reff_ice_lw_band   &   
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        do n=1, Cldrad_control%nlwcldb
        Cld_spec%crndlw_band(:,:,:,n) = 0.0E+00
        Cld_spec%nrndlw_band(:,:,n) = 0
        Cld_spec%cld_thickness_lw_band(:,:,:,n) = 0              
        Cld_spec%lwp_lw_band(:,:,:,n)   = 0.0
        Cld_spec%iwp_lw_band(:,:,:,n)   = 0.0
        Cld_spec%reff_liq_lw_band(:,:,:,n)      = 10.0
        Cld_spec%reff_ice_lw_band (:,:,:,n)     = 30.0
        end do
      endif

!--------------------------------------------------------------------
!    allocate and initialize various arrays that are used by one or
!    another cloud scheme to specify the cloud locations and amounts.
!    initialization provides values consistent with the absence of
!    cloud, with the exception of the particle size fields which are
!    set to small, non-zero values.
!---------------------------------------------------------------------
      allocate (Cld_spec%hi_cloud       (ix, jx, kx) )
      allocate (Cld_spec%mid_cloud      (ix, jx, kx) )
      allocate (Cld_spec%low_cloud      (ix, jx, kx) )
      allocate (Cld_spec%ice_cloud      (ix, jx, kx) )
      allocate (Cld_spec%iwp            (ix, jx, kx) )
      allocate (Cld_spec%lwp            (ix, jx, kx) )
      allocate (Cld_spec%reff_liq       (ix, jx, kx) )
      allocate (Cld_spec%reff_ice       (ix, jx, kx) )
      allocate (Cld_spec%reff_liq_micro (ix, jx, kx) )
      allocate (Cld_spec%reff_ice_micro (ix, jx, kx) )
      allocate (Cld_spec%tau            (ix, jx, kx, num_slingo_bands) )
      allocate (Cld_spec%liq_frac       (ix, jx, kx) )
      allocate (Cld_spec%cld_thickness  (ix, jx, kx) )
      allocate (Cld_spec%cloud_water    (ix,jx,kx) )
      allocate (Cld_spec%cloud_ice      (ix,jx,kx)   )
      allocate (Cld_spec%cloud_area     (ix,jx,kx)  )

      Cld_spec%hi_cloud (:,:,:)     = .false.
      Cld_spec%mid_cloud(:,:,:)     = .false.
      Cld_spec%low_cloud(:,:,:)     = .false.
      Cld_spec%ice_cloud(:,:,:)     = .false.
      Cld_spec%lwp(:,:,:)    = 0.0
      Cld_spec%iwp(:,:,:)    = 0.0
      Cld_spec%reff_liq(:,:,:)      = 10.0
      Cld_spec%reff_ice(:,:,:)      = 30.0
      Cld_spec%reff_liq_micro(:,:,:) = 10.0
      Cld_spec%reff_ice_micro(:,:,:) = 30.0
      Cld_spec%liq_frac(:,:,:)      = 0.0
      Cld_spec%cld_thickness(:,:,:) = 0              
      Cld_spec%cloud_water(:,:,:)   = 0.
      Cld_spec%cloud_ice(:,:,:)     = 0.
      Cld_spec%cloud_area(:,:,:)    = 0.
      do n=1, num_slingo_bands
      Cld_spec%tau(:,:,:,n)  = 0.0
      end do

!---------------------------------------------------------------------



end subroutine initialize_cldamts              



!###################################################################
! <SUBROUTINE NAME="combine_cloud_properties">
!  <OVERVIEW>
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!  </OVERVIEW>
!  <DESCRIPTION>
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call combine_cloud_properties (Lsc_microphys, Meso_microphys,  &
!                                     Cell_microphys, Cld_spec)
!  </TEMPLATE>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </IN>
! </SUBROUTINE>
! 
subroutine combine_cloud_properties (Lsc_microphys, Meso_microphys,  &
                                     Cell_microphys, Cld_spec)

!----------------------------------------------------------------------
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!----------------------------------------------------------------------

type(microphysics_type),        intent(in)    :: Lsc_microphys, &
                                                 Meso_microphys, &
                                                 Cell_microphys
type(cld_specification_type), intent(inout)   :: Cld_spec

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       Lsc_microphys  microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!       Meso_microphys microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!       Cell_microphys microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!
!    intent(inout) variables:
!
!       Cld_spec       variables on the model grid which define all or
!                      some of the following, dependent on the specific
!                      cloud parameterization: cloud optical paths, 
!                      particle sizes, cloud fractions, cloud thickness,
!                      number of clouds in a column, and /or cloud type 
!                      (high/mid/low, ice/liq or random/max overlap)
!                      [ cld_specification_type ]
!
!---------------------------------------------------------------------

      integer   :: nb

!---------------------------------------------------------------------
!    total-cloud specification properties need be defined only when
!    strat_cloud and/or donner_deep clouds are active. 
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds .and.    &
          Cldrad_control%do_donner_deep_clouds) then

!----------------------------------------------------------------------
!    if strat_cloud and donner_deep are both active, define the random
!    overlap cloud fraction as the sum of the fractions of the large-
!    scale, meso-scale and cell-scale clouds.
!---------------------------------------------------------------------
        Cld_spec%crndlw = Lsc_microphys%cldamt +    &
                          Cell_microphys%cldamt + Meso_microphys%cldamt
        if (Cldrad_control%do_stochastic_clouds) then
          do nb=1, Cldrad_control%nlwcldb
            Cld_spec%crndlw_band(:,:,:,nb) =   &
                          Lsc_microphys%lw_stoch_cldamt(:,:,:,nb) +    &
                          Cell_microphys%cldamt(:,:,:) +     &
                          Meso_microphys%cldamt(:,:,:)
          end do
        endif
        if (Cldrad_control%do_stochastic_clouds) then
          do nb=1, Solar_spect%nbands
            Cld_spec%camtsw_band(:,:,:,nb) =   &
                         Lsc_microphys%sw_stoch_cldamt(:,:,:,nb) +    &
                         Cell_microphys%cldamt(:,:,:) +     &
                         Meso_microphys%cldamt(:,:,:)
          end do
        endif

!----------------------------------------------------------------------
!    if strat cloud is activated but donner_deep is not, define the 
!    total-cloud amount to be the large scale cloud amount.
!----------------------------------------------------------------------
      else if (Cldrad_control%do_strat_clouds) then
        Cld_spec%crndlw = Lsc_microphys%cldamt
        if (Cldrad_control%do_stochastic_clouds) then
          do nb=1, Cldrad_control%nlwcldb
            Cld_spec%crndlw_band(:,:,:,nb) =   &
                           Lsc_microphys%lw_stoch_cldamt(:,:,:,nb)
          end do
        endif
        if (Cldrad_control%do_stochastic_clouds) then
          do nb=1, Solar_spect%nbands
            Cld_spec%camtsw_band(:,:,:,nb) =   &
                           Lsc_microphys%sw_stoch_cldamt(:,:,:,nb)
          end do
        endif

!---------------------------------------------------------------------
!    if donner_deep is active but strat cloud is not, then the mesoscale
!    and cell-scale cloud amounts are combined to define the total
!    random-overlap cloud fraction.
!----------------------------------------------------------------------
      else if (Cldrad_control%do_donner_deep_clouds) then
        Cld_spec%crndlw = Cell_microphys%cldamt + Meso_microphys%cldamt
      endif

!---------------------------------------------------------------------
!    randomly-overlapped clouds are being assumed for donner_deep and 
!    strat cloud module clouds. set the max overlap cloud fraction to 
!    zero, be certain that the random overlap fraction is .le. 1. after
!    the summing of the component cloud fractions, and define the total
!    cloud fraction to be used by the sw code.
!---------------------------------------------------------------------
      Cld_spec%cmxolw = 0.0
      Cld_spec%crndlw = MIN (Cld_spec%crndlw, 1.00)
      if (Cldrad_control%do_stochastic_clouds) then
        Cld_spec%crndlw_band = MIN (Cld_spec%crndlw_band, 1.00)
        Cld_spec%camtsw_band = MIN (Cld_spec%camtsw_band, 1.00)
      endif
      Cld_spec%camtsw = Cld_spec%crndlw

!--------------------------------------------------------------------


end subroutine combine_cloud_properties 



!###################################################################
! <SUBROUTINE NAME="microphs_presc_conc">
!  <OVERVIEW>
!   Subroutine to determine water droplet and ice crystal based on
!   prescribed microphysics model.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine uses prescribed microphysics model to determine
!   concentrations of water droplets and ice crystals. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_presc_conc (is, ie, js, je, deltaz, temp,      &
!                                 Cld_spec, Lsc_microphys)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting indice of the x dimension in the physics domain
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   ending indice of the x dimension in the physics domain
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting indice of the y dimension in the physics domain
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   ending indice of the y dimension in the physics domain 
!  </IN>
!  <IN NAME="deltaz" TYPE="real">
!   Height of each pressure layers.
!  </IN>
!  <IN NAME="temp" TYPE="real">
!   Temperatures of pressure levels
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
! </SUBROUTINE>
!
subroutine microphys_presc_conc (is, ie, js, je, deltaz, temp,      &
                                 Cld_spec, Lsc_microphys)

!---------------------------------------------------------------------
!    microphys_presc_conc defines microphysical properties based on the
!    assumption of specified total water paths for high, middle and low 
!    clouds.
!---------------------------------------------------------------------

integer,                      intent(in)     :: is, ie, js, je
real, dimension(:,:,:),       intent(in)     :: deltaz, temp  
type(cld_specification_type), intent(in)     :: Cld_spec
type(microphysics_type),      intent(inout)  :: Lsc_microphys

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      deltaz         model vertical grid separation that is to be used
!                     for cloud calculations
!                     [meters]
!      temp           temperature at model levels (1:nlev) that is to
!                     be used in cloud calculations
!                     [ deg K ]
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution
!
!   intent(inout) variables:
!
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------

      real,    dimension(size(temp,1), size(temp,2), size(temp,3)) :: &
                                                       conc

      integer, dimension(size(temp,1), size(temp,2)) :: &
                                                       nhi_clouds, &
                                                       nmid_clouds, &
                                                       nlow_clouds

      integer  :: i,j,k

!--------------------------------------------------------------------
!  local variables:
!
!      conc             droplet concentration  [ g / m**3 ]
!      nhi_clouds       number of layers with high clouds
!      nmid_clouds      number of layers with middle clouds
!      nlow_clouds      number of layers with low clouds
!      i,j,k            do-loop indices
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!! RSH NOTE:
!
!    THE FOLLOWING treatment of diag_cloud_mod is here as an INITIAL 
!    IMPLEMENTATION to allow compilation and model execution, and 
!    provide "reasonable ?? " values.
! 
!    Code developed but NOT YET ADDED HERE reflects a later approach. 
!    That code is available under the fez release, and will be added to
!    the repository when upgrades to the cloud-radiation modules are 
!    completed.
!
!    obtain drop and ice size and concentrations, consistent with 
!    the diag_cloud scheme. As a test case, the following is a simple 
!    specification of constant concentration and size in all boxes 
!    defined as cloudy, attempting to come close to the prescribed 
!    values used for other cloud schemes. assume ice cld thickness 
!    = 2.0 km; then conc_ice=10.0E-03 => iwp = 20 g/m^2, similar to that
!    prescribed in microphys_presc_conc. assume water cld thickness 
!    = 3.5 km; then conc_drop = 20E-03 => lwp = 70 g / m^2, similar to 
!    that prescribed in microphys_presc_conc.  use sizes as used in 
!    microphys_presc_conc (50 and 20 microns). when done, radiative 
!    boundary fluxes are "similar" to non-microphysical results
!    for test case done here, and shows reasonable sensitivity to
!    variations in concentrations.
!    AGAIN, THIS IS AN INITIAL IMPLEMENTATION FOR TESTING ONLY !!!!
!
!!! RSH
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!---------------------------------------------------------------------

      if (Cldrad_control%do_diag_clouds) then
!---------------------------------------------------------------------
!    define concentrations and sizes of ice at those points which have
!    condensate and which were previously determined to be ice points
!    (Cld_spec%ice_cloud), and define concentrations and sizes of
!    liquid droplets at those points with condensate that had been 
!    determined to support liquid clouds.
!---------------------------------------------------------------------
        do k=1,size(Cld_spec%camtsw,3)
          do j=1,size(Cld_spec%camtsw,2)
            do i=1,size(Cld_spec%camtsw,1)
              if (Cld_spec%camtsw(i,j,k) > 0.0) then
                if (Cld_spec%ice_cloud(i,j,k)) then
                  Lsc_microphys%conc_ice(i,j,k) = 10.0E-03  
                  Lsc_microphys%size_ice(i,j,k) = 50.     
                else
                  Lsc_microphys%conc_drop(i,j,k) = 20.0E-03
                  Lsc_microphys%size_drop(i,j,k) = 20.
                endif
              endif
            end do
          end do
        end do
!----------------------------------------------------------------------
!    for the non-diag_cloud_mod cases, assume that the water path is 
!    preset at fixed values (lwpath_hi, _mid, _low) for "high", "mid", 
!    "low" clouds. the lwpath in each cloud layer within "hi", "mid" 
!    "low" pressure intervals is that lwpath_... divided by the number 
!    of clouds present in that pressure interval.
!----------------------------------------------------------------------
      else

!----------------------------------------------------------------------
!    define the number of high, middle, low clouds according to
!    Wetherald's criterion.
!----------------------------------------------------------------------
        do j=1,size(Cld_spec%camtsw,2)
          do i=1,size(Cld_spec%camtsw,1)
            nhi_clouds(i,j)  = 0
            nmid_clouds(i,j) = 0
            nlow_clouds(i,j) = 0
            do k=1,size(Cld_spec%camtsw,3)
              if (Cld_spec%hi_cloud(i,j,k)) &
                               nhi_clouds(i,j)  =  nhi_clouds(i,j)  + 1
              if (Cld_spec%mid_cloud(i,j,k)) &
                               nmid_clouds(i,j) =  nmid_clouds(i,j) + 1
              if (Cld_spec%low_cloud(i,j,k))  &
                               nlow_clouds(i,j) =  nlow_clouds(i,j) + 1
            end do
          end do
        end do

!----------------------------------------------------------------------
!    compute the water substance concentration in each layer 
!    (as water path / layer geometric path).
!----------------------------------------------------------------------
        conc(:,:,:) = 0.0E+00
        do j=1,size(Cld_spec%camtsw,2)
          do i=1,size(Cld_spec%camtsw,1)
            do k=1,size(Cld_spec%camtsw,3)
              if (Cld_spec%hi_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_hi/   &
                              (nhi_clouds(i,j)*deltaz(i,j,k))
              endif
              if (Cld_spec%mid_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_mid/    &
                              (nmid_clouds(i,j)*deltaz(i,j,k))
              endif
              if (Cld_spec%low_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_low    /                   &
                              (nlow_clouds(i,j)*deltaz(i,j,k))
              endif
            end do
          end do
        end do

!----------------------------------------------------------------------
!    split conc into conc_ice and conc_drop, depending on temperature
!    criterion (T < 273.16). assume that rain and / or snow are not
!    present.
!----------------------------------------------------------------------
        do k=1,size(Cld_spec%camtsw,3)
          do j=1,size(Cld_spec%camtsw,2)
            do i=1,size(Cld_spec%camtsw,1)
              if (temp(i,j,k) .LT. 273.16) then
                Lsc_microphys%conc_ice(i,j,k) = conc(i,j,k)
              else
                Lsc_microphys%conc_drop(i,j,k) = conc(i,j,k)
              endif
            end do
          end do
        end do

!----------------------------------------------------------------------
!    define sizes of microphysical species, using namelist values. note
!    that namelist drop and rain sizes are radii, so multiply by 2 to 
!    produce diameter, as desired for the %size_ arrays.
!----------------------------------------------------------------------
        Lsc_microphys%size_drop(:,:,:) = 2.0*wtr_cld_reff
        Lsc_microphys%size_rain(:,:,:) = 2.0*rain_reff
        Lsc_microphys%size_ice (:,:,:) = ice_cld_reff
      endif ! (do_diag_clouds)

!--------------------------------------------------------------------


end subroutine microphys_presc_conc


!#################################################################



                       end module cloud_spec_mod

