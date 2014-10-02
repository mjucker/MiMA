module atmos_model_mod
!<CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman  
!</CONTACT>
! <REVIEWER EMAIL="Zhi.Liang@noaa.gov">
!  Zhi Liang
! </REVIEWER>
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, mpp_error
use mpp_domains_mod,    only: domain2d
use fms_mod,            only: file_exist, error_mesg, field_size, FATAL, NOTE
use fms_mod,            only: close_file,  write_version_number, stdlog
use fms_mod,            only: read_data, write_data, clock_flag_default
use fms_mod,            only: open_restart_file, open_namelist_file, check_nml_error
use fms_io_mod,         only: get_restart_io_mode
use time_manager_mod,   only: time_type, operator(+), get_time
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: register_tracers
use diag_integral_mod,  only: diag_integral_init, diag_integral_end
use diag_integral_mod,  only: diag_integral_output
use atmosphere_mod,     only: atmosphere_up, atmosphere_down, atmosphere_init
use atmosphere_mod,     only: atmosphere_end, get_bottom_mass, get_bottom_wind
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_boundary, get_atmosphere_axes
use atmosphere_mod,     only: surf_diff_type

!-----------------------------------------------------------------------

implicit none
private

public update_atmos_model_down, update_atmos_model_up
public atmos_model_init, atmos_model_end, atmos_data_type
public land_ice_atmos_boundary_type, land_atmos_boundary_type
public ice_atmos_boundary_type
!-----------------------------------------------------------------------

!<PUBLICTYPE >
 type atmos_data_type
     type (domain2d)               :: domain             ! domain decomposition
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real, pointer, dimension(:)   :: glon_bnd => NULL() ! global longitude axis grid box boundaries in radians.
     real, pointer, dimension(:)   :: glat_bnd => NULL() ! global latitude axis grid box boundaries in radians.
     real, pointer, dimension(:)   :: lon_bnd  => NULL() ! local longitude axis grid box boundaries in radians.
     real, pointer, dimension(:)   :: lat_bnd  => NULL() ! local latitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: t_bot    => NULL() ! temperature at lowest model level
     real, pointer, dimension(:,:) :: q_bot    => NULL() ! specific humidity at lowest model level
     real, pointer, dimension(:,:) :: z_bot    => NULL() ! height above the surface for the lowest model level
     real, pointer, dimension(:,:) :: p_bot    => NULL() ! pressure at lowest model level
     real, pointer, dimension(:,:) :: u_bot    => NULL() ! zonal wind component at lowest model level
     real, pointer, dimension(:,:) :: v_bot    => NULL() ! meridional wind component at lowest model level
     real, pointer, dimension(:,:) :: p_surf   => NULL() ! surface pressure 
     real, pointer, dimension(:,:) :: gust     => NULL() ! gustiness factor
     real, pointer, dimension(:,:) :: coszen   => NULL() ! cosine of the zenith angle
     real, pointer, dimension(:,:) :: flux_sw  => NULL() ! net shortwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: flux_sw_dir            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_dif            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dir   =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dif   =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dir =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dif =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis_dir        =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis_dif        =>NULL()
     real, pointer, dimension(:,:) :: flux_lw  => NULL() ! net longwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: lprec    => NULL() ! mass of liquid precipitation since last time step (Kg/m2)
     real, pointer, dimension(:,:) :: fprec    => NULL() ! ass of frozen precipitation since last time step (Kg/m2)
     type (surf_diff_type)         :: Surf_diff          ! store data needed by the multi-step version of the diffusion algorithm
     type (time_type)              :: Time               ! current time
     type (time_type)              :: Time_step          ! atmospheric time step.
     type (time_type)              :: Time_init          ! reference time.
     integer, pointer              :: pelist(:) =>NULL() ! pelist where atmosphere is running.
     logical                       :: pe                 ! current pe.
 end type
!</PUBLICTYPE >

!<PUBLICTYPE >
type land_ice_atmos_boundary_type
   ! variables of this type are declared by coupler_main, allocated by flux_exchange_init.
!quantities going from land+ice to atmos
   real, dimension(:,:),   pointer :: t              =>NULL() ! surface temperature for radiation calculations
   real, dimension(:,:),   pointer :: albedo         =>NULL() ! surface albedo for radiation calculations
   real, dimension(:,:),   pointer :: albedo_vis_dir =>NULL()
   real, dimension(:,:),   pointer :: albedo_nir_dir =>NULL()
   real, dimension(:,:),   pointer :: albedo_vis_dif =>NULL()
   real, dimension(:,:),   pointer :: albedo_nir_dif =>NULL()
   real, dimension(:,:),   pointer :: land_frac      =>NULL() ! fraction amount of land in a grid box 
   real, dimension(:,:),   pointer :: dt_t           =>NULL() ! temperature tendency at the lowest level
   real, dimension(:,:),   pointer :: dt_q           =>NULL() ! specific humidity tendency at the lowest level
   real, dimension(:,:),   pointer :: u_flux         =>NULL() ! zonal wind stress
   real, dimension(:,:),   pointer :: v_flux         =>NULL() ! meridional wind stress
   real, dimension(:,:),   pointer :: dtaudu         =>NULL() ! derivative of zonal wind stress w.r.t. the lowest zonal level wind speed
   real, dimension(:,:),   pointer :: dtaudv         =>NULL() ! derivative of meridional wind stress w.r.t. the lowest meridional level wind speed
   real, dimension(:,:),   pointer :: u_star         =>NULL() ! friction velocity
   real, dimension(:,:),   pointer :: b_star         =>NULL() ! bouyancy scale
   real, dimension(:,:),   pointer :: q_star         =>NULL() ! moisture scale
   real, dimension(:,:),   pointer :: rough_mom      =>NULL() ! surface roughness (used for momentum)
   real, dimension(:,:,:), pointer :: data           =>NULL() !collective field for "named" fields above
   integer                         :: xtype                   !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
type :: land_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from land alone to atmos (none at present)
end type land_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
!quantities going from ice alone to atmos (none at present)
type :: ice_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from ice alone to atmos (none at present)
end type ice_atmos_boundary_type
!</PUBLICTYPE >

!Balaji
integer :: atmClock
!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmos_model.f90,v 12.0 2005/04/14 15:35:34 fms Exp $'
character(len=128) :: tagname = '$Name: lima $'

!-----------------------------------------------------------------------
character(len=80) :: restart_format = 'atmos_coupled_mod restart format 01'
!-----------------------------------------------------------------------
logical           :: do_netcdf_restart = .true.
logical           :: restart_tbot_qbot = .false.
namelist /atmos_model_nml/ do_netcdf_restart, restart_tbot_qbot  

contains

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_down">
!
! <OVERVIEW>
!   compute the atmospheric tendencies for dynamics, radiation, 
!   vertical diffusion of momentum, tracers, and heat/moisture.
! </OVERVIEW>
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down". 
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_model_down( Surface_boundary, Atmos )
!   </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

subroutine update_atmos_model_down( Surface_boundary, Atmos )
!
!-----------------------------------------------------------------------
  type(land_ice_atmos_boundary_type), intent(inout) :: Surface_boundary
  type (atmos_data_type), intent(inout) :: Atmos
                                      
!-----------------------------------------------------------------------
  call mpp_clock_begin(atmClock)

    call atmosphere_down (Atmos%Time, Surface_boundary%land_frac,        &
                          Surface_boundary%t,  Surface_boundary%albedo,  &
                          Surface_boundary%albedo_vis_dir,   &
                          Surface_boundary%albedo_nir_dir,   &
                          Surface_boundary%albedo_vis_dif,   &
                          Surface_boundary%albedo_nir_dif,   &
                          Surface_boundary%rough_mom,   &
                          Surface_boundary%u_star,      &
                          Surface_boundary%b_star,      &
                          Surface_boundary%q_star, &
                          Surface_boundary%dtaudu,      &
                          Surface_boundary%dtaudv,      &
                          Surface_boundary%u_flux,      &
                          Surface_boundary%v_flux,      &
                          Atmos%gust,                   &
                          Atmos%coszen,                 &
                          Atmos%flux_sw,                &
                          Atmos%flux_sw_dir,            &
                          Atmos%flux_sw_dif,            &
                          Atmos%flux_sw_down_vis_dir,   &
                          Atmos%flux_sw_down_vis_dif,   &
                          Atmos%flux_sw_down_total_dir, &
                          Atmos%flux_sw_down_total_dif, &
                          Atmos%flux_sw_vis,            &
                          Atmos%flux_sw_vis_dir,        &
                          Atmos%flux_sw_vis_dif,        &
                          Atmos%flux_lw,                &
                          Atmos%Surf_diff               )

!-----------------------------------------------------------------------

  call mpp_clock_end(atmClock)
 end subroutine update_atmos_model_down
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_up">
!
!-----------------------------------------------------------------------
! <OVERVIEW>
!   upward vertical diffusion of heat/moisture and moisture processes
! </OVERVIEW>

!<DESCRIPTION>
!   Called every time step as the atmospheric driver to finish the upward
!   sweep of the tridiagonal elimination for heat/moisture and compute the
!   convective and large-scale tendencies.  The atmospheric variables are
!   advanced one time step and tendencies set back to zero. 
!</DESCRIPTION>

! <TEMPLATE>
!     call  update_atmos_model_up( Surface_boundary, Atmos )
! </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

 subroutine update_atmos_model_up( Surface_boundary, Atmos )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
type (atmos_data_type), intent(inout) :: Atmos
                                      
!-----------------------------------------------------------------------
  call mpp_clock_begin(atmClock)


    Atmos%Surf_diff%delta_t = Surface_boundary%dt_t
    Atmos%Surf_diff%delta_q = Surface_boundary%dt_q

    call atmosphere_up (Atmos%Time,  Surface_boundary%land_frac, Atmos%Surf_diff, &
                        Atmos%lprec, Atmos%fprec, Atmos%gust)

!   --- advance time ---

    Atmos % Time = Atmos % Time + Atmos % Time_step


    call get_bottom_mass (Atmos % t_bot, Atmos % q_bot,  &
                          Atmos % p_bot, Atmos % z_bot,  &
                          Atmos % p_surf                 )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)


!------ global integrals ------

    call diag_integral_output (Atmos % Time)

!-----------------------------------------------------------------------
  call mpp_clock_end(atmClock)

end subroutine update_atmos_model_up
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!     This routine allocates storage and returns a variable of type
!     atmos_boundary_data_type, and also reads a namelist input and restart file. 
! </DESCRIPTION>

! <TEMPLATE>
!     call atmos_model_init (Atmos, Time_init, Time, Time_step)
! </TEMPLATE>

! <IN NAME="Time_init" TYPE="type(time_type)" >
!   The base (or initial) time of the experiment.
! </IN>

! <IN NAME="Time" TYPE="type(time_type)" >
!   The current time.
! </IN>

! <IN NAME="Time_step" TYPE="type(time_type)" >
!   The atmospheric model/physics time step.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)


type (atmos_data_type), intent(inout) :: Atmos
type (time_type), intent(in) :: Time_init, Time, Time_step

  integer :: unit, ntrace, ntprog, ntdiag, ntfamily, i, j
  integer :: mlon, mlat, nlon, nlat, sec, day, ipts, jpts, dt, dto
  real    :: r_ipts, r_jpts, r_dto
  character(len=80) :: control
  integer :: ierr, io, siz(4)
!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step

   if ( file_exist('input.nml')) then
      unit = open_namelist_file ( )
      ierr=1
      do while (ierr /= 0)
         read  (unit, nml=atmos_model_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'atmos_model_nml')
      enddo
 10     call close_file (unit)
   endif
   call get_restart_io_mode(do_netcdf_restart)

!-----------------------------------------------------------------------
! how many tracers have been registered?
!  (will print number below)
   call register_tracers ( MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfamily )
   if ( ntfamily > 0 ) call error_mesg ('atmos_model', 'ntfamily > 0', FATAL)

!-----------------------------------------------------------------------
!  ----- initialize atmospheric model -----

    call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                          Atmos%Surf_diff )
                           
!-----------------------------------------------------------------------
!---- allocate space ----

    call atmosphere_resolution (mlon, mlat, global=.true.)
    call atmosphere_resolution (nlon, nlat, global=.false.)
    call atmosphere_domain     (Atmos%domain)

    allocate ( Atmos % glon_bnd (mlon+1),    &
               Atmos % glat_bnd (mlat+1),    &
               Atmos %  lon_bnd (nlon+1),    &
               Atmos %  lat_bnd (nlat+1),    &
               Atmos % t_bot    (nlon,nlat), &
               Atmos % q_bot    (nlon,nlat), &
               Atmos % z_bot    (nlon,nlat), &
               Atmos % p_bot    (nlon,nlat), &
               Atmos % u_bot    (nlon,nlat), &
               Atmos % v_bot    (nlon,nlat), &
               Atmos % p_surf   (nlon,nlat), &
               Atmos % gust     (nlon,nlat), &
               Atmos % flux_sw  (nlon,nlat), &
               Atmos % flux_sw_dir (nlon,nlat), &
               Atmos % flux_sw_dif (nlon,nlat), &
               Atmos % flux_sw_down_vis_dir (nlon,nlat), &
               Atmos % flux_sw_down_vis_dif (nlon,nlat), &
               Atmos % flux_sw_down_total_dir (nlon,nlat), &
               Atmos % flux_sw_down_total_dif (nlon,nlat), &
               Atmos % flux_sw_vis (nlon,nlat), &
               Atmos % flux_sw_vis_dir (nlon,nlat), &
               Atmos % flux_sw_vis_dif(nlon,nlat), &
               Atmos % flux_lw  (nlon,nlat), &
               Atmos % coszen   (nlon,nlat), &
               Atmos % lprec    (nlon,nlat), &
               Atmos % fprec    (nlon,nlat)  )

    do j = 1, nlat
       do i = 1, nlon    
          Atmos % flux_sw(i,j)                 = 0.0
          Atmos % flux_lw(i,j)                 = 0.0    
          Atmos % flux_sw_dir (i,j)            = 0.0
          Atmos % flux_sw_dif (i,j)            = 0.0 
          Atmos % flux_sw_down_vis_dir (i,j)   = 0.0 
          Atmos % flux_sw_down_vis_dif (i,j)   = 0.0 
          Atmos % flux_sw_down_total_dir (i,j) = 0.0 
          Atmos % flux_sw_down_total_dif (i,j) = 0.0 
          Atmos % flux_sw_vis (i,j)            = 0.0 
          Atmos % flux_sw_vis_dir (i,j)        = 0.0 
          Atmos % flux_sw_vis_dif(i,j)         = 0.0 
          Atmos % coszen(i,j)                  = 0.0 
       enddo
    enddo
!-----------------------------------------------------------------------
!------ get initial state for dynamics -------

    call get_atmosphere_axes ( Atmos % axes )

    call atmosphere_boundary ( Atmos % glon_bnd, Atmos % glat_bnd, &
                               global=.true. )
    call atmosphere_boundary ( Atmos %  lon_bnd, Atmos %  lat_bnd, &
                               global=.false. )

    call get_bottom_mass (Atmos % t_bot, Atmos % q_bot,  &
                          Atmos % p_bot, Atmos % z_bot,  &
                          Atmos % p_surf                 )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)

!-----------------------------------------------------------------------
!---- print version number to logfile ----

   call write_version_number ( version, tagname )
!  write the namelist to a log file
   if( mpp_pe()==0 ) then
      unit = stdlog( )
      write (unit, nml=atmos_model_nml)
      call close_file (unit)
   endif

!  number of tracers
   if (mpp_pe() == mpp_root_pe()) then
        write (stdlog(), '(a,i3)') 'Number of tracers =', ntrace
        write (stdlog(), '(a,i3)') 'Number of prognostic tracers =', ntprog
        write (stdlog(), '(a,i3)') 'Number of diagnostic tracers =', ntdiag
   endif

!------ read initial state for several atmospheric fields ------

   if ( file_exist('INPUT/atmos_coupled.res.nc') ) then
       if(mpp_pe() == mpp_root_pe() ) call mpp_error ('atmos_model_mod', &
                   'Reading netCDF formatted restart file: INPUT/atmos_coupled.res.nc', NOTE)
       call read_data('INPUT/atmos_coupled.res.nc', 'glon_bnd', ipts, no_domain=.true.)
       call read_data('INPUT/atmos_coupled.res.nc', 'glat_bnd', jpts, no_domain=.true.)

       if (ipts /= mlon .or. jpts /= mlat) call error_mesg &
               ('coupled_atmos_init', 'incorrect resolution on restart file', FATAL)

       call read_data('INPUT/atmos_coupled.res.nc', 'dt', dto, no_domain=.true.)
       call read_data('INPUT/atmos_coupled.res.nc', 'lprec', Atmos % lprec, Atmos%domain)
       call read_data('INPUT/atmos_coupled.res.nc', 'fprec', Atmos % fprec, Atmos%domain)
       call read_data('INPUT/atmos_coupled.res.nc', 'gust', Atmos % gust, Atmos%domain)

       if (restart_tbot_qbot) then
          call read_data('INPUT/atmos_coupled.res.nc', 't_bot', Atmos%t_bot, Atmos%domain)
          call read_data('INPUT/atmos_coupled.res.nc', 'q_bot', Atmos%q_bot, Atmos%domain)
       endif 

       call get_time (Atmos % Time_step, sec, day)
       dt = sec + 86400*day  ! integer seconds
       if (dto /= dt) then
          Atmos % lprec = Atmos % lprec * real(dto)/real(dt)
          Atmos % fprec = Atmos % fprec * real(dto)/real(dt)
          if (mpp_pe() == mpp_root_pe()) write (stdlog(),50)
       endif
   else if (file_exist('INPUT/atmos_coupled.res')) then
          if(mpp_pe() == mpp_root_pe() ) call mpp_error ('atmos_model_mod', &
                   'Reading native formatted restart file: INPUT/atmos_coupled.res', NOTE)
          unit = open_restart_file ('INPUT/atmos_coupled.res', 'read')
          !--- check version number (format) of restart file ---
          read  (unit) control
          if (trim(control) /= trim(restart_format)) call error_mesg &
               ('coupled_atmos_init', 'invalid restart format', FATAL)
          !--- check resolution and time step ---
          read  (unit) ipts,jpts,dto
          if (ipts /= mlon .or. jpts /= mlat) call error_mesg &
               ('coupled_atmos_init', 'incorrect resolution on restart file', FATAL)

          !--- read data ---
          call read_data ( unit, Atmos % lprec )
          call read_data ( unit, Atmos % fprec )
          call read_data ( unit, Atmos % gust  )
          if (restart_tbot_qbot) then
             call read_data ( unit, Atmos % t_bot  )
             call read_data ( unit, Atmos % q_bot )
          endif
          call close_file (unit)

          !---- if the time step has changed then convert ----
       !        tendency to conserve mass of water
          call get_time (Atmos % Time_step, sec, day)
          dt = sec + 86400*day  ! integer seconds
          if (dto /= dt) then
             Atmos % lprec = Atmos % lprec * real(dto)/real(dt)
             Atmos % fprec = Atmos % fprec * real(dto)/real(dt)
             if (mpp_pe() == mpp_root_pe()) write (stdlog(),50)
 50         format (/,'The model time step changed .... &
                      &modifying precipitation tendencies')
          endif
   else
        Atmos % lprec = 0.0
        Atmos % fprec = 0.0
        Atmos % gust  = 1.0
   endif
  
!------ initialize global integral package ------

    call diag_integral_init (Atmos % Time_init, Atmos % Time,  &
                             Atmos % lon_bnd,   Atmos % lat_bnd)

!-----------------------------------------------------------------------
atmClock = mpp_clock_id( 'Atmosphere', flags=clock_flag_default, grain=CLOCK_COMPONENT )
end subroutine atmos_model_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_end (Atmos)

type (atmos_data_type), intent(inout) :: Atmos
integer :: unit, sec, day, dt
character(len=64) :: fname = 'RESTART/atmos_coupled.res.nc'
!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----
                                              
  call atmosphere_end (Atmos % Time)

!------ global integrals ------

  call diag_integral_end (Atmos % Time)

!---- compute integer time step (in seconds) ----
  call get_time (Atmos % Time_step, sec, day)
  dt = sec + 86400*day

!------ write several atmospheric fields ------
!        also resolution and time step

  if( do_netcdf_restart) then
     if(mpp_pe() == mpp_root_pe()) then
        call mpp_error ('atmos_model_mod', 'Writing netCDF formatted restart file.', NOTE)
     endif
     call write_data(fname, 'glon_bnd', real( size(Atmos%glon_bnd(:))-1), no_domain=.true. )
     call write_data(fname, 'glat_bnd', real( size(Atmos%glat_bnd(:))-1), no_domain=.true. )
     call write_data(fname, 'dt', real( dt), no_domain=.true. )
     call write_data(fname, 'lprec', Atmos%lprec, Atmos%domain)
     call write_data(fname, 'fprec', Atmos%fprec, Atmos%domain)
     call write_data(fname, 'gust', Atmos%gust, Atmos%domain)
     if(restart_tbot_qbot) then
        call write_data(fname, 't_bot', Atmos%t_bot, Atmos%domain)
        call write_data(fname, 'q_bot', Atmos%q_bot, Atmos%domain)
     endif
  else
     unit = open_restart_file ('RESTART/atmos_coupled.res', 'write')
     if (mpp_pe() == mpp_root_pe()) then
        write (unit) restart_format
        write (unit) size(Atmos%glon_bnd(:))-1, size(Atmos%glat_bnd(:))-1, dt
     endif
     call write_data ( unit, Atmos % lprec )
     call write_data ( unit, Atmos % fprec )
     call write_data ( unit, Atmos % gust  )
     if(restart_tbot_qbot) then
        call write_data ( unit, Atmos % t_bot  )
        call write_data ( unit, Atmos % q_bot  )
     endif
     call close_file (unit)
  endif

!-------- deallocate space --------

  deallocate ( Atmos % glon_bnd , &
               Atmos % glat_bnd , &
               Atmos %  lon_bnd , &
               Atmos %  lat_bnd , &
               Atmos % t_bot    , &
               Atmos % q_bot    , &
               Atmos % z_bot    , &
               Atmos % p_bot    , &
               Atmos % u_bot    , &
               Atmos % v_bot    , &
               Atmos % p_surf   , &
               Atmos % gust     , &
               Atmos % flux_sw  , &
               Atmos % flux_sw_dir  , &
               Atmos % flux_sw_dif  , &
               Atmos % flux_sw_down_vis_dir  , &
               Atmos % flux_sw_down_vis_dif  , &
               Atmos % flux_sw_down_total_dir  , &
               Atmos % flux_sw_down_total_dif  , &
               Atmos % flux_sw_vis  , &
               Atmos % flux_sw_vis_dir  , &
               Atmos % flux_sw_vis_dif  , &
               Atmos % flux_lw  , &
               Atmos % coszen   , &
               Atmos % lprec    , &
               Atmos % fprec      )

!-----------------------------------------------------------------------

end subroutine atmos_model_end
! </SUBROUTINE>

!#######################################################################

end module atmos_model_mod

