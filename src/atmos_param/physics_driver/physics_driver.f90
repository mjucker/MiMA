

module physics_driver_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!     Provides high level interfaces for calling the entire
!     FMS atmospheric physics package.
!
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
! </OVERVIEW>
! <DESCRIPTION>
!     This version of physics_driver_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Radiation, Rayleigh damping, gravity wave drag, vertical diffusion of
!     momentum and tracers, and the downward pass of vertical diffusion for
!     temperature and specific humidity are performed in the down routine.
!     The up routine finishes the vertical diffusion and computes moisture
!     related terms (convection,large-scale condensation, and precipitation).
! </DESCRIPTION>
! <DIAGFIELDS>
! </DIAGFIELDS>
! <DATASET NAME="physics_driver.res">
! native format restart file
! </DATASET>
!
! <DATASET NAME="physics_driver.res.nc">
! netcdf format restart file
! </DATASET>


! <INFO>

!   <REFERENCE>            </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!   </NOTE>
!   <FUTURE> Deal with conservation of total energy?              </FUTURE>

! </INFO>
!   shared modules:

use time_manager_mod,        only: time_type, get_time, operator (-), &
                                   time_manager_init
use field_manager_mod,       only: field_manager_init, MODEL_ATMOS
use tracer_manager_mod,      only: tracer_manager_init, &
                                   get_number_tracers

use atmos_tracer_driver_mod, only: atmos_tracer_driver_init,    &
                                   atmos_tracer_driver,  &
                                   atmos_tracer_driver_end
use fms_mod,                 only: mpp_clock_id, mpp_clock_begin,   &
                                   mpp_clock_end, CLOCK_MODULE_DRIVER, &
                                   MPP_CLOCK_SYNC,  fms_init,  &
                                   open_namelist_file, stdlog, &
                                   write_version_number, &
                                   file_exist, error_mesg, FATAL,   &
                                   WARNING, NOTE, check_nml_error, &
                                   open_restart_file, read_data, &
                                   close_file, mpp_pe, mpp_root_pe, &
                                   write_data, mpp_error, mpp_chksum
use fms_io_mod,              only: get_restart_io_mode

use constants_mod,           only: cp_air !mj for rrtmg

!    shared radiation package modules:

use rad_utilities_mod,       only: aerosol_type, radiative_gases_type, &
                                   rad_utilities_init, rad_output_type,&
                                   cld_specification_type,   &
                                   surface_type, &
                                   atmos_input_type, microphysics_type

!    component modules:

use  moist_processes_mod,    only: moist_processes,    &
                                   moist_processes_init,  &
                                   moist_processes_end,  &
                                   doing_strat

use vert_turb_driver_mod,    only: vert_turb_driver,  &
                                   vert_turb_driver_init,  &
                                   vert_turb_driver_end

use vert_diff_driver_mod,    only: vert_diff_driver_down,  &
                                   vert_diff_driver_up,    &
                                   vert_diff_driver_init,  &
                                   vert_diff_driver_end,   &
                                   surf_diff_type

use radiation_driver_mod,    only: radiation_driver_init,    &
                                   define_rad_times, define_surface,   &
                                   define_atmos_input_fields, &
                                   radiation_driver,  &
                                   atmos_input_dealloc,    &
                                   surface_dealloc, &
                                   radiation_driver_end
  
use cloud_spec_mod,          only: cloud_spec_init, cloud_spec, &
                                   cloud_spec_dealloc, cloud_spec_end

use aerosol_mod,             only: aerosol_init, aerosol_driver, &
                                   aerosol_dealloc, aerosol_end
 
use radiative_gases_mod,     only: radiative_gases_init,   &
                                   define_radiative_gases, &
                                   radiative_gases_dealloc, &
                                   radiative_gases_end

use damping_driver_mod,      only: damping_driver,      &
                                   damping_driver_init, &
                                   damping_driver_end

use grey_radiation_mod, only: grey_radiation_init, grey_radiation, grey_radiation_end

use local_heating_mod, only:  local_heating_init,local_heating

!mj: RRTM radiative scheme
use rrtmg_lw_init
use rrtmg_lw_rad
use rrtmg_sw_init
use rrtmg_sw_rad
use rrtm_radiation
use rrtm_vars     
!-----------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id: physics_driver.f90,v 12.0.6.2 2005/05/16 13:56:54 pjp Exp $'
character(len=128) :: tagname = '$Name:  $'


!---------------------------------------------------------------------
!-------  interfaces --------

public  physics_driver_init, physics_driver_down,   &
        physics_driver_up, physics_driver_end, &
        do_moist_in_phys_up, get_diff_t, &
        get_radturbten, zero_radturbten, &
        do_local_heating

private          &

!  called from physics_driver_init:
         read_restart_file, read_restart_nc,    &

!  called from physics_driver_down:
         check_args, &

!  called from check_args:
         check_dim

interface check_dim
     module procedure check_dim_2d, check_dim_3d, check_dim_4d
end interface


!---------------------------------------------------------------------
!------- namelist ------

logical :: do_moist_processes = .true.
                               ! call moist_processes routines
real    :: tau_diff = 3600.    ! time scale for smoothing diffusion 
                               ! coefficients
logical :: do_radiation = .false.
                               ! calculating radiative fluxes and
                               ! heating rates?

logical :: do_grey_radiation = .false.

logical :: do_rrtm_radiation = .true.

logical :: do_damping = .true.

logical :: do_local_heating = .false.

real    :: diff_min = 1.e-3    ! minimum value of a diffusion 
                               ! coefficient beneath which the
                               ! coefficient is reset to zero
logical :: diffusion_smooth = .true.
                               ! diffusion coefficients should be 
                               ! smoothed in time?
logical :: do_netcdf_restart = .true.              
! <NAMELIST NAME="physics_driver_nml">
!  <DATA NAME="do_netcdf_restart" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
! netcdf/native format restart file
!  </DATA>
!  <DATA NAME="do_radiation" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!calculating radiative fluxes and
! heating rates?
!  </DATA>
!  <DATA NAME="do_moist_processes" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!call moist_processes routines
!  </DATA>
!  <DATA NAME="tau_diff" UNITS="" TYPE="real" DIM="" DEFAULT="3600.">
!time scale for smoothing diffusion 
! coefficients
!  </DATA>
!  <DATA NAME="diff_min" UNITS="" TYPE="real" DIM="" DEFAULT="1.e-3">
!minimum value of a diffusion 
! coefficient beneath which the
! coefficient is reset to zero
!  </DATA>
!  <DATA NAME="diffusion_smooth" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!diffusion coefficients should be 
! smoothed in time?
!  </DATA>
! </NAMELIST>
!
namelist / physics_driver_nml / do_netcdf_restart, do_radiation, &
                                do_moist_processes, tau_diff,      &
                                diff_min, diffusion_smooth, &
                                do_grey_radiation, do_rrtm_radiation, &
                                do_damping, do_local_heating

!---------------------------------------------------------------------
!------- public data ------
! <DATA NAME="surf_diff_type" UNITS="" TYPE="surf_diff_type" DIM="" DEFAULT="">
! Defined in vert_diff_driver_mod, republished here. See vert_diff_mod for details.
! </DATA>

public  surf_diff_type   ! defined in  vert_diff_driver_mod, republished
                         ! here
 
!---------------------------------------------------------------------
!------- private data ------

!--------------------------------------------------------------------
! list of restart versions readable by this module:
!
! version 1: initial implementation 1/2003, contains diffusion coef-
!            ficient contribution from cu_mo_trans_mod. This variable
!            is generated in physics_driver_up (moist_processes) and
!            used on the next step in vert_diff_down, necessitating
!            its storage.
!
! version 2: adds pbltop as generated in vert_turb_driver_mod. This 
!            variable is then used on the next timestep by topo_drag
!            (called from damping_driver_mod), necessitating its 
!            storage.
!
! version 3: adds the diffusion coefficients which are passed to 
!            vert_diff_driver.  These diffusion are saved should
!            smoothing of vertical diffusion coefficients be turned
!            on.
!
! version 4: adds a logical variable, convect, which indicates whether
!            or not the grid column is convecting. This diagnostic is
!            needed by the entrain_module in vert_turb_driver.
!
! version 5: adds radturbten when strat_cloud_mod is active, adds 
!            lw_tendency when edt_mod or entrain_mod is active.
!
!---------------------------------------------------------------------
integer, dimension(5) :: restart_versions = (/ 1, 2, 3, 4, 5 /)

!--------------------------------------------------------------------
!    the following allocatable arrays are either used to hold physics 
!    data between timesteps when required, or hold physics data between
!    physics_down and physics_up.
!  
!    diff_cu_mo     contains contribution to difusion coefficient
!                   coming from cu_mo_trans_mod (called from 
!                   moist_processes in physics_driver_up) and then used 
!                   as input on the next time step to vert_diff_down 
!                   called in physics_driver_down.
!    diff_t         vertical diffusion coefficient for temperature
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    diff_m         vertical diffusion coefficient for momentum
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    radturbten     the sum of the radiational and turbulent heating,
!                   generated in both physics_driver_down (radiation)
!                   and physics_driver_up (turbulence) and then used
!                   in moist_processes
!    lw_tendency    longwave heating rate, generated in radiation and
!                   needed in vert_turb_driver when either edt_mod
!                   or entrain_mod is active. must be saved because
!                   radiation is not calculated on each step.
!    pbltop         top of boundary layer obtained from vert_turb_driver
!                   and then used on the next timestep in topo_drag_mod
!                   called from damping_driver_down        
!    convect        flag indicating whether convection is occurring in
!                   a grid column. generated in physics_driver_up and
!                   then used in vert_turb_driver called from 
!                   physics_driver_down on the next step.
!----------------------------------------------------------------------
real,    dimension(:,:,:), allocatable :: diff_cu_mo, diff_t, diff_m
real,    dimension(:,:,:), allocatable :: radturbten, lw_tendency
real,    dimension(:,:)  , allocatable :: pbltop     
logical, dimension(:,:)  , allocatable :: convect
   
!---------------------------------------------------------------------
!    internal timing clock variables:
!---------------------------------------------------------------------
integer :: radiation_clock, damping_clock, turb_clock,   &
           tracer_clock, diff_up_clock, diff_down_clock, &
           moist_processes_clock

!--------------------------------------------------------------------
!    miscellaneous control variables:
!---------------------------------------------------------------------
logical   :: do_check_args = .true.   ! argument dimensions should 
                                      ! be checked ?
logical   :: module_is_initialized = .false.
                                      ! module has been initialized ?
logical   :: doing_edt                ! edt_mod has been activated ?
logical   :: doing_entrain            ! entrain_mod has been activated ?
integer   :: nt                       ! total no. of tracers
integer   :: ntp                      ! total no. of prognostic tracers
!---------------------------------------------------------------------
!---------------------------------------------------------------------



                            contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="physics_driver_init">
!  <OVERVIEW>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_init (Time, lonb, latb, axes, pref, &
!                             trs, Surf_diff, phalf, mask, kbot  )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   reference prssure profiles
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   array of model latitudes at cell boundaries [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   array of model longitudes at cell boundaries [radians]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!  </IN>
!  <INOUT NAME="trs" TYPE="real">
!   atmospheric tracer fields
!  </INOUT>
!  <INOUT NAME="Surf_diff" TYPE="surf_diff_type">
!   surface diffusion derived type
!  </INOUT>
!  <IN NAME="phalf" TYPE="real">
!   pressure at model interface levels
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! <ERROR MSG="physics_driver_init must be called first" STATUS="FATAL">
! </ERROR>
! </SUBROUTINE>
!
subroutine physics_driver_init (Time, lonb, latb, axes, pref, &
                                trs, Surf_diff, phalf, mask, kbot, &
                                diffm, difft  )

!---------------------------------------------------------------------
!    physics_driver_init is the constructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type),         intent(in)              :: Time
real,dimension(:),       intent(in)              :: lonb, latb
integer,dimension(4),    intent(in)              :: axes
real,dimension(:,:),     intent(in)              :: pref
real,dimension(:,:,:,:), intent(inout)           :: trs
type(surf_diff_type),    intent(inout)           :: Surf_diff
real,dimension(:,:,:),   intent(in)              :: phalf
real,dimension(:,:,:),   intent(in),   optional  :: mask
integer,dimension(:,:),  intent(in),   optional  :: kbot
real, dimension(:,:,:),  intent(out),  optional  :: diffm, difft

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     Time       current time (time_type)
!     lonb       longitude of the grid box edges [ radians ]
!     latb       latitude of the grid box edges [ radians ]
!     axes       axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!     pref       two reference profiles of pressure at nlev+1 levels
!                pref(nlev+1,1)=101325. and pref(nlev+1,2)=81060.
!     phalf      pressure at model interface levels
!                [ Pa ]
!
!   intent(inout) variables:
!
!     trs        atmosperic tracer fields
!     Surf_diff  surface diffusion derived type variable
!
!   intent(in), optional variables:
!
!        mask    present when running eta vertical coordinate,
!                mask to remove points below ground
!        kbot    present when running eta vertical coordinate,
!                index of lowest model level above ground
!   
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(lonb(:))-1, size(latb(:))-1) :: sgsmtn
      character(len=64), dimension(:), pointer :: aerosol_names => NULL()
      character(len=64), dimension(:), pointer :: aerosol_family_names => NULL()
      integer          ::  id, jd, kd
      integer          ::  ierr, io, unit
      integer          ::  ndum

!---------------------------------------------------------------------
!  local variables:
!
!       sgsmtn        sgs orography obtained from mg_drag_mod;
!                     appears to not be currently used
!       aerosol_names names associated with the activated aerosols
!                     that will be seen by the radiation package
!       aerosol_family_names
!              names associated with the activated aerosol
!              families that will be seen by the radiation package
!       id,jd,kd      model dimensions on the processor  
!       ierr          error code
!       io            io status returned from an io call
!       unit          unit number used for an i/ operation

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that the modules used by this module that are not called 
!    later in this subroutine have already been initialized.
!---------------------------------------------------------------------
      call fms_init
 
!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=physics_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'physics_driver_nml')
        enddo
10      call close_file (unit)
      endif

      if(do_radiation .and. do_grey_radiation) &
        call error_mesg('physics_driver_init','do_radiation and do_grey_radiation cannot both be .true.',FATAL)
      if(do_radiation .and. do_rrtm_radiation) &
        call error_mesg('physics_driver_init','do_radiation and do_rrtm_radiation cannot both be .true.',FATAL)
      if(do_grey_radiation .and. do_rrtm_radiation) &
        call error_mesg('physics_driver_init','do_grey_radiation and do_rrtm_radiation cannot both be .true.',FATAL)

      if(do_radiation) call rad_utilities_init
      call time_manager_init
      call tracer_manager_init
      call field_manager_init (ndum)
      call get_restart_io_mode(do_netcdf_restart)

!--------------------------------------------------------------------
!    write version number and namelist to log file.
!--------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
               write(stdlog(), nml=physics_driver_nml)
 
!---------------------------------------------------------------------
!    define the model dimensions on the local processor.
!---------------------------------------------------------------------
      id = size(lonb(:))-1 
      jd = size(latb(:))-1 
      kd = size(trs,3)
      call get_number_tracers (MODEL_ATMOS, num_tracers=nt, &
                               num_prog=ntp)

!-----------------------------------------------------------------------
        call  moist_processes_init (id, jd, kd, lonb, latb, pref(:,1),&
                                    axes, Time)
     
!-----------------------------------------------------------------------
!    initialize damping_driver_mod.
!-----------------------------------------------------------------------
      if(do_damping) &
      call damping_driver_init (lonb, latb, pref(:,1), axes, Time, &
                                sgsmtn)

!-----------------------------------------------------------------------
!    initialize vert_turb_driver_mod.
!-----------------------------------------------------------------------
      call vert_turb_driver_init (lonb, latb, id, jd, kd, axes, Time, &
                                  doing_edt, doing_entrain)

!-----------------------------------------------------------------------
!    initialize vert_diff_driver_mod.
!-----------------------------------------------------------------------
      call vert_diff_driver_init (Surf_diff, id, jd, kd, axes, Time )

     if (do_radiation) then
!-----------------------------------------------------------------------
!    initialize cloud_spec_mod.
!-----------------------------------------------------------------------
      call cloud_spec_init (pref, lonb, latb, axes, Time)
 
!-----------------------------------------------------------------------
!    initialize aerosol_mod.     
!-----------------------------------------------------------------------
      call aerosol_init (lonb, latb, aerosol_names, aerosol_family_names)
 
!-----------------------------------------------------------------------
!    initialize radiative_gases_mod.
!-----------------------------------------------------------------------
      call radiative_gases_init (pref, latb, lonb)
 
!-----------------------------------------------------------------------
!    initialize radiation_driver_mod.
!-----------------------------------------------------------------------
      call radiation_driver_init (lonb, latb, pref, axes, time,  &
                                  aerosol_names, aerosol_family_names)

!---------------------------------------------------------------------
!    deallocate space for local pointers.
!---------------------------------------------------------------------
      deallocate (aerosol_names, aerosol_family_names)

      endif ! do_radiation

      if(do_grey_radiation) call grey_radiation_init(axes, Time)

      if(do_rrtm_radiation) then
         call rrtmg_lw_ini(cp_air)
         call rrtmg_sw_ini(cp_air)
         call rrtm_radiation_init(axes,Time,id*jd,kd,lonb,latb)
      endif
      
      if(do_local_heating) call local_heating_init(axes, Time)

!-----------------------------------------------------------------------
!    initialize atmos_tracer_driver_mod.
!-----------------------------------------------------------------------
      call atmos_tracer_driver_init (lonb, latb, trs, axes, time,  &
                                     phalf, mask)




!---------------------------------------------------------------------
!    initialize  various clocks used to time the physics components.
!---------------------------------------------------------------------
      radiation_clock       =       &
                mpp_clock_id( '   Physics_down: Radiation',    &
                   grain=CLOCK_MODULE_DRIVER, flags = MPP_CLOCK_SYNC )
      damping_clock         =     &
                mpp_clock_id( '   Physics_down: Damping',    &
                  grain=CLOCK_MODULE_DRIVER, flags = MPP_CLOCK_SYNC )
      turb_clock            =      &
                mpp_clock_id( '   Physics_down: Vert. Turb.', &
                  grain=CLOCK_MODULE_DRIVER, flags = MPP_CLOCK_SYNC )
      tracer_clock          =      &
                mpp_clock_id( '   Physics_down: Tracer',    &
                 grain=CLOCK_MODULE_DRIVER, flags = MPP_CLOCK_SYNC )
      diff_down_clock       =     &
                mpp_clock_id( '   Physics_down: Vert. Diff.',   &
                 grain=CLOCK_MODULE_DRIVER, flags = MPP_CLOCK_SYNC )
      diff_up_clock         =     &
                mpp_clock_id( '   Physics_up: Vert. Diff.',     &
                grain=CLOCK_MODULE_DRIVER, flags = MPP_CLOCK_SYNC )
      moist_processes_clock =      &
                mpp_clock_id( '   Physics_up: Moist Processes', &
                grain=CLOCK_MODULE_DRIVER, flags = MPP_CLOCK_SYNC )

!---------------------------------------------------------------------
!    allocate space for the module variables.
!---------------------------------------------------------------------
      allocate ( diff_t     (id, jd, kd) )
      allocate ( diff_m     (id, jd, kd) )
      allocate ( diff_cu_mo (id, jd, kd) )
      allocate ( pbltop     (id, jd) )
      allocate ( convect    (id, jd) )
      allocate ( radturbten (id, jd, kd))
      allocate ( lw_tendency(id, jd, kd))
       
!--------------------------------------------------------------------
!    call read_restart_file to obtain initial values for the module
!    variables.
!--------------------------------------------------------------------
      if(file_exist('INPUT/physics_driver.res.nc')) then
         call read_restart_nc
      else
         call read_restart_file ! Also handles initialization when no restart data exists
      endif

!---------------------------------------------------------------------
!    if desired, define variables to return diff_m and diff_t.
!---------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t
      endif
      if (present(diffm)) then
        diffm = diff_m
      endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------------



 end subroutine physics_driver_init


!######################################################################
! <SUBROUTINE NAME="physics_driver_down">
!  <OVERVIEW>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.    
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down (is, ie, js, je,                       &
!                                Time_prev, Time, Time_next,           &
!                                lat, lon, area,                       &
!                                p_half, p_full, z_half, z_full,       &
!                                u, v, t, q, r, um, vm, tm, qm, rm,    &
!                                frac_land, rough_mom,                 &
!                                albedo,    t_surf_rad,                &
!                                u_star,    b_star, q_star,            &
!                                dtau_du,  dtau_dv,  tau_x,  tau_y,    &
!                                udt, vdt, tdt, qdt, rdt,              &
!                                flux_sw,  flux_lw,  coszen,  gust,    &
!                                Surf_diff,                            &
!                                mask, kbot
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <INOUT NAME="rd" TYPE="real">
!   multiple 3d diagnostic tracer fields 
!  </INOUT>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <IN NAME="rough_mom" TYPE="real">
!   boundary layer roughness
!  </IN>
!  <IN NAME="albedo" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="t_surf_rad" TYPE="real">
!   surface radiative temperature
!  </IN>
!  <IN NAME="u_star" TYPE="real">
!   boundary layer wind speed (frictional speed)
!  </IN>
!  <IN NAME="b_star" TYPE="real">
!   ???
!  </IN>
!  <IN NAME="q_star" TYPE="real">
!   boundary layer specific humidity
!  </IN>
!  <IN NAME="dtau_du" TYPE="real">
!   derivative of zonal surface stress w.r.t zonal wind speed
!  </IN>
!  <IN NAME="dtau_dv" TYPE="real">
!   derivative of meridional surface stress w.r.t meridional wind speed
!  </IN>
!  <INOUT NAME="tau_x" TYPE="real">
!   boundary layer meridional component of wind shear
!  </INOUT>
!  <INOUT NAME="tau_y" TYPE="real">
!   boundary layer zonal component of wind shear
!  </INOUT>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="flux_sw" TYPE="real">
!   Shortwave flux from radiation package
!  </OUT>
!  <OUT NAME="flux_lw" TYPE="real">
!   Longwave flux from radiation package
!  </OUT>
!  <OUT NAME="coszen" TYPE="real">
!   cosine of zenith angle
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
!
!  <IN NAME="diff_cum_mom" TYPE="real">
!   OPTIONAL: present when do_moist_processes=.false.
!    cu_mo_trans diffusion coefficients, which are passed through to vert_diff_down.
!    Should not be present when do_moist_processes=.true., since these
!    values are passed out from moist_processes.
!  </IN>
!
!  <IN NAME="moist_convect" TYPE="real">
!   OPTIONAL: present when do_moist_processes=.false.
!    Should not be present when do_moist_processes=.true., since these
!    values are passed out from moist_processes.
!  </IN>
! </SUBROUTINE>
!
subroutine physics_driver_down (is, ie, js, je,                       &
                                Time_prev, Time, Time_next,           &
                                lat, lon, area,                       &
                                p_half, p_full, z_half, z_full,       &
                                u, v, t, q, r, um, vm, tm, qm, rm,    &
                                frac_land, rough_mom,                 &
                                albedo, albedo_vis_dir, albedo_nir_dir,&
                                albedo_vis_dif, albedo_nir_dif,       &
                                t_surf_rad,                           &
                                u_star,    b_star, q_star,            &
                                dtau_du, dtau_dv,  tau_x,  tau_y,     &
                                udt, vdt, tdt, qdt, rdt,              &
                                flux_sw,                              &
                                flux_sw_dir,                          &
                                flux_sw_dif,                          &
                                flux_sw_down_vis_dir,                 &
                                flux_sw_down_vis_dif,                 &
                                flux_sw_down_total_dir,               &
                                flux_sw_down_total_dif,               &
                                flux_sw_vis,                          &
                                flux_sw_vis_dir,                      &
                                flux_sw_vis_dif,                      &
                                flux_lw,  coszen,  gust,              &
                                Surf_diff,                            &
                                mask, kbot, diff_cum_mom,             &
                                moist_convect, diffm, difft  )

!---------------------------------------------------------------------
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!-----------------------------------------------------------------------

integer,                 intent(in)             :: is, ie, js, je
type(time_type),         intent(in)             :: Time_prev, Time,  &
                                                   Time_next
real,dimension(:,:),     intent(in)             :: lat, lon, area
real,dimension(:,:,:),   intent(in)             :: p_half, p_full,   &
                                                   z_half, z_full,   &
                                                   u , v , t , q ,   &
                                                   um, vm, tm, qm
real,dimension(:,:,:,:), intent(inout)          :: r
real,dimension(:,:,:,:), intent(inout)          :: rm
real,dimension(:,:),     intent(in)             :: frac_land,   &
                                                   rough_mom, &
                                                   albedo, t_surf_rad, &
                                                   albedo_vis_dir, albedo_nir_dir, &
                                                   albedo_vis_dif, albedo_nir_dif, &
                                                   u_star, b_star,    &
                                                   q_star, dtau_du, dtau_dv
real,dimension(:,:),     intent(inout)          :: tau_x,  tau_y
real,dimension(:,:,:),   intent(inout)          :: udt,vdt,tdt,qdt
real,dimension(:,:,:,:), intent(inout)          :: rdt
real,dimension(:,:),     intent(out)            :: flux_sw,  &
                                                   flux_sw_dir, &
                                                   flux_sw_dif, flux_lw,  &
                                                   coszen,  gust, &
                                                   flux_sw_down_vis_dir, &
                                                   flux_sw_down_vis_dif, &
                                                   flux_sw_down_total_dir, &
                                                   flux_sw_down_total_dif, &
                                                   flux_sw_vis, &
                                                   flux_sw_vis_dir, & 
                                                   flux_sw_vis_dif 
type(surf_diff_type),    intent(inout)          :: Surf_diff
real,dimension(:,:,:),   intent(in)   ,optional :: mask
integer, dimension(:,:), intent(in)   ,optional :: kbot
real,  dimension(:,:,:), intent(in)   ,optional :: diff_cum_mom
logical, dimension(:,:), intent(in)   ,optional :: moist_convect
real,  dimension(:,:,:), intent(out)  ,optional :: diffm, difft 

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      frac_land
!      rough_mom
!      albedo
!      albedo_vis_dir surface visible direct albedo [ dimensionless ]
!      albedo_nir_dir surface nir direct albedo [ dimensionless ]
!      albedo_vis_dif surface visible diffuse albedo [ dimensionless ]
!      albedo_nir_dif surface nir diffuse albedo [ dimensionless ]
!      t_surf_rad
!      u_star
!      b_star
!      q_star
!      dtau_du
!      dtau_dv
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!      rd             multiple 3d diagnostic tracer fields 
!                     [ unit / unit / sec ]
!      Surf_diff      surface_diffusion_type variable
!
!   intent(out) variables:
!
!      flux_sw
!      flux_sw_dir            net shortwave surface flux (down-up) [ w / m^2 ]
!      flux_sw_dif            net shortwave surface flux (down-up) [ w / m^2 ]
!      flux_sw_down_vis_dir   downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_down_vis_dif   downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_down_total_dir total downward shortwave surface flux [ w / m^2 ]
!      flux_sw_down_total_dif total downward shortwave surface flux [ w / m^2 ]
!      flux_sw_vis            net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_vis_dir        net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_vis_dif        net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_lw
!      coszen
!      gust
!
!   intent(in), optional variables:
!
!       mask        mask that designates which levels do not have data
!                   present (i.e., below ground); 0.=no data, 1.=data
!       kbot        lowest level which has data
!                   note:  both mask and kbot must be present together.
!
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      real, dimension(size(u,1),size(u,2),size(u,3)) :: diff_t_vert, &
                                                        diff_m_vert
      real, dimension(size(u,1),size(u,2))           :: z_pbl 
      type(aerosol_type)                             :: Aerosol
      type(cld_specification_type)                   :: Cld_spec
      type(radiative_gases_type)                     :: Rad_gases
      type(atmos_input_type)                         :: Atmos_input
      type(surface_type)                             :: Surface
      type(rad_output_type)                          :: Radiation
      type(time_type)                                :: Rad_time
      type(microphysics_type)                        :: Lsc_microphys, &
                                                        Meso_microphys,&
                                                        Cell_microphys
      integer          ::    sec, day
      real             ::    dt, alpha, dt2
      logical          ::    need_aerosols, need_clouds, need_gases,   &
                             need_basic

!---------------------------------------------------------------------
!   local variables:
!
!      diff_t_vert     vertical diffusion coefficient for temperature
!                      calculated on the current step
!      diff_m_vert     vertical diffusion coefficient for momentum   
!                      calculated on the current step
!      z_pbl           height of planetary boundary layer
!      Aerosol         aerosol_type variable describing the aerosol
!                      fields to be seen by the radiation package
!      Cld_spec        cld_specification_type variable describing the
!                      cloud field to be seen by the radiation package 
!      Rad_gases       radiative_gases_type variable describing the
!                      radiatively-active gas distribution to be seen 
!                      by the radiation package
!      Atmos_input     atmos_input_type variable describing the atmos-
!                      pheric state to be seen by the radiation package
!      Surface         surface_type variable describing the surface
!                      characteristics to be seen by the radiation 
!                      package
!      Radiation       rad_output_type variable containing the variables
!                      output from the radiation package, for passage
!                      to other modules
!      Rad_time        time at which the radiation calculation is to
!                      apply [ time_type ]
!      Lsc_microphys   microphysics_type variable containing the micro-
!                      physical characteristics of the large-scale
!                      clouds to be seen by the radiation package 
!      Meso_microphys  microphysics_type variable containing the micro-
!                      physical characteristics of the mesoscale
!                      clouds to be seen by the radiation package 
!      Cell_microphys  microphysics_type variable containing the micro-
!                      physical characteristics of the cell-scale
!                      clouds to be seen by the radiation package 
!      sec, day        second and day components of the time_type 
!                      variable
!      dt              model physics time step [ seconds ]
!      alpha           ratio of physics time step to diffusion-smoothing
!                      time scale
!      need_aerosols   need to obtain aerosol data on this time step
!                      to input to the radiation package ?
!      need_clouds     need to obtain cloud data on this time step
!                      to input to the radiation package ?
!      need_gases      need to obtain radiative gas data on this time 
!                      step to input to the radiation package ?
!      need_basic      need to obtain atmospheric state variables on
!                      this time step to input to the radiation package?
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
                         'module has not been initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check the size of the input arguments. this is only done on the
!    first call to physics_driver_down.
!---------------------------------------------------------------------
      if (do_check_args) call check_args  &
                   (lat, lon, area, p_half, p_full, z_half, z_full, &
                    u, v, t, q, r, um, vm, tm, qm, rm,              &
                    udt, vdt, tdt, qdt, rdt)

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next - Time_prev, sec, day)
      dt = real(sec + day*86400)

     if (do_radiation) then
!----------------------------------------------------------------------
!    prepare to calculate radiative forcings. obtain the valid time
!    at which the radiation calculation is to apply, the needed atmos-
!    pheric fields, and any needed inputs from other physics modules.
!---------------------------------------------------------------------
      call mpp_clock_begin ( radiation_clock )
 
!----------------------------------------------------------------------
!    call define_rad_times to obtain the time to be used in the rad-
!    iation calculation (Rad_time) and to determine which, if any, 
!    externally-supplied inputs to radiation_driver must be obtained on 
!    this timestep.  logical flags are returned indicating the need or 
!    lack of need for the aerosol fields, the cloud fields, the rad-
!    iative gas fields, and the basic atmospheric variable fields.
!----------------------------------------------------------------------
      call define_rad_times (Time, Time_next, Rad_time, &
                             need_aerosols, need_clouds, &
                             need_gases, need_basic)

!---------------------------------------------------------------------
!    call define_surface to define a surface_type variable containing
!    the surface albedoes and land fractions for each grid box. this
!    variable must be provided on all timesteps for use in generating
!    netcdf output.
!---------------------------------------------------------------------
      call define_surface (is, ie, js, je, albedo, albedo_vis_dir,  &
                           albedo_nir_dir, albedo_vis_dif, &
                           albedo_nir_dif, frac_land, Surface)

!---------------------------------------------------------------------
!    if the basic atmospheric input variables to the radiation package
!    are needed, pass the model pressure (p_full, p_half), temperature 
!    (t, t_surf_rad) and specific humidity (q) to subroutine
!    define_atmos_input_fields, which will put these fields and some
!    additional auxiliary fields into the form desired by the radiation
!    package and store them as components of the derived-type variable 
!    Atmos_input.
!---------------------------------------------------------------------
      if (need_basic) then
        call define_atmos_input_fields     &
                              (is, ie, js, je, p_full, p_half, t, q,  &
                               t_surf_rad, Atmos_input, kbot=kbot)
      endif

!---------------------------------------------------------------------
!    if the aerosol fields are needed as input to the radiation_package,
!    call aerosol_driver to access the aerosol data and place it into 
!    an aerosol_type derived-type variable Aerosol.
!---------------------------------------------------------------------
      if (need_aerosols) then
        call aerosol_driver (is, js, Rad_time, Atmos_input%pflux, &
                             Aerosol)
      endif
 
!---------------------------------------------------------------------
!    if the cloud fields are needed, call cloud_spec to retrieve bulk
!    cloud data and place it into a cld_specification_type derived-type 
!    variable Cld_spec and retrieve microphysical data which is returned
!    in microphysics_type variables Lsc_microphys, Meso_microphys and 
!    Cell_microphys, when applicable. 
!---------------------------------------------------------------------
      if (need_clouds) then
        if (present(kbot) ) then
          call cloud_spec (is, ie, js, je, lat,              &
                           z_half, z_full, Rad_time,   &
                           Atmos_input, Surface, Cld_spec,   &
                           Lsc_microphys, Meso_microphys,    &
                           Cell_microphys, r=r(:,:,:,1:ntp), &
                           kbot=kbot, mask=mask)
        else
          call cloud_spec (is, ie, js, je, lat,              &
                           z_half, z_full, Rad_time,   &
                           Atmos_input, Surface, Cld_spec,   &
                           Lsc_microphys, Meso_microphys,    &
                           Cell_microphys, r=r(:,:,:,1:ntp))
        endif
      endif

!---------------------------------------------------------------------
!    if the radiative gases are needed, call define_radiative_gases to 
!    obtain the values to be used for the radiatively-active gases and 
!    place them in radiative_gases_type derived-type variable Rad_gases.
!---------------------------------------------------------------------
      if (need_gases) then
        call define_radiative_gases (is, ie, js, je, Rad_time, lat, &
                                     Atmos_input, r, Time_next, Rad_gases)
      endif

!---------------------------------------------------------------------
!    allocate the components of a rad_output_type variable which will
!    be used to return the output from radiation_driver_mod that is
!    needed by other modules.
!---------------------------------------------------------------------
      allocate (Radiation%tdt_rad               (size(q,1), size(q,2),size(q,3)))
      allocate (Radiation%flux_sw_surf          (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_surf_dir      (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_surf_dif      (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_down_vis_dir  (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_down_vis_dif  (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_down_total_dir(size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_down_total_dif(size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_vis           (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_vis_dir       (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_sw_vis_dif       (size(q,1), size(q,2)          ))
      allocate (Radiation%flux_lw_surf          (size(q,1), size(q,2)          ))
      allocate (Radiation%coszen_angle          (size(q,1), size(q,2)          ))
      allocate (Radiation%tdtlw                 (size(q,1), size(q,2),size(q,3)))

!--------------------------------------------------------------------
!    call  radiation_driver to perform the radiation calculation.
!--------------------------------------------------------------------
      !call radiation_driver (is, ie, js, je, Time, Time_next, lat,  &
      !                       lon, Surface, Atmos_input, Aerosol, &
      !                       Cld_spec, Rad_gases, Lsc_microphys, &
      !                       Meso_microphys, Cell_microphys,&
      !                       Radiation=Radiation, mask=mask, kbot=kbot)

!-------------------------------------------------------------------
!    process the variables returned from radiation_driver_mod. the 
!    radiative heating rate is added to the accumulated physics heating
!    rate (tdt). net surface lw and sw fluxes and the cosine of the 
!    zenith angle are placed in locations where they can be exported
!    for use in other component models. the lw heating rate is stored
!    in a module variable for potential use in other physics modules.
!    the radiative heating rate is also added to a variable which is
!    accumulating the radiative and turbulent heating rates, and which
!    is needed by strat_cloud_mod.
!-------------------------------------------------------------------
      tdt     = tdt + Radiation%tdt_rad(:,:,:)
      flux_sw = Radiation%flux_sw_surf
      flux_sw_dir            = Radiation%flux_sw_surf_dir
      flux_sw_dif            = Radiation%flux_sw_surf_dif
      flux_sw_down_vis_dir   = Radiation%flux_sw_down_vis_dir
      flux_sw_down_vis_dif   = Radiation%flux_sw_down_vis_dif
      flux_sw_down_total_dir = Radiation%flux_sw_down_total_dir
      flux_sw_down_total_dif = Radiation%flux_sw_down_total_dif
      flux_sw_vis            = Radiation%flux_sw_vis
      flux_sw_vis_dir        = Radiation%flux_sw_vis_dir
      flux_sw_vis_dif        = Radiation%flux_sw_vis_dif
      flux_lw = Radiation%flux_lw_surf
      coszen  = Radiation%coszen_angle
      lw_tendency(is:ie,js:je,:) = Radiation%tdtlw(:,:,:)
      radturbten (is:ie,js:je,:) = radturbten(is:ie,js:je,:) + &
                                   Radiation%tdt_rad(:,:,:)

!--------------------------------------------------------------------
!    deallocate the arrays used to return the radiation_driver_mod 
!    output.
!--------------------------------------------------------------------
      deallocate ( Radiation%tdt_rad      )
      deallocate ( Radiation%flux_sw_surf )
      deallocate ( Radiation%flux_sw_surf_dir )
      deallocate ( Radiation%flux_sw_surf_dif )
      deallocate ( Radiation%flux_sw_down_vis_dir )
      deallocate ( Radiation%flux_sw_down_vis_dif )
      deallocate ( Radiation%flux_sw_down_total_dir )
      deallocate ( Radiation%flux_sw_down_total_dif )
      deallocate ( Radiation%flux_sw_vis )
      deallocate ( Radiation%flux_sw_vis_dir )
      deallocate ( Radiation%flux_sw_vis_dif )
      deallocate ( Radiation%flux_lw_surf )
      deallocate ( Radiation%coszen_angle )
      deallocate ( Radiation%tdtlw        )
 
!---------------------------------------------------------------------
!    call routines to deallocate the components of the derived type 
!    arrays input to radiation_driver.
!---------------------------------------------------------------------
      if (need_gases) then
        call radiative_gases_dealloc (Rad_gases)
      endif
      if (need_clouds) then
        call cloud_spec_dealloc (Cld_spec, Lsc_microphys,   &
                                 Meso_microphys, Cell_microphys)
      endif
      if (need_aerosols) then
        call aerosol_dealloc (Aerosol)
      endif
      if (need_basic) then
        call atmos_input_dealloc (Atmos_input)
      endif
      call surface_dealloc (Surface)
      call mpp_clock_end ( radiation_clock )
      else
        flux_sw = 0.0
        flux_sw_dir = 0.0
        flux_sw_dif = 0.0
        flux_sw_down_vis_dir = 0.0
        flux_sw_down_vis_dif = 0.0
        flux_sw_down_total_dir = 0.0
        flux_sw_down_total_dif = 0.0
        flux_sw_vis = 0.0
        flux_sw_vis_dir = 0.0
        flux_sw_vis_dif = 0.0
        flux_lw = 0.0
        coszen  = 0.0
        lw_tendency(is:ie,js:je,:) = 0.0
      endif ! do_radiation

      if(do_grey_radiation) then
         call grey_radiation(is, js, Time_next, lat, lon, p_half, albedo, t_surf_rad, t, tdt, flux_sw, flux_lw)
         coszen = 1.0
      endif
      if(do_rrtm_radiation) then
         !need t at half grid
         call interp_temp(z_full,z_half,t_surf_rad,t)
         call run_rrtmg(is,js,Time,lat,lon,p_full,p_half,albedo,q,t,t_surf_rad,tdt,coszen,flux_sw,flux_lw)
      endif
!----------------------------------------------------------------------
!    artificial local heating if required
!----------------------------------------------------------------------
      if(do_local_heating) then
        call local_heating(is,js,Time,lon,lat,p_full,tdt)
      endif

!----------------------------------------------------------------------
!    call damping_driver to calculate the various model dampings that
!    are desired. 
!----------------------------------------------------------------------
      z_pbl(:,:) = pbltop(is:ie,js:je) 
      if(do_damping) then
        call mpp_clock_begin ( damping_clock )
        call damping_driver (is, js, lat, Time_next, dt,           &
                             p_full, p_half, z_full, z_half,          &
                             um, vm, tm, qm, rm(:,:,:,1:ntp), &
                             udt, vdt, tdt, qdt, rdt,&
                             z_pbl , mask=mask, kbot=kbot)
       call mpp_clock_end ( damping_clock )
     endif

!---------------------------------------------------------------------
!    If moist_processes is not called in physics_driver_down then values
!    of convect must be passed in via the optional argument "moist_convect".
!---------------------------------------------------------------------
      if(.not.do_moist_processes) then
        if(present(moist_convect)) then
          convect(is:ie,js:je) = moist_convect
        else
          call error_mesg('physics_driver_down', &
          'moist_convect be present when do_moist_processes=.false.',FATAL) 
        endif
      endif

!---------------------------------------------------------------------
!    call vert_turb_driver to calculate diffusion coefficients. save
!    the planetary boundary layer height on return.
!---------------------------------------------------------------------
      call mpp_clock_begin ( turb_clock )
      call vert_turb_driver (is, js, Time, Time_next, dt,            &
                             lw_tendency(is:ie,js:je,:), frac_land,  &
                             p_half, p_full, z_half, z_full, u_star, &
                             b_star, q_star, rough_mom, lat,         &
                             convect(is:ie,js:je),                   &
                             u, v, t, q, r(:,:,:,1:ntp), um, vm,     &
                             tm, qm, rm(:,:,:,1:ntp),                &
                             udt, vdt, tdt, qdt, rdt,                &
                             diff_t_vert, diff_m_vert, gust, z_pbl,  &
                             mask=mask, kbot=kbot             )
     call mpp_clock_end ( turb_clock )
     pbltop(is:ie,js:je) = z_pbl(:,:)

     
!-----------------------------------------------------------------------
!    process any tracer fields.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( tracer_clock )
      call atmos_tracer_driver (is, ie, js, je, Time, lon, lat,  &
                                frac_land, p_half, p_full, r, u, v, t, &
                                q, u_star, rdt, rm, dt, z_half,   &
                                z_full, t_surf_rad, albedo, coszen,  &
                                Time_next, kbot)
      call mpp_clock_end ( tracer_clock )

!-----------------------------------------------------------------------
!    If moist_processes is not called in physics_driver_down then values
!    of the cu_mo_trans diffusion coefficients must be passed in via
!    the optional argument "diff_cum_mom".
!-----------------------------------------------------------------------
      if(.not.do_moist_processes) then
        if(present(diff_cum_mom)) then
          diff_cu_mo(is:ie,js:je,:) = diff_cum_mom          
        else
          call error_mesg('physics_driver_down', &
          'diff_cum_mom must be present when do_moist_processes=.false.',FATAL)
        endif
      endif

!-----------------------------------------------------------------------
!    optionally use an implicit calculation of the vertical diffusion 
!    coefficients.
!
!    the vertical diffusion coefficients are solved using an implicit
!    solution to the following equation:
!
!    dK/dt   = - ( K - K_cur) / tau_diff
!
!    where K         = diffusion coefficient
!          K_cur     = diffusion coefficient diagnosed from current 
!                      time steps' state
!          tau_diff  = time scale for adjustment
!
!    in the code below alpha = dt / tau_diff
!---------------------------------------------------------------------
      if (diffusion_smooth) then
        call get_time (Time_next - Time, sec, day)
        dt2 = real(sec + day*86400)
        alpha = dt2/tau_diff
        diff_m(is:ie,js:je,:) = (diff_m(is:ie,js:je,:) +       &
                                 alpha*(diff_m_vert(:,:,:) +  &
                                 diff_cu_mo(is:ie,js:je,:)) )/&
                                 (1. + alpha)
        where (diff_m(is:ie,js:je,:) < diff_min)
          diff_m(is:ie,js:je,:) = 0.0
        end where
        diff_t(is:ie,js:je,:) = (diff_t(is:ie,js:je,:) +      &
                                 alpha*diff_t_vert(:,:,:) )/  &
                                 (1. + alpha)
        where (diff_t(is:ie,js:je,:) < diff_min)
          diff_t(is:ie,js:je,:) = 0.0
        end where
      else
        diff_t(is:ie,js:je,:) = diff_t_vert
        diff_m(is:ie,js:je,:) = diff_m_vert + diff_cu_mo(is:ie, js:je,:)
      end if

!-----------------------------------------------------------------------
!    call vert_diff_driver_down to calculate the first pass atmos-
!    pheric vertical diffusion.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( diff_down_clock )
      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) - tdt(:,:,:)
      call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
                                  p_full, z_full,   &
                                  diff_m(is:ie,js:je,:),         &
                                  diff_t(is:ie,js:je,:),         &
                                  um ,vm ,tm ,qm ,rm(:,:,:,1:ntp), &
                                  dtau_du, dtau_dv, tau_x, tau_y,  &
                                  udt, vdt, tdt, qdt, rdt,       &
                                  Surf_diff,                     &
                                  mask=mask, kbot=kbot           )

!---------------------------------------------------------------------
!    if desired, return diff_m and diff_t to calling routine.
!-----------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t(is:ie,js:je,:)
      endif
      if (present(diffm)) then
        diffm = diff_m(is:ie,js:je,:)
      endif
      
     call mpp_clock_end ( diff_down_clock )

 end subroutine physics_driver_down



!#######################################################################
! <SUBROUTINE NAME="physics_driver_up">
!  <OVERVIEW>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_up (is, ie, js, je,                    &
!                               Time_prev, Time, Time_next,        &
!                               lat, lon, area,                    &
!                               p_half, p_full, z_half, z_full,    & 
!                               omega,                             &
!                               u, v, t, q, r, um, vm, tm, qm, rm, &
!                               frac_land,                         &
!                               udt, vdt, tdt, qdt, rdt,           &
!                               Surf_diff,                         &
!                               lprec,   fprec, gust,              &
!                               mask, kbot    )
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="omega" TYPE="real">
!   Veritical pressure tendency
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="lprec" TYPE="real">
!  </OUT>
!  <OUT NAME="fprec" TYPE="real">
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
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
 subroutine physics_driver_up (is, ie, js, je,                    &
                               Time_prev, Time, Time_next,        &
                               lat, lon, area,                    &
                               p_half, p_full, z_half, z_full,    & 
                               omega,                             &
                               u, v, t, q, r, um, vm, tm, qm, rm, &
                               frac_land,                         &
                               udt, vdt, tdt, qdt, rdt,           &
                               Surf_diff,                         &
                               lprec,   fprec, gust,              &
                               mask, kbot                         )

!----------------------------------------------------------------------
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!---------------------------------------------------------------------

integer,                intent(in)             :: is, ie, js, je
type(time_type),        intent(in)             :: Time_prev, Time,   &
                                                  Time_next
real,dimension(:,:),    intent(in)             :: lat, lon, area
real,dimension(:,:,:),  intent(in)             :: p_half, p_full,   &
                                                  omega,  &
                                                  z_half, z_full,     &
                                                  u , v , t , q ,    &
                                                  um, vm, tm, qm
real,dimension(:,:,:,:),intent(in)             :: r,rm
real,dimension(:,:),    intent(in)             :: frac_land
real,dimension(:,:,:),  intent(inout)          :: udt,vdt,tdt,qdt
real,dimension(:,:,:,:),intent(inout)          :: rdt
type(surf_diff_type),   intent(inout)          :: Surf_diff
real,dimension(:,:),    intent(out)            :: lprec, fprec
real,dimension(:,:),    intent(inout)          :: gust
real,dimension(:,:,:),  intent(in),   optional :: mask
integer,dimension(:,:), intent(in),   optional :: kbot

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      omega
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      frac_land
!      rough_mom
!      albedo
!      t_surf_rad
!      u_star
!      b_star
!      q_star
!      dtau_du
!      dtau_dv
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!      Surf_diff      surface_diffusion_type variable
!      gust
!
!   intent(out) variables:
!
!      lprec     
!      fprec       
!
!   intent(in), optional variables:
!
!       mask        mask that designates which levels do not have data
!                   present (i.e., below ground); 0.=no data, 1.=data
!       kbot        lowest level which has data
!                   note:  both mask and kbot must be present together.
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!   local variables:

      real, dimension(size(u,1), size(u,2), size(u,3)) :: diff_cu_mo_loc
      real, dimension(size(u,1), size(u,2))            :: gust_cv
      integer :: sec, day
      real    :: dt
   
!---------------------------------------------------------------------
!   local variables:
!
!        diff_cu_mo_loc   diffusion coefficient contribution due to 
!                         cumulus momentum transport
!        gust_cv
!        sec, day         second and day components of the time_type 
!                         variable
!        dt               physics time step [ seconds ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
             'module has not been initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next-Time_prev, sec, day)
      dt = real(sec+day*86400)

!------------------------------------------------------------------
!    call vert_diff_driver_up to complete the vertical diffusion
!    calculation.
!------------------------------------------------------------------
      call mpp_clock_begin ( diff_up_clock )
! XXX df's version of vert_diff_driver_up requires t for one of his new diagnostic fields
        call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                  Surf_diff, tdt, qdt, mask=mask,  &
                                  kbot=kbot, t=t)
      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + tdt(:,:,:)
      call mpp_clock_end ( diff_up_clock )

!-----------------------------------------------------------------------
!    if the fms integration path is being followed, call moist processes
!    to compute moist physics, including convection and processes 
!    involving condenstion.
!-----------------------------------------------------------------------
      if (do_moist_processes) then
        call mpp_clock_begin ( moist_processes_clock )
        call moist_processes (is, ie, js, je, Time_next, dt, frac_land, &
                              p_half, p_full, z_half, z_full, omega,    &
                              diff_t(is:ie,js:je,:),                    &
                              radturbten(is:ie,js:je,:),                &
                              t, q, r, u, v, tm, qm, rm, um, vm,        &
                              tdt, qdt, rdt, udt, vdt, diff_cu_mo_loc , &
                              convect(is:ie,js:je), lprec, fprec,       &
                              gust_cv, area, lat, mask=mask, kbot=kbot)
        call mpp_clock_end ( moist_processes_clock )
        diff_cu_mo(is:ie, js:je,:) = diff_cu_mo_loc(:,:,:)
        radturbten(is:ie,js:je,:) = 0.0

!---------------------------------------------------------------------
!    add the convective gustiness effect to that previously obtained 
!    from non-convective parameterizations.
!---------------------------------------------------------------------
        gust = sqrt( gust*gust + gust_cv*gust_cv)
      endif ! do_moist_processes

!-----------------------------------------------------------------------


 end subroutine physics_driver_up


!#######################################################################
! <SUBROUTINE NAME="physics_driver_end">
!  <OVERVIEW>
!   physics_driver_end is the destructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_end is the destructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_end (Time)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
! </SUBROUTINE>
!
subroutine physics_driver_end (Time)

!---------------------------------------------------------------------
!    physics_driver_end is the destructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      Time      current time [ time_type(days, seconds) ]
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variable:

     integer :: unit    ! unit number for restart file
     character(len=64)  :: fname='RESTART/physics_driver.res.nc'
     real,dimension(size(convect, 1), size(convect, 2))    :: r_convect
!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
              'module has not been initialized', FATAL)
      endif

      if(do_netcdf_restart) then
         if (mpp_pe() == mpp_root_pe() ) then
            call error_mesg('physics_driver_mod', 'Writing netCDF formatted restart file: RESTART/physics_driver.res.nc', NOTE)
         endif
         call write_data(fname, 'vers', real(restart_versions(size(restart_versions(:)))), no_domain =.true. )
         if(doing_strat()) then 
            call write_data(fname, 'doing_strat', 1.0, no_domain=.true.)
         else
            call write_data(fname, 'doing_strat', 0.0, no_domain=.true.)
         endif
         if(doing_edt) then 
            call write_data(fname, 'doing_edt', 1.0, no_domain=.true.)
         else
            call write_data(fname, 'doing_edt', 0.0, no_domain=.true.)
         endif
         if(doing_entrain) then 
            call write_data(fname, 'doing_entrain', 1.0, no_domain=.true.)
         else
            call write_data(fname, 'doing_entrain', 0.0, no_domain=.true.)
         endif
         !--------------------------------------------------------------------
         !    write out the data fields that are relevant for this experiment.
         !--------------------------------------------------------------------
         call write_data (fname, 'diff_cu_mo', diff_cu_mo)
         call write_data (fname, 'pbltop', pbltop)
         call write_data (fname, 'diff_t', diff_t)
         call write_data (fname, 'diff_m', diff_m)
         r_convect = 0.
         where(convect)
            r_convect = 1.0
         end where
         call write_data (fname, 'convect', r_convect)
         if (doing_strat()) then
            call write_data (fname, 'radturbten', radturbten)
         endif
         if (doing_edt .or. doing_entrain) then
            call write_data (fname, 'lw_tendency', lw_tendency)
         endif
      else
         if (mpp_pe() == mpp_root_pe() ) then
            call error_mesg('physics_driver_mod', 'Writing native formatted restart file.', NOTE)
         endif
         !---------------------------------------------------------------------
         !    open a unit for the restart file.
         !---------------------------------------------------------------------
         unit = open_restart_file ('RESTART/physics_driver.res', 'write')
         
         !---------------------------------------------------------------------
         !    write the header records, indicating restart version number and
         !    which variables fields are present.
         !---------------------------------------------------------------------
         if (mpp_pe() == mpp_root_pe() ) then
            write (unit) restart_versions(size(restart_versions(:)))
            write (unit) doing_strat(), doing_edt, doing_entrain
         endif
         
         !--------------------------------------------------------------------
         !    write out the data fields that are relevant for this experiment.
         !--------------------------------------------------------------------
         call write_data (unit, diff_cu_mo)
         call write_data (unit, pbltop    )
         call write_data (unit, diff_t)
         call write_data (unit, diff_m)
         call write_data (unit, convect)
         if (doing_strat()) then
            call write_data (unit, radturbten)
         endif
         if (doing_edt .or. doing_entrain) then
            call write_data (unit, lw_tendency)
         endif
         
         !--------------------------------------------------------------------
         !    close the restart file unit.
         !--------------------------------------------------------------------
         call close_file (unit)
      endif
!--------------------------------------------------------------------
!    call the destructor routines for those modules who were initial-
!    ized from this module.
!--------------------------------------------------------------------
      call vert_turb_driver_end
      call vert_diff_driver_end
      if (do_radiation) then
        call radiation_driver_end
        call radiative_gases_end
        call cloud_spec_end
        call aerosol_end
      endif
      if(do_grey_radiation) call grey_radiation_end
      if(do_rrtm_radiation) call rrtm_radiation_end
      call moist_processes_end
      call atmos_tracer_driver_end
      if(do_damping) call damping_driver_end

!---------------------------------------------------------------------
!    deallocate the module variables.
!---------------------------------------------------------------------
      deallocate (diff_cu_mo, diff_t, diff_m, pbltop, convect,   &
                  radturbten, lw_tendency)
 
!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.


!-----------------------------------------------------------------------

 end subroutine physics_driver_end



!#######################################################################
! <FUNCTION NAME="do_moist_in_phys_up">
!  <OVERVIEW>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </OVERVIEW>
!  <DESCRIPTION>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </DESCRIPTION>
!  <TEMPLATE>
!   logical = do_moist_in_phys_up()
!  </TEMPLATE>
! </FUNCTION>
!
function do_moist_in_phys_up()

!--------------------------------------------------------------------
!    do_moist_in_phys_up returns the value of do_moist_processes
!----------------------------------------------------------------------

logical :: do_moist_in_phys_up

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('do_moist_in_phys_up',  &
              'module has not been initialized', FATAL)
      endif
 
!-------------------------------------------------------------------
!    define output variable.
!-------------------------------------------------------------------
      do_moist_in_phys_up = do_moist_processes

 
end function do_moist_in_phys_up

!#####################################################################
! <FUNCTION NAME="get_diff_t">
!  <OVERVIEW>
!    returns the values of array diff_t
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array diff_t
!  </DESCRIPTION>
!  <TEMPLATE>
!   diff_t(:,:,:) = get_diff_t()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_diff_t() result(diff_t_out)
real, dimension(size(diff_t,1),size(diff_t,2),size(diff_t,3)) :: diff_t_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_diff_t','module has not been initialized', FATAL)
  endif

  diff_t_out = diff_t

end function get_diff_t

!#####################################################################
! <FUNCTION NAME="get_radturbten">
!  <OVERVIEW>
!    returns the values of array radturbten
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array radturbten
!  </DESCRIPTION>
!  <TEMPLATE>
!   radturbten(:,:,:) = get_radturbten()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_radturbten() result(radturbten_out)
real, dimension(size(radturbten,1),size(radturbten,2),size(radturbten,3)) :: radturbten_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_radturbten','module has not been initialized', FATAL)
  endif

  radturbten_out = radturbten

end function get_radturbten
!#####################################################################
! <SUBROUTINE NAME="zero_radturbten">
!  <OVERVIEW>
!    sets all values of array radturbten to zero
!  </OVERVIEW>
!  <DESCRIPTION>
!    sets all values of array radturbten to zero
!  </DESCRIPTION>
!  <TEMPLATE>
!   call zero_radturbten()
!  </TEMPLATE>
! </SUBROUTINE>
!
!#####################################################################
subroutine zero_radturbten()

  if ( .not. module_is_initialized) then
    call error_mesg ('zero_radturbten','module has not been initialized', FATAL)
  endif

  radturbten = 0.0

end subroutine zero_radturbten
!#####################################################################



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
               
     
!#####################################################################
! <SUBROUTINE NAME="read_restart_file">
!  <OVERVIEW>
!    read_restart_file will read the physics_driver.res file and process
!    its contents. if no restart data can be found, the module variables
!    are initialized to flag values.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_restart_file will read the physics_driver.res file and process
!    its contents. if no restart data can be found, the module variables
!    are initialized to flag values.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_restart_file
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine read_restart_file                                     

!---------------------------------------------------------------------
!    read_restart_file will read the physics_driver.res file and process
!    its contents. if no restart data can be found, the module variables
!    are initialized to flag values.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      integer  :: ierr, io, unit
      integer  :: vers, vers2
      character(len=8) :: chvers
      logical  :: was_doing_strat, was_doing_edt, was_doing_entrain
      logical  :: success = .false.

!--------------------------------------------------------------------
!   local variables:
!
!      ierr              error code
!      io                error status returned from i/o operation
!      unit              io unit number for reading restart file
!      vers              restart version number if that is contained in 
!                        file; otherwise the first word of first data 
!                        record of file
!      vers2             second word of first data record of file
!      was_doing_strat   logical indicating if strat_cloud_mod was 
!                        active in job which wrote restart file
!      was_doing_edt     logical indicating if edt_mod was active
!                        in job which wrote restart file
!      was_doing_entrain logical indicating if entrain_mod was active
!                        in job which wrote restart file
!      success           logical indicating that restart data has been
!                        processed
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    obtain values for radturbten, either from physics_driver.res, if
!    reading a newer version of the file which contains it, or from 
!    strat_cloud.res when an older version of physics_driver.res is
!    being read.
!--------------------------------------------------------------------
      if(mpp_pe() == mpp_root_pe()) call mpp_error ('physics_driver_mod', &
            'Reading native formatted restart file.', NOTE)
      if (file_exist('INPUT/physics_driver.res')) then
        unit = open_restart_file ('INPUT/physics_driver.res', 'read')

!--------------------------------------------------------------------
!    read restart file version number.
!--------------------------------------------------------------------
        read (unit) vers
        if ( .not. any(vers ==restart_versions) ) then
          write (chvers, '(i4)') vers
          call error_mesg ('physics_driver_mod', &
            'restart version ' //chvers// ' cannot be read by this'//&
                                              'module version', FATAL)
        endif

!--------------------------------------------------------------------
!    starting with v5,  logicals are written indicating which variables
!    are present.
!--------------------------------------------------------------------
        if (vers >= 5 ) then
          read (unit) was_doing_strat, was_doing_edt, was_doing_entrain
        endif

!---------------------------------------------------------------------
!    read the contribution to diffusion coefficient from cumulus
!    momentum transport.
!---------------------------------------------------------------------
        call read_data (unit, diff_cu_mo)

!---------------------------------------------------------------------
!    pbl top is present in file versions 2 and up. if not present,
!    set a flag.
!---------------------------------------------------------------------
        if (vers >= 2) then
          call read_data (unit, pbltop)
        else
          pbltop     = -999.0
        endif

!---------------------------------------------------------------------
!    the temperature and momentum diffusion coefficients are present
!    beginning with v3. if not prsent, set to 0.0.
!---------------------------------------------------------------------
        if (vers >= 3) then
          call read_data (unit, diff_t)
          call read_data (unit, diff_m)
        else
          diff_t = 0.0
          diff_m = 0.0
        end if 

!---------------------------------------------------------------------
!    a flag indicating columns in which convection is occurring is
!    present beginning with v4. if not present, set it to .false.
!---------------------------------------------------------------------
        if (vers >= 4) then
          call read_data (unit, convect)
        else
          convect = .false.
        end if 

!---------------------------------------------------------------------
!    radturbten may be present in versions 5 onward, if strat_cloud_mod
!    was active in the job writing the .res file.
!---------------------------------------------------------------------
        if (vers >= 5) then

!--------------------------------------------------------------------
!    if radturbten was written, read it.
!--------------------------------------------------------------------
          if (was_doing_strat) then
            call read_data (unit, radturbten)

!---------------------------------------------------------------------
!    if strat_cloud_mod was not active in the job which wrote the 
!    restart file but it is active in the current job, initialize
!    radturbten to 0.0 and put a message in the output file.  
!---------------------------------------------------------------------
          else
            if (doing_strat()) then
              radturbten = 0.0
              call error_mesg ('physics_driver_mod', &
              ' initializing radturbten to 0.0, since it not present'//&
                            ' in physics_driver.res file', NOTE)
            endif
          endif

!--------------------------------------------------------------------
!    if lw_tendency was written, read it.
!--------------------------------------------------------------------
          if (was_doing_edt .or. was_doing_entrain) then
            call read_data (unit, lw_tendency)

!---------------------------------------------------------------------
!    if edt_mod or entrain_mod was not active in the job which wrote the
!    restart file but it is active in the current job, initialize
!    lw_tendency to 0.0 and put a message in the output file.  
!---------------------------------------------------------------------
          else
            if (doing_edt .or. doing_entrain) then
              lw_tendency = 0.0
              call error_mesg ('physics_driver_mod', &
             ' initializing lw_tendency to 0.0, since it not present'//&
                  ' in physics_driver.res file', NOTE)
            endif
          endif

!---------------------------------------------------------------------
!    close the io unit associated with physics_driver.res. set flag
!    to indicate that the restart data has been processed. 
!---------------------------------------------------------------------
          call close_file (unit)
          success = .true.
        endif  ! (vers >=5)

!---------------------------------------------------------------------
!    if there is no physics_driver.res, set the remaining module
!    variables to 0.0
!---------------------------------------------------------------------
      else
        diff_t = 0.0
        diff_m = 0.0
        diff_cu_mo = 0.0
        pbltop     = -999.0
        convect = .false.
      endif  ! present(.res)

!--------------------------------------------------------------------
!    if a version of physics_driver.res containing the needed data is
!    not present, check for the presence of the radturbten data in 
!    strat_cloud.res.
!--------------------------------------------------------------------
      if ( .not. success) then
        if (doing_strat()) then
          if (file_exist('INPUT/strat_cloud.res')) then
            unit = open_restart_file ('INPUT/strat_cloud.res', 'read')
            read (unit, iostat=io, err=142) vers, vers2

!----------------------------------------------------------------------
!    if an i/o error does not occur, then the strat_cloud.res file 
!    contains the variable radturbten. rewind and read. close file upon
!    completion.
!----------------------------------------------------------------------
142         continue
            if (io == 0) then
              call error_mesg ('physics_driver_mod',  &
                'reading pre-version number strat_cloud.res file, '//&
                 'reading  radturbten', NOTE)
              rewind (unit)
              call read_data (unit, radturbten)
              call close_file (unit)

!---------------------------------------------------------------------
!    if the eor was reached (io /= 0), then the strat_cloud.res file
!    does not contain the radturbten data.  set values to 0.0 and
!    put a note in the output file.
!---------------------------------------------------------------------
            else
              radturbten = 0.0
              call error_mesg ('physics_driver_mod',  &
                  'neither strat_cloud.res nor physics_driver.res '//&
                   'contain the radturbten data, setting it to 0.0', &
                                                                NOTE)
            endif

!----------------------------------------------------------------------
!    if strat_cloud.res is not present, set radturbten to 0.0.
!----------------------------------------------------------------------
          else
            radturbten = 0.0
            call error_mesg ('physics_driver_mod',  &
              'setting radturbten to zero, no strat_cloud.res '//&
               'file present, data not in physics_driver.res', NOTE)
          endif
        endif

!--------------------------------------------------------------------
!    check if the lw_tendency data is in edt.res.
!--------------------------------------------------------------------
        if (doing_edt) then
          if (file_exist('INPUT/edt.res')) Then
            unit = open_restart_file ('INPUT/edt.res', 'read')
            read (unit, iostat=io, err=143) vers, vers2

!----------------------------------------------------------------------
!    if an i/o error does not occur, then the edt.res file 
!    contains the variable lw_tendency. rewind and read. close file 
!    upon completion.
!----------------------------------------------------------------------
143         continue
            if (io == 0) then
              call error_mesg ('physics_driver_mod',  &
                'reading pre-version number edt.res file, &
                 &reading  lw_tendency', NOTE)
              rewind (unit)
              call read_data (unit, lw_tendency)
              call close_file (unit)

!---------------------------------------------------------------------
!    if the eor was reached (io /= 0), then the edt.res file 
!    does not contain the lw_tendency data.  set values to 0.0 and
!    put a note in the output file.
!---------------------------------------------------------------------
            else
              lw_tendency = 0.0
              call error_mesg ('physics_driver_mod',  &
                  'neither edt.res nor physics_driver.res &
                   &contain the lw_tendency data, setting it to 0.0', &
                                                                NOTE)
            endif

!----------------------------------------------------------------------
!    if edt.res is not present, set lw_tendency to 0.0.
!----------------------------------------------------------------------
          else
            lw_tendency = 0.0
            call error_mesg ('physics_driver_mod',  &
               'setting lw_tendency to zero, no edt.res &
               &file present, data not in physics_driver.res', NOTE)
          endif
        endif

!--------------------------------------------------------------------
!    check if the lw_tendency data is in entrain.res. only 1 form of
!    entrain.res has ever existed, containing only the lw_tendency
!    variable, so it can be read without further checking.
!--------------------------------------------------------------------
        if (doing_entrain) then
          if (file_exist('INPUT/entrain.res')) Then
            unit = open_restart_file ('INPUT/entrain.res', 'read')
            call read_data (unit, lw_tendency)
            call close_file (unit)

!----------------------------------------------------------------------
!    if entrain.res is not present, set lw_tendency to 0.0.
!----------------------------------------------------------------------
          else
            lw_tendency = 0.0
            call error_mesg ('physics_driver_mod',  &
              'setting lw_tendency to zero, no entrain.res &
               &file present, data not in physics_driver.res', NOTE)
          endif
        endif
      endif  ! (.not. success)

!----------------------------------------------------------------------


 end subroutine read_restart_file     

!#####################################################################
! <SUBROUTINE NAME="read_restart_nc">
!  <OVERVIEW>
!    read_restart_nc will read the physics_driver.res file and process
!    its contents. if no restart data can be found, the module variables
!    are initialized to flag values.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_restart_nc will read the physics_driver.res file and process
!    its contents. if no restart data can be found, the module variables
!    are initialized to flag values.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_restart_nc
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine read_restart_nc

!---------------------------------------------------------------------
!    read_restart_file will read the physics_driver.res file and process
!    its contents. if no restart data can be found, the module variables
!    are initialized to flag values.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real  :: vers, vers2
      character(len=8) :: chvers
      real  :: was_doing_strat=0., was_doing_edt=0., was_doing_entrain=0.
      logical  :: success = .false.
      character(len=64) :: fname = 'INPUT/physics_driver.res.nc'
      real, dimension(size(convect,1), size(convect,2)) :: r_convect
!--------------------------------------------------------------------
!   local variables:
!
!      vers              restart version number if that is contained in 
!                        file; otherwise the first word of first data 
!                        record of file
!      vers2             second word of first data record of file
!      was_doing_strat   logical indicating if strat_cloud_mod was 
!                        active in job which wrote restart file
!      was_doing_edt     logical indicating if edt_mod was active
!                        in job which wrote restart file
!      was_doing_entrain logical indicating if entrain_mod was active
!                        in job which wrote restart file
!      success           logical indicating that restart data has been
!                        processed
!
!---------------------------------------------------------------------      
                    
      if(file_exist(fname)) then
         if(mpp_pe() == mpp_root_pe()) call mpp_error ('physics_driver_mod', &
            'Reading NetCDF formatted restart file: INPUT/physics_driver.res.nc', NOTE)
         call read_data(fname, 'vers', vers, no_domain=.true.)
         call read_data(fname, 'doing_strat', was_doing_strat, no_domain=.true.)
         call read_data(fname, 'doing_edt', was_doing_edt, no_domain=.true.)
         call read_data(fname, 'doing_entrain', was_doing_entrain, no_domain=.true.)
!---------------------------------------------------------------------
!    read the contribution to diffusion coefficient from cumulus
!    momentum transport.
!---------------------------------------------------------------------
         call read_data (fname, 'diff_cu_mo', diff_cu_mo)
         
!---------------------------------------------------------------------
!    pbl top is present in file versions 2 and up. if not present,
!    set a flag.
!---------------------------------------------------------------------
         call read_data (fname, 'pbltop', pbltop)

!---------------------------------------------------------------------
!    the temperature and momentum diffusion coefficients are present
!    beginning with v3. if not prsent, set to 0.0.
!---------------------------------------------------------------------
         call read_data (fname, 'diff_t', diff_t)
         call read_data (fname, 'diff_m', diff_m)

!---------------------------------------------------------------------
!    a flag indicating columns in which convection is occurring is
!    present beginning with v4. if not present, set it to .false.
!---------------------------------------------------------------------
         convect = .false.
         r_convect = 0.
         call read_data (fname, 'convect', r_convect)
         where(r_convect .GT. 0.) 
            convect = .true.
         end where
         
!---------------------------------------------------------------------
!    radturbten may be present in versions 5 onward, if strat_cloud_mod
!    was active in the job writing the .res file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    if radturbten was written, read it.
!--------------------------------------------------------------------
          if (was_doing_strat .GT. 0.) then
            call read_data (fname, 'radturbten', radturbten)

!---------------------------------------------------------------------
!    if strat_cloud_mod was not active in the job which wrote the 
!    restart file but it is active in the current job, initialize
!    radturbten to 0.0 and put a message in the output file.  
!---------------------------------------------------------------------
          else
            if (doing_strat()) then
              radturbten = 0.0
              call error_mesg ('physics_driver_mod', &
              ' initializing radturbten to 0.0, since it not present'//&
                            ' in physics_driver.res.nc file', NOTE)
            endif
          endif

!--------------------------------------------------------------------
!    if lw_tendency was written, read it.
!--------------------------------------------------------------------
          if (was_doing_edt .GT. 0. .or. was_doing_entrain .GT. 0.) then
            call read_data (fname, 'lw_tendency', lw_tendency)

!---------------------------------------------------------------------
!    if edt_mod or entrain_mod was not active in the job which wrote the
!    restart file but it is active in the current job, initialize
!    lw_tendency to 0.0 and put a message in the output file.  
!---------------------------------------------------------------------
          else
            if (doing_edt .or. doing_entrain ) then
              lw_tendency = 0.0
              call error_mesg ('physics_driver_mod', &
             ' initializing lw_tendency to 0.0, since it not present'//&
                  ' in physics_driver.res.nc file', NOTE)
            endif
          endif
       endif
!----------------------------------------------------------------------
     
     
end subroutine read_restart_nc


!#####################################################################
! <SUBROUTINE NAME="check_args">
!  <OVERVIEW>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
!                        u, v, t, q, r, um, vm, tm, qm, rm,             &
!                        udt, vdt, tdt, qdt, rdt, mask, kbot)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="udt" TYPE="real">
!   zonal wind tendency
!  </IN>
!  <IN NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </IN>
!  <IN NAME="tdt" TYPE="real">
!   temperature tendency
!  </IN>
!  <IN NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </IN>
!  <IN NAME="rdt" TYPE="real">
!   multiple tracer tendencies
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
subroutine check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
                        u, v, t, q, r, um, vm, tm, qm, rm,             &
                        udt, vdt, tdt, qdt, rdt, mask, kbot)

!----------------------------------------------------------------------
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!-----------------------------------------------------------------------

real,    dimension(:,:),    intent(in)          :: lat, lon, area
real,    dimension(:,:,:),  intent(in)          :: p_half, p_full,   &
                                                   z_half, z_full,   &
                                                   u, v, t, q, um, vm, &
                                                   tm, qm
real,    dimension(:,:,:,:),intent(in)          :: r, rm
real,    dimension(:,:,:),  intent(in)          :: udt, vdt, tdt, qdt
real,    dimension(:,:,:,:),intent(in)          :: rdt
real,    dimension(:,:,:),  intent(in),optional :: mask
integer, dimension(:,:),    intent(in),optional :: kbot

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!
!   intent(in), optional:
!
!       mask        mask that designates which levels do not have data
!                   present (i.e., below ground); 0.=no data, 1.=data
!       kbot        lowest level which has data
!                   note:  both mask and kbot must be present together.
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer ::  id, jd, kd  ! model dimensions on the processor  
      integer ::  ierr        ! error flag

!--------------------------------------------------------------------
!    define the sizes that the arrays should be.
!--------------------------------------------------------------------
      id = size(u,1) 
      jd = size(u,2) 
      kd = size(u,3) 

!--------------------------------------------------------------------
!    check the dimensions of each input array. if they are incompat-
!    ible in size with the standard, the error flag is set to so
!    indicate.
!--------------------------------------------------------------------
      ierr = 0
      ierr = ierr + check_dim (lat, 'lat',  id,jd)
      ierr = ierr + check_dim (lon, 'lon',  id,jd)
      ierr = ierr + check_dim (area,'area', id,jd)

      ierr = ierr + check_dim (p_half,'p_half', id,jd,kd+1)
      ierr = ierr + check_dim (p_full,'p_full', id,jd,kd)
      ierr = ierr + check_dim (z_half,'z_half', id,jd,kd+1)
      ierr = ierr + check_dim (z_full,'z_full', id,jd,kd)

      ierr = ierr + check_dim (u, 'u',  id,jd,kd)
      ierr = ierr + check_dim (v, 'v',  id,jd,kd)
      ierr = ierr + check_dim (t, 't',  id,jd,kd)
      ierr = ierr + check_dim (q, 'q',  id,jd,kd)
      ierr = ierr + check_dim (um,'um', id,jd,kd)
      ierr = ierr + check_dim (vm,'vm', id,jd,kd)
      ierr = ierr + check_dim (tm,'tm', id,jd,kd)
      ierr = ierr + check_dim (qm,'qm', id,jd,kd)

      ierr = ierr + check_dim (udt,'udt', id,jd,kd)
      ierr = ierr + check_dim (vdt,'vdt', id,jd,kd)
      ierr = ierr + check_dim (tdt,'tdt', id,jd,kd)
      ierr = ierr + check_dim (qdt,'qdt', id,jd,kd)

      if (nt > 0) then
        ierr = ierr + check_dim (r,  'r',   id,jd,kd,nt)
        ierr = ierr + check_dim (rm, 'rm',  id,jd,kd,nt)
      endif
      if (ntp > 0) then
        ierr = ierr + check_dim (rdt,'rdt', id,jd,kd,ntp)
      endif

!--------------------------------------------------------------------
!    if any problems were detected, exit with an error message.
!--------------------------------------------------------------------
      if (ierr > 0) then
        call error_mesg ('physics_driver_mod', 'bad dimensions', FATAL)
      endif

!--------------------------------------------------------------------
!    set a flag to indicate that this check was done and need not be
!    done again.
!--------------------------------------------------------------------
      do_check_args = .false.

!-----------------------------------------------------------------------


      end subroutine check_args


!#######################################################################
! <FUNCTION NAME="check_dim_2d">
!  <OVERVIEW>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_2d (data,name,id,jd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd" TYPE="integer">
!   expected i and j dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_2d (data,name,id,jd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:) :: data
character(len=*), intent(in)        :: name
integer, intent(in)                 :: id, jd
integer                             :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd      expected i and j dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
             'dimension 1 of argument ' //  &
              name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

!----------------------------------------------------------------------

      end function check_dim_2d

!#######################################################################
! <FUNCTION NAME="check_dim_3d">
!  <OVERVIEW>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_3d (data,name,id,jd, kd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd" TYPE="integer">
!   expected i, j and k dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_3d (data,name,id,jd,kd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_3d compares the size of thr1eedimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:,:) :: data
character(len=*), intent(in)          :: name
integer, intent(in)                   :: id, jd, kd
integer  ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd,kd   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
              'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_3d


!#######################################################################
! <FUNCTION NAME="check_dim_4d">
!  <OVERVIEW>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_4d (data,name,id,jd, kd, nt) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd, nt" TYPE="integer">
!   expected i, j, k and 4th dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_4d (data,name,id,jd,kd,nt) result (ierr)

!--------------------------------------------------------------------
!    check_dim_4d compares the size of four dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------
real,    intent(in), dimension(:,:,:,:) :: data
character(len=*), intent(in)            :: name
integer, intent(in)                     :: id, jd, kd, nt
integer                                 :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data          array to be checked
!     name          name associated with array to be checked
!     id,jd,kd,nt   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr          set to 0 if ok, otherwise is a count of the number
!                   of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,4) /= nt) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 4 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_4d



!#######################################################################
 
 
 
                end module physics_driver_mod
