                       module donner_deep_mod

use time_manager_mod,       only: time_manager_init, time_type, &
                                  set_date, get_time,   &
                                  get_calendar_type, &
                                  operator(>=), operator (<)
use diag_manager_mod,       only: register_diag_field, send_data, &
                                  diag_manager_init
use field_manager_mod,      only: MODEL_ATMOS, field_manager_init

!  tracer_manager_init not yet available:
!use tracer_manager_mod,     only: tracer_manager_init, get_tracer_names
use tracer_manager_mod,     only:                      get_tracer_names

use sat_vapor_pres_mod,     only: lookup_es, sat_vapor_pres_init
use fms_mod,                only: fms_init, mpp_pe, mpp_root_pe,  &
                                  file_exist,  check_nml_error,  &
                                  error_mesg, FATAL, WARNING, NOTE,  &
                                  close_file, open_namelist_file,    &
                                  stdlog, write_version_number,  &
                                  read_data, write_data,    &
                                  open_restart_file
use constants_mod,          only: constants_init, DENS_H2O, RDGAS,   &
                                  GRAV, CP_AIR, pie=>PI
use column_diagnostics_mod, only: column_diagnostics_init, &
                                  initialize_diagnostic_columns, &
                                  column_diagnostics_header, &
                                  close_column_diagnostics_units


implicit none
private

!--------------------------------------------------------------------
!         module to compute the effects of deep convection
!
!         Primary remaining work:
!               optimization / cleanup of code -- current code adds 
!               ~ 40% to total N30L40 FMS model time, after donner_deep
!               is spunup.
!     
!
!--------------------------------------------------------------------





!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: donner_deep.f90,v 12.0 2005/04/14 15:40:45 fms Exp $'
character(len=128)  :: tagname =  '$Name: lima $'


!--------------------------------------------------------------------
!---interfaces------

public   &
        donner_deep_init, donner_deep, donner_deep_avg, donner_deep_end

private   &
!  routines called by donner_deep_init:
         register_fields, initialize_variables, allocate_variables, &
         read_restart_donner_deep, &
!  routines called by donner_deep:
         initialize_local_variables, donner_convection_driver, &
         define_output_fields, cape_calculation_driver,   &
         donner_deep_netcdf,  deallocate_local_variables, &
!  lower-level routines: 
         cupar, generate_cape_sounding, polat_vect, calculate_cape,  &
         ieq_x, satad, lcl, mulsub_vect, cloudm_vect, simult, meens, &
         mesub, micro,  prean, andge, cell_liquid_size_comp, &
         strat_cloud_donner_tend , donner_deep_sum, &
!  column diagnostics routines
        donner_column_init, donner_column_control,    &
        donner_column_input_fields, donner_column_cape_call, &
        donner_column_end_of_step

!----------------------------------------------------------------------
!   define a private derived type donner_conv_type containing
!   the following components:
!
!   cecon          normalized cell condensation/deposition
!                  [ K/sec ]
!   ceefc          normalized cell entropy-flux convergence [ K/sec ]
!                  (excludes convergence of surface flux) Entropy-flux
!                  convergence divided by (p0/p)**(rd/cp).
!   cell_ice_geneff_diam  
!                  cell ice generalized effective diameter.
!                  default is smallest data value given in andge
!                  subroutine (micrometers)
!   cell_liquid_eff_diam  
!                  cell liquid effective diameter. defined 
!                  using Bower parameterization or from namelist
!                  input. (micrometers)
!   cememf_mod     normalized moisture forcing, cells+meso, after 
!                  any modifications needed to prevent the creation
!                  of negative values [ kg(H2O)/kg/sec ]
!   cemfc          normalized cell moisture-flux convergence
!                  (excludes convergence of surface moisture flux)
!                  [ kg(H2O)/kg/sec ]
!   cmus           normalized mesoscale-updraft deposition
!                  [ kg(H2O)/kg/sec]
!   cual           cloud fraction, cells+meso, normalized by a(1,p_b)
!   cuqi           cell ice content (kg(ice)/kg)
!   cuql           cell liquid content (kg{water)/kg)
!   dgeice         mesoscale ice generalized effective size, defined 
!                  as in Fu  (1996, J. Clim.) (micrometers)
!   dmeml          mass flux in mesoscale downdraft [ kg/((m**2) s) ]
!                  (normalized by a(1,p_b)) 
!   ecds           normalized convective downdraft evaporation
!                  [ kg(H2O)/kg/sec ]
!   eces           normalzed convective-updraft evporation/sublimation
!                  [ kg(H2O)/kg/sec ]
!   elt            normalized melting [ K/sec ]
!   emds           normalized mesoscale-downdraft sublimation
!                  [ kg(H2O)/kg/sec ]
!   emes           normalized mesoscale-updraft sublimation
!                  [ kg(H2O)/kg/sec ]
!   fre            normalized freezing [ K/sec ]
!   pb_v           pressure at base of cumulus updrafts (Pa)
!   pmd_v          pressure at top of mesoscale downdraft (Pa)
!   pztm_v         pressure at top of mesoscale updraft (Pa)
!   qmes           normalized mesoscale moisture-flux convergence
!                  [ kg(H2O)/kg/sec ]
!   tmes           normalized mesoscale entropy-flux convergence
!                  [ K/sec ]
!                  Entropy-flux convergence is mesoscale component
!                  of second term in expression for cumulus thermal
!                  forcing in Fig. 3 of Donner (1993, JAS).
!   uceml          normalized mass fluxes in cell updrafts
!                  [ kg/((m**2)*s ] 
!   umeml          mass flux in mesoscale updraft [ kg/((m**2) s) ]
!                  (normalized by a(1,p_b)) 
!   wmms           normalized mesoscale deposition of water vapor from
!                  cells [ kg(H2O)/kg/sec ]
!   wmps           normalized mesoscale redistribution of water vapor
!                  from cells [ kg(H2O)/kg/sec ]
!   xice           mesoscale ice mass mixing ratio (kg(ice)/kg)
!   qtceme         tracer tendencies due to donner_deep_mod 
!                  [ kg/kg/s ]
!   xgcm1          tracer fields that are subjected to donner_deep 
!                  transport
!                  [ kg/kg ]
!   qtren1         tracer time tendency due to cell scale motions
!                  [ kg/kg/s ]
!   qtmes1         tracer time tendency due to mesoscale  motions
!                  [ kg/kg/s ]
!   wtp1           redistribution of tracer from cellscale to mesoscale
!                  [ kg/kg/s ]
!   a1             fractional area of index-1 cu subensemble
!   amax           maximum value for a_1(p_b)
!                  See "a Bounds 6/7/97" notes
!   amos           upper limit on cloud fractional area based on
!                  moisture constraint See "Moisture Constraint," 
!                  8/8/97.
!   ampta1         area weighted mesoscale cloud fraction, normal-
!                  ized by a(1,p_b)
!   contot         ratio of convective to total precipitation
!   dcape          time rate of change of cape [ J/(kg s) ]
!   emdi_v         vertical integral of mesoscale-downdraft 
!                  sublimation
!   rcoa1          area weighted convective precipitation rate
!                  [ mm/day ]
!
!---------------------------------------------------------------------

private donner_conv_type

type donner_conv_type
     real, dimension(:,:,:), pointer  ::   &
          cecon=>NULL(), ceefc=>NULL(), cell_ice_geneff_diam=>NULL(),&
          cell_liquid_eff_diam=>NULL(), &
          cememf_mod=>NULL(), cemfc=>NULL(), cmus=>NULL(),  &
          cual=>NULL(), cuqi=>NULL(), cuql=>NULL(), dgeice=>NULL(), &
          dmeml=>NULL(), &
          ecds=>NULL(), eces=>NULL(), elt=>NULL(), emds=>NULL(),   &
          emes=>NULL(), fre=>NULL(), qmes=>NULL(), tmes=>NULL(),   &
          uceml=>NULL(), umeml=>NULL(), &
          wmps=>NULL(), wmms=>NULL(), xice=>NULL()
     real, dimension(:,:,:,:), pointer ::   &
           qtceme=>NULL(), xgcm1=>NULL(), qtren1=>NULL(), &
           qtmes1=>NULL(), wtp1=>NULL()
     real, dimension(:,:),   pointer  ::    &
          a1=>NULL(), amax=>NULL(), amos=>NULL(), ampta1=>NULL(), &
         contot=>NULL(), dcape=>NULL(), emdi_v=>NULL(), rcoa1=>NULL(),&
          pb_v=>NULL(), pmd_v=>NULL(), pztm_v=>NULL()
end type donner_conv_type

type(donner_conv_type), save    :: Don_conv

!----------------------------------------------------------------------
!   define a private derived type donner_cape_type containing
!   the following components:
!
!   coin           convective inhibition 
!                  energy required to lift parcel from level istart
!                  to level of free convection. [ J/kg ]
!   plcl           pressure at lifting condensation level [ Pa ]
!   plfc           pressure at level of free convection [ Pa ]
!                  height of plfc .le. height of plcl. if parcel 
!                  becomes buoyant below plcl, cin can be .lt. 0
!   plzb           pressure at level of zero buoyancy [ Pa ]
!   qint           vertically integrated column moisture
!                  [ kg (h20)/(m**2) ]
!   xcape          convective available potential energy. energy 
!                  released as parcel moves from level of free 
!                  convection to level of zero buoyancy, calculated on 
!                  convection step [ J/kg ].

!   the following variables are on the enhanced cape vertical grid
!   (ncap levels), with k index 1 closest to the ground:

!   parcel_r       parcel mixing ratio                [ kg/kg ]
!   parcel_t       parcel temperature                 [ K ]
!   cape_p         pressure levels of cape grid       [ hPa ]
!   env_r          environmental mixing ratio profile [ kg/kg ]
!   env_t          environmental temperature profile  [ K ]

!   the following variables are on model levels:

!   model_p        pressure profile used to define pressure 
!                  profile used in cape calculation [ hPa ]
!   model_r        mixing ratio profile used to define
!                  moisture profile used in cape calculation [ kg/kg ]
!   model_t        temperature profile used to define
!                  temperature profile used in cape calculation [ K ]
!       
!----------------------------------------------------------------------

private donner_cape_type

type donner_cape_type
     real, dimension(:,:), pointer ::    &
                         coin=>NULL(), plcl=>NULL(), plfc=>NULL(), &
                         plzb=>NULL(), qint=>NULL(), xcape=>NULL()
     real, dimension(:,:,:), pointer ::    &
                         parcel_r=>NULL(), parcel_t=>NULL(),   &
                         cape_p=>NULL(), env_r=>NULL(), env_t=>NULL()
     real, dimension (:,:,:), pointer ::  &
                         model_p=>NULL(), model_r=>NULL(),  &
                         model_t=>NULL()
end type donner_cape_type

type(donner_cape_type), save    :: Don_cape

!---------------------------------------------------------------------
!---namelist----

integer :: model_levels_in_sfcbl=2 ! number of levels at which the
                                   ! temperature and vapor profiles are 
                                   ! not allowed to change from lag 
                                   ! value when calculating the time 
                                   ! tendency of cape
integer :: donner_deep_freq = 1800 ! frequency of calling donner_deep 
                                   ! [ sec ]; must be <= 86400 
character(len=16) :: cell_liquid_size_type = 'bower' 
                                   ! choose either 'input' or 'bower'
character(len=16) :: cell_ice_size_type = 'default' 
                                   ! choose either 'input' or 'default'
real :: cell_liquid_eff_diam_input = -1.0 
                                   ! input cell droplet 
                                   ! effective diameter [ micrometers];
                                   ! needed when cell_liquid_size_type 
                                   ! == 'input'
real :: cell_ice_geneff_diam_input = -1.0 
                                   ! input cell ice generalized
                                   ! effective diameter [ micrometers ];
                                   ! needed when cell_ice_size_type 
                                   ! == 'input'
logical :: do_average = .false.    ! time-average donner cloud proper-
                                   ! ties for use by radiation package?

!   the following variables are used when the column diagnostics 
!   option is activated:

integer, parameter  :: MAX_PTS = 20 
                                   ! max nunber of diagnostic columns
real :: diagnostics_pressure_cutoff =  5000.   
                                   ! column data will be 
                                   ! printed on model levels with 
                                   ! pressures greater than this value 
                                   ! [ Pa ]
integer, dimension(6) :: diagnostics_start_time= (/ 0,0,0,0,0,0 /)
                                   ! integer specification of time for 
                                   ! column diagnostics to be activated
                                   ! [year, month, day, hour, min, sec ]

integer :: num_diag_pts_ij = 0     ! number of diagnostic columns 
                                   ! specified by global(i,j) 
                                   ! coordinates
integer :: num_diag_pts_latlon = 0 ! number of diagnostic columns 
                                   ! specified by lat-lon coordinates
integer, dimension(MAX_PTS) :: i_coords_gl = -100
                                   ! global i coordinates for ij 
                                   ! diagnostic columns
integer, dimension(MAX_PTS) :: j_coords_gl = -100
                                   ! global j coordinates for ij 
                                   ! diagnostic columns
real, dimension(MAX_PTS) :: lat_coords_gl = -999.
                                   ! latitudes for lat-lon  diagnostic
                                   ! columns [ degrees, -90. -> 90. ]
real, dimension(MAX_PTS) :: lon_coords_gl = -999.
                                   ! longitudes for lat-lon  diagnostic
                                   ! columns [ degrees, 0. -> 360. ]

namelist / donner_deep_nml /      &
                            model_levels_in_sfcbl, &
                            donner_deep_freq,  &
                            cell_liquid_size_type, cell_ice_size_type, &
                            cell_liquid_eff_diam_input, &
                            cell_ice_geneff_diam_input, &
                            do_average,         &
                            diagnostics_pressure_cutoff, &
                            diagnostics_start_time, &
                            num_diag_pts_ij, num_diag_pts_latlon, &
                            i_coords_gl, j_coords_gl, &
                            lat_coords_gl, lon_coords_gl


!--------------------------------------------------------------------
!--- public data ----------




!--------------------------------------------------------------------
!----private data-----------


!--------------------------------------------------------------------
!   list of restart versions readable by this module:
!
!   version 1 is original form  (pre 6/13/2001) 
!   version 2 adds the tprea1 array to the restart file (6/13/2001)
!   version 3 contains the time remaining until the next donner step,
!             rather than the actual time of the next donner step 
!             (2/4/2002)
!   version 4 contains the variables tempbl and ratpbl, needed for
!             the catendb closure (6/29/2002)
!   version 5 has tempbl and ratpbl dimensioned by model_levels_in_sfcbl
!             rather than nlev; also adds mtot, delta_ql, delta_qi and 
!             delta_qa as variables to be saved across timesteps (needed
!             when donner step does not equal the physics step). it also
!             has only conv_alarm, donner_frequency as control info, 
!             and does not print out calendar type, compatible with 
!             newer versions of other restart files.  (12/31/2002)
!   version 6 has added the tracer tendencies due to donner_deep mod,
!             needed in case donner_deep_mod is not called on every
!             timestep. (9/9/03)

integer, dimension(6)  :: restart_versions = (/ 1, 2, 3, 4, 5, 6 /)

!--------------------------------------------------------------------
!   these arrays contain information used to provide the variables
!   needed by moist_processes_mod on every timestep, so they must be
!   preserved across timesteps, in case donner_deep is not called on 
!   every step:
!      
!     cememf         normalized moisture forcing, cells+meso 
!                    [ kg(H2O)/kg/sec ]
!     cemetf         normalized thermal forcing, cells+meso [ K/sec ]
!                    (excludes convergence of surface heat flux)
!                    Cumulus thermal forcing defined as in Fig. 3 of 
!                    Donner (1993, JAS).
!     tprea1         area weighted total normalized precipitation 
!                    [ mm/day ]
!     mass_flux      total mass flux at model levels, convective +
!                    mesoscale  [ kg/(m**2)/sec) ]
!     delta_ql       increment in cloud liquid over the physics timestep
!                    resulting from donner convection [ kg/kg/sec ]
!     delta_qi       increment in cloud ice over the physics timestep
!                    resulting from donner convection [ kg/kg/sec ]
!     delta_qa       increment in cloud area over the physics timestep
!                    resulting  from donner convection [ 1/sec ]
!     tracer_tends   tendencies to tracer fields due to donner_deep
!                    mod [ kg/kg/sec ]
!
!   these fields are needed by the radiation package and so must be 
!   saved between time steps:
!
!     cell_cloud_frac     fractional area of convective cells in grid 
!                         box [ dimensionless ] 
!     cell_liquid_amt     liquid water content of convective cells 
!                         [ kg(h2o)/kg(air) ] 
!     cell_liquid_size    assumed effective size of cell liquid drops
!                         [ micrometers ]
!     cell_ice_amt        ice water content of cells 
!                         [ kg(h2o)/kg(air) ]
!     cell_ice_size       generalized effective diameter for ice in
!                         convective cells [ micrometers ]
!     meso_cloud_frac     fractional area of mesoscale clouds in grid 
!                         box [ dimensionless ]
!     meso_liquid_amt     liquid water content in mesoscale clouds
!                         currently set to 0.0
!                         [ kg(h2o)/kg(air) ]
!     meso_liquid_size    assumed effective size of mesoscale drops
!                         currently set to 0.0 [ micrometers ]
!     meso_ice_amt        ice water content of mesoscale elements 
!                         [ kg(h2o)/kg(air) ]
!     meso_ice_size       generalized ice effective size for anvil ice
!                         [ micrometers ]
!     nsum                number of time levels of data contained in
!                         the accumulation arrays; needed when time-
!                         averaging of cloud properties to be used in
!                         radiation package is desired
!
!--------------------------------------------------------------------

real,    dimension(:,:,:), allocatable  :: cemetf        
real,    dimension(:,:,:), allocatable  :: cememf    
real,    dimension(:,:),   allocatable  :: tprea1     
real,    dimension(:,:,:), allocatable  :: mass_flux               
real,    dimension(:,:,:), allocatable  :: delta_ql                
real,    dimension(:,:,:), allocatable  :: delta_qi                
real,    dimension(:,:,:), allocatable  :: delta_qa               
real,    dimension(:,:,:,:), allocatable  :: tracer_tends           
real,    dimension(:,:,:), allocatable  :: cell_cloud_frac
real,    dimension(:,:,:), allocatable  :: cell_liquid_amt
real,    dimension(:,:,:), allocatable  :: cell_liquid_size
real,    dimension(:,:,:), allocatable  :: cell_ice_amt
real,    dimension(:,:,:), allocatable  :: cell_ice_size
real,    dimension(:,:,:), allocatable  :: meso_cloud_frac
real,    dimension(:,:,:), allocatable  :: meso_liquid_amt
real,    dimension(:,:,:), allocatable  :: meso_liquid_size
real,    dimension(:,:,:), allocatable  :: meso_ice_amt
real,    dimension(:,:,:), allocatable  :: meso_ice_size
integer, dimension(:,:)  , allocatable  :: nsum
 
!--------------------------------------------------------------------
!   this field contains a sum over time, and so must be preserved
!   across timesteps:
!
!     omint_acc        time-integrated low level displacement [ Pa ]
!
!--------------------------------------------------------------------
real, dimension(:,:), allocatable  :: omint_acc   

!--------------------------------------------------------------------
!   these fields are used in defining time derivatives, and so must 
!   be preserved across timesteps:
!
!     qint_lag         vertically integrated column moisture,
!                      calculated on step prior to convection step
!                      [ kg (h20)/(m**2) ]
!     xcape_lag        convective available potential energy. energy 
!                      released as parcel moves from level of free 
!                      convection to level of zero buoyancy [ J/kg ]
!                      calculated on step prior to convection step
!
!--------------------------------------------------------------------
real, dimension(:,:), allocatable  :: xcape_lag  
real, dimension(:,:), allocatable  :: qint_lag 

!--------------------------------------------------------------------
!   these fields are used to retain the low-level temperature and mix-
!   ing ratio profiles from the lag step, so that in the the cape cal-
!   culation on the mid step, the lowest model_levels_in_sfcbl levels 
!   of these fields are the same, eliminating a source of temporal 
!   noise in the cape time-tendency field (catendb closure), and so 
!   must be preserved across timesteps:
!
!     tempbl    temperature profile to be used in cape calculation at
!               the lowest model_levels_in_sfcbl levels [ deg K ]
!     ratpbl    mixing ratio profile to be used in cape calculation at
!               the lowest model_levels_in_sfcbl levels [ kg/kg ]
!
!--------------------------------------------------------------------
real, dimension(:,:,:), allocatable  :: tempbl     
real, dimension(:,:,:), allocatable  :: ratpbl 
      
!------------------------------------------------------------------
!   module loop and dimension variables

integer :: idf, jdf      ! processor subdomain dimensions
integer :: nlev          ! number of model layers
integer :: isize, jsize  ! physics window sizes

!------------------------------------------------------------------
!   module control variables

logical :: coldstart=.false.            ! is this a donner_deep 
                                        ! coldstart ?
logical :: module_is_initialized =  &
                              .false.   ! has this module been 
                                        ! initialized ?
logical :: first_call = .true.          ! is this the first call to
                                        ! donner_deep_mod in this job ?
logical :: conv_calc_on_this_step =  &
                               .false.  ! is this a step on which to 
                                        ! calculate convection ?
logical  :: conv_calc_on_next_step = &
                               .false.  ! is the next step a convection
                                        ! calculation step, i.e., is 
                                        ! this a lag-time cape calcul-
                                        ! ation step?
logical  :: do_input_cell_liquid_size =   &
                               .false.  ! cell liquid size to be input ?
logical  :: do_bower_cell_liquid_size =   &
                               .false.  ! cell liquid size from bower 
                                        ! calculation ?
logical  :: do_input_cell_ice_size =   &
                               .false.  ! cell ice size to be input ?
logical  :: do_default_cell_ice_size =  &
                               .false.  ! use default cell ice size ?
integer  :: total_pts                   ! total number of points in 
                                        ! the processor's subdomain
integer  :: pts_processed_conv          ! number of points processed 
                                        ! during current convection 
                                        ! calculation
integer  :: pts_processed_cape          ! number of points processed
                                        ! during current lag-time cape
                                        ! calculation
integer  :: conv_alarm                  ! time remaining until next
                                        ! convection calculation [ sec ]
integer  :: cape_alarm                  ! time remaining until next
                                        ! lag-time cape calculation 
                                        ! [ sec ]

!---------------------------------------------------------------------
!  constants required for deep cumulus parameterization. these will
!  ultimately be "used" from a constants module.
 
real, parameter   ::  rgas = 8.314320E+03    ! universal gas constant
                                             ! [ J/(kg K) ]
real, parameter   ::  wtmair = 2.896440E+01  ! molecular wt of air 
real, parameter   ::  wtmh2o = 1.801534E+01  ! molecular wt of water
real, parameter   ::  cpi    = 1.004840E+03  ! specific heat of dry air 
                                             ! at constant pressure 
                                             ! [ J/(kg K) ]
real, parameter   ::  gravm  = 9.806650      ! acceleration of gravity
                                             ! [ m/(sec**2) ]
real, parameter   ::  cpv    = 1850.         ! specific heat of water
                                             ! vapor at constant pres-
                                             ! sure [ J/(kg K) ]
real              ::  rair   = rgas/wtmair   ! gas constant for dry air
                                             ! [ J/(kg K) ]
real, parameter   ::  rvap   = rgas/wtmh2o   ! gas constant for water
                                             ! vapor [ J/(kg K) ]
real              ::  rh2o   = 461.          ! gas constant for water
                                             ! vapor [ J/(kg K) ]
real, parameter   ::  LATICE=3.336E05        ! latent heat of fusion
                                             ! [ J/kg ]
real, parameter   ::  latvap = 2.5104e06     ! latent heat of 
                                             ! vaporization [ J/kg ]
real, parameter   ::  rocp   = 0.622         ! ratio of molecular 
                                             ! weights of water vapor 
                                             ! and dry air
real              ::  epsilo = 0.622         ! ratio of molecular 
                                             ! weights of water vapor 
                                             ! and dry air

!--------------------------------------------------------------------
!   module internal parameters

integer, parameter :: kpar=7             !  number of cumulus    
                                         !  subensembles
integer, parameter :: ncap=100           !  number of levels in cloud 
                                         !  model
real, parameter    :: pdeep_cv = 500.e02 !  required pressure difference
                                         !  between GCM level closest to
                                         !  ground and pressure at level
                                         !  of zero buoyancy for deep 
                                         !  convection  to occur (Pa).
!real, parameter    :: cdeep_cv = 10.    !  maximum value of convective 
real, parameter    :: cdeep_cv = 100.    !  maximum value of convective
                                         !  inhibition (J/kg) that 
                                         !  allows convection. Value of 
                                         !  10 suggested by Table 2 in 
                                         !  Thompson et al. (1979, JAS).

real, parameter    :: pdeep_mc = 200.e02 !  pressure thickness (Pa) 
                                         !  required for mesoscale circ-
                                         !  ulation. It refers to the 
                                         !  least penetrative ensemble 
                                         !  member. For this check to 
                                         !  function properly, the en-
                                         !  trainment coefficient 
                                         !  in Cloudm for kou=1 must be
                                         !  the largest entrainment
                                         !  coefficient. 

logical            :: do_donner_tracer=    &
                                 .false. !  tracers are to be trans-
                                         !  ported by the donner con-
                                         !  vection scheme ?

integer            :: ncont = 0          !  number of tracers being 
                                         !  transported by donner
                                         !  convection
character(len=32), dimension(:), allocatable  ::  &
                      tracername        !  names of the tracers

!--------------------------------------------------------------------
!   module internal parameters related to drop and ice particle sizes
!   and concentrations
 
real, parameter :: r_conv_land  = 10.0      ! radius conv cld  
                                            ! over land (micrometers)
real, parameter :: r_conv_ocean = 16.0      ! radius conv cld  
                                            ! over ocean (micrometers)
real, parameter :: N_land       = 600*1.0e6 ! droplet number conc 
                                            ! over land (m**-3)
real, parameter :: N_ocean      = 150*1.0e6 ! droplet number conc 
                                            ! over ocean (m**-3)
real, parameter :: delz_land    = 500.0     ! cloud land vert depth (m) 
real, parameter :: delz_ocean   = 1500.0    ! cloud ocean vert depth (m)
real, parameter :: cell_liquid_eff_diam_def =   &
                                  15.0      ! default cell liquid eff 
                                            ! diameter (micrometers)
real, parameter :: cell_ice_geneff_diam_def =    &
                                  13.3      ! default cell ice general-
                                            ! ized effective diameter 
                                            ! (micrometers)

!--------------------------------------------------------------------
!   variables for netcdf diagnostics:

integer    :: id_cemetf_deep, id_ceefc_deep, id_cecon_deep, &
              id_cemfc_deep, id_cememf_deep, id_cememf_mod_deep, &
              id_cual_deep, id_fre_deep, id_elt_deep, &
              id_cmus_deep, id_ecds_deep, id_eces_deep, &
              id_emds_deep, id_emes_deep, id_qmes_deep,&
              id_wmps_deep, id_wmms_deep, id_tmes_deep,&
              id_dmeml_deep, id_uceml_deep, id_umeml_deep, &
              id_xice_deep,  id_dgeice_deep, id_dgeliq_deep,  &
              id_cuqi_deep, id_cuql_deep, &
              id_plcl_deep, id_plfc_deep, id_plzb_deep, &
              id_xcape_deep, id_coin_deep,  &
              id_dcape_deep, id_qint_deep, id_a1_deep, &
              id_amax_deep, id_amos_deep, &
              id_tprea1_deep, id_ampta1_deep, &
              id_omint_deep, id_rcoa1_deep

integer, dimension(:), allocatable :: id_qtren1, id_qtmes1, &
                                      id_wtp1, id_qtceme
integer, dimension(:), allocatable :: id_qtren1_col, id_qtmes1_col, &
                                      id_wtp1_col, id_qtceme_col


real              :: missing_value = -999.
character(len=16) :: mod_name = 'donner_deep'

!--------------------------------------------------------------------
!   variables for column diagnostics option
!
!   arrays dimensioned by number of diagnostic columns:
!    col_diag_unit     unit numbers for each column's output file 
!    col_diag_lon      each column's longitude 
!                      [ degrees, 0 < lon < 360 ]
!    col_diag_lat      each column's latitude
!                      [degrees, -90 < lat < 90 ]
!    col_diag_j        each column's j index (processor coord-
!                      inates)
!    col_diag_i        each column's i index (processor coord-
!                      inates) 
!
!    arrays dimensioned by MAX_PTS, but with values defined only for the
!    number of diagnostic columns in the current physics window 
!    (ncols_in_window):
!    i_dc              column's i index (window coordinates)
!    j_dc              column's j index (window coordinates)
!    unit_dc           column's output file unit number
!    igl_dc            column's i index (processor coordinates)
!    jgl_dc            column's j index (processor coordinates)
!
!    array dimensioned by the number of jrows on the processor:
!    do_column_diagnostics   is there a diagnostic column on this jrow ? 
integer, dimension(:), allocatable :: col_diag_unit
real   , dimension(:), allocatable :: col_diag_lon, col_diag_lat   
integer, dimension(:), allocatable :: col_diag_j, col_diag_i        

integer, dimension(MAX_PTS) :: i_dc=-99                     
integer, dimension(MAX_PTS) :: j_dc=-99                      
integer, dimension(MAX_PTS) :: unit_dc=-1     
integer, dimension(MAX_PTS) :: igl_dc=-99  
integer, dimension(MAX_PTS) :: jgl_dc=-99

logical, dimension(:), allocatable :: do_column_diagnostics

!   scalar variables for column diagnostics
!
type(time_type) :: Time_col_diagnostics  ! time in model simulation at 
                                         ! which to activate column 
                                         ! diagnostics 
logical         :: in_diagnostics_window = .false.
                                         ! are column diagnostics 
                                         ! desired anywhere in current 
                                         ! window ?  
integer         :: num_diag_pts = 0      ! total number of activated 
                                         ! diagnostics columns
integer         :: ncols_in_window=0     ! number of activated diag-
                                         ! nostic columns in this window
logical         :: column_diagnostics_desired = .false. 
                                         ! are column diagnostics 
                                         ! requested in any column ?
integer         :: kstart_diag=-99       ! output array elements for 
                                         ! model levels with k index 
                                         ! > kstart_diag will be written
                                         ! to the output file

!---------------------------------------------------------------------



                        contains




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################

subroutine donner_deep_init (lonb, latb, pref, axes, Time,  &
                             tracers_in_donner)

!---------------------------------------------------------------------
!    donner_deep_init is the constructor for donner_deep_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real,    dimension(:), intent(in)            :: lonb, latb, pref
integer, dimension(4), intent(in)            :: axes
type(time_type),       intent(in)            :: Time
logical, dimension(:), intent(in), optional  :: tracers_in_donner

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      lonb         array of model longitudes on cell boundaries  
!                   [ radians ]
!      latb         array of model latitudes on cell boundaries
!                   [ radians ]
!      pref         array of reference pressures at full levels (plus 
!                   surface value at nlev+1), based on 1013.25 hPa pstar
!                   [ Pa ]
!      axes         data axes for diagnostics
!      Time         current time [ time_type ]
!
!      tracers_in_donner 
!                   logical array indicating which of the activated 
!                   tracers are to be transported by donner_deep_mod
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      integer                 :: unit, ierr, io
      integer                 :: secs, days
      integer                 :: dum
  
!-------------------------------------------------------------------
!  local variables:
!
!     unit                unit number for nml file
!     ierr                error return flag
!     io                  error return code
!     secs                seconds component of time_type variable Time
!     days                days component of time_type variable Time
!                         
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
      if (module_is_initialized) return
      if (.not. module_is_initialized .and. &
          .not. present(tracers_in_donner)) then  
        call error_mesg ('donner_deep_mod', &
           ' must have first call to donner_deep_init provide '// &
            'tracers_in_donner as argument', FATAL)
      endif
      

!---------------------------------------------------------------------
!    verify that all modules used by this module are initialized.
!---------------------------------------------------------------------
      call fms_init
      call constants_init
      call diag_manager_init
      call field_manager_init (dum)
!  not yet existent
!     call tracer_manager_init
      call time_manager_init
      call column_diagnostics_init (lonb, latb)
      call sat_vapor_pres_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=donner_deep_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'donner_deep_nml')
        enddo
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() )    &
                                 write (stdlog(), nml=donner_deep_nml)

!-------------------------------------------------------------------
!    define the grid dimensions. idf and jdf are the (i,j) dimensions of
!    the domain on this processor, nlev is the number of model layers.
!-------------------------------------------------------------------
      nlev = size(pref(:)) - 1
      idf  = size(lonb(:)) - 1
      jdf  = size(latb(:)) - 1

!---------------------------------------------------------------------
!    test that donner_deep_freq has a valid value.
!---------------------------------------------------------------------
      if (donner_deep_freq > 86400) then
        call error_mesg ( 'donner_deep_mod', &
         ' donner convection must be called at least once per day', &
                                                           FATAL)
      else if (donner_deep_freq <= 0) then
        call error_mesg ( 'donner_deep_mod', &
          ' a positive value must be assigned to donner_deep_freq', &
                                                            FATAL)
      endif

!---------------------------------------------------------------------
!    test that cell_liquid_size_type has been validly specified.
!---------------------------------------------------------------------
      if (trim(cell_liquid_size_type) == 'input') then
        do_input_cell_liquid_size = .true.
      else if (trim(cell_liquid_size_type) == 'bower') then
        do_bower_cell_liquid_size = .true.
      else
        call error_mesg ( 'donner_deep_mod', &
           ' cell_liquid_size_type must be input or bower', &
                                                  FATAL)
      endif

!---------------------------------------------------------------------
!    test that cell_ice_size_type has been validly specified, and if
!    specified as 'input', that cell_ice_geneff_diam_input has also 
!    been defined.
!---------------------------------------------------------------------
      if (trim(cell_ice_size_type) == 'input') then
        do_input_cell_ice_size = .true.
        if (cell_ice_geneff_diam_input <= 0.0) then
          call error_mesg ('donner_deep_mod', &
               ' must define a nonnegative generalized effective '//&
       'diameter for ice when cell_ice_size_type is input', &
                                                           FATAL)
        endif
      else if (trim(cell_ice_size_type) == 'default') then
        do_default_cell_ice_size = .true.
      else
        call error_mesg ( 'donner_deep_init', &
             ' cell_ice_size_type must be input or default', &
                                                     FATAL)
      endif

!---------------------------------------------------------------------
!    if column diagnostics are desired from this module, call 
!    donner_column_init to do necessary initialization.
!---------------------------------------------------------------------
      num_diag_pts = num_diag_pts_ij + num_diag_pts_latlon
      if (num_diag_pts > 0) then
        call donner_column_init (pref, Time)
      endif

!---------------------------------------------------------------------
!    if tracers_in_donner has been passed in, determine if and how many 
!    tracers are to be transported by donner_deep convection. if not
!    present, no tracers are transported (default).
!---------------------------------------------------------------------
      ncont = count(tracers_in_donner)
      if (ncont > 0) then
        do_donner_tracer = .true.
      else
        do_donner_tracer = .false.
      endif
      
!--------------------------------------------------------------------
!    allocate module variables that will be saved across timesteps.
!--------------------------------------------------------------------
      call allocate_variables

!--------------------------------------------------------------------
!    activate the netcdf diagnostic fields.
!-------------------------------------------------------------------
      call register_fields (Time, axes, tracers_in_donner)

!--------------------------------------------------------------------
!    if a restart file is present, call read_restart_donner_deep to 
!    read it.
!--------------------------------------------------------------------
      if (file_exist ('INPUT/donner_deep.res') ) then
        call read_restart_donner_deep 

!--------------------------------------------------------------------
!    if a restart file is not present, define the time until the calcul-
!    ation of convection by donner_deep_mod. it is given
!    by the time remaining until the next integral multiple of 
!    donner_deep_freq from the start of the day. set a flag to indicate 
!    that this module is being coldstarted, and call 
!    initialize_variables to provide values for the module variables 
!    until that first call is made.
!--------------------------------------------------------------------
      else
        call get_time (Time, secs, days)
        if (secs == 0) then    ! i.e., 00Z
          conv_alarm     = donner_deep_freq
        else 
          conv_alarm     = donner_deep_freq -   &
                    MOD (secs, donner_deep_freq)
        endif
        coldstart = .true.
        call initialize_variables
      endif

!---------------------------------------------------------------------
!    initialize the points processed counters. define the total number 
!    of columns present on the processor. 
!---------------------------------------------------------------------
      pts_processed_conv   = 0
      pts_processed_cape   = 0
      total_pts = idf*jdf

!--------------------------------------------------------------------
!    set flag to indicate that donner_deep_mod has been initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------



end subroutine donner_deep_init



!###################################################################

subroutine donner_deep (is, ie, js, je, temp, mixing_ratio, pfull,  &
                        phalf, omega, dt, land, Time, &
                        ttnd, qtnd, precip, ahuco, qrat,  &
                        kbot, cf, qlin, qiin, &
                        delta_qa, delta_ql, delta_qi, mtot, &
                        tracers, qtrceme)
                        
!-------------------------------------------------------------------
!   donner_deep is the prognostic subroutine of donner_deep_mod. it 
!   returns the tendencies of temperature and mixing ratio (and changes
!   in cloudwater, cloudice and cloud area when run with a prognostic 
!   cloud scheme), and the precipitation and the mass flux (if needed)
!   produced by deep convection, as determined by this parameterization.
!-------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                 intent(in)            :: is, ie, js, je
real, dimension(:,:,:),  intent(in)            :: temp, mixing_ratio
real, dimension(:,:,:),  intent(in)            :: pfull, phalf
real, dimension(:,:,:),  intent(in)            :: omega
real,                    intent(in)            :: dt
real, dimension(:,:),    intent(in)            :: land
type(time_type),         intent(in)            :: Time
real, dimension(:,:,:),  intent(out)           :: ttnd, qtnd
real, dimension(:,:),    intent(out)           :: precip      
real, dimension(:,:,:),  intent(inout)         :: ahuco, qrat
real, dimension(:,:,:,:), intent(in), optional :: tracers 
real, dimension(:,:,:,:), intent(out),optional :: qtrceme
integer, dimension(:,:), intent(in),  optional :: kbot
real, dimension(:,:,:),  intent(in),  optional :: cf
real, dimension(:,:,:),  intent(in),  optional :: qlin,qiin
real, dimension(:,:,:),  intent(out), optional :: delta_qa, delta_ql, &
                                                  delta_qi
real, dimension(:,:,:),  intent(out), optional :: mtot
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   mixing ratio field at model levels [ kg/kg ]
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field at half-levels 1:nlev+1  [ Pa ]
!     omega          model omega field [ Pa/sec ]
!     dt             physics time step [ sec ]
!     land           fraction of surface (grid box) covered by land
!                    [ nondimensional ]
!     Time           current time (time_type)
!
!   intent(out) variables:
!
!     ttnd           time tendency of temperature due to deep convection
!                    [ deg K / sec ]
!     qtnd           time tendency of mixing ratio due to deep
!                    convection [ (kg/kg) / sec ]
!     precip         precipitation generated by deep convection
!                    [ kg / m**2 ]
!     ahuco          cell + meso fraction for specific humidity
!                    (mesoscale downdraft does not contain cloud
!                    water, so this fraction should not be used
!                    for radiation)
!                    index 1 nearest ground
!     qrat           ratio of large-scale specific humidity to
!                    specific humidity in environment outside convective
!                    system
!                    index 1 nearest ground
!
!   intent(in), optional variables:
!
!     kbot           index of lowest model level (only needed for eta
!                    coordinate case)

!     These variables are present when a prognostic cloud scheme is 
!     employed:

!     cf             large-scale cloud fraction  
!                    [ nondimensional ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(liquid) / kg ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(ice) / kg ]
!
!     This variable is present when tracers are transported by donner
!     convection:
!
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg(air) ]
!
!   intent(out), optional variables:
!
!     These variables are present when a prognostic cloud scheme is 
!     employed:
!
!     delta_qa       change in cloud area due to deep convection
!                    during the time step [ dimensionless ]
!     delta_ql       change in cloud water due to deep convection 
!                    during the time step [ (kg/kg) ]
!     delta_qi       change in cloud ice due to deep convection 
!                    during the time step [ (kg/kg) ]
!     mtot           mass flux at model levels, convective plus meso-
!                    scale, due to donner_deep_mod [ (kg/m**2) / sec ]
!
!     This variable is present when tracers are transported by donner
!     convection:
!
!     qtrceme        tracer time tendencies due to donner_deep_mod
!                    [ kg / kg / sec ]
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!   local variables:
      real, dimension (size(temp,1), size(temp,2))  ::  omint
      real, dimension (size(temp,1), size(temp,2),              &
                                     size(temp,3))  ::  &
                                     cape_input_temp, cape_input_vapor,&
                                     pmass
      integer    :: idt
      integer    :: i, j, k        
      integer    :: itest,jtest,ktest,kcont

!--------------------------------------------------------------------
!   local variables:
!
!     omint                low-level parcel displacement to use in cape
!                          calculation on this step [ Pa ]
!     cape_input_temp      temperature profile to use in cape  cal-
!                          culation [ deg K ]
!     cape_input_vapor     mixing ratio profile to use in cape cal-
!                          culation [ kg/kg ]
!     idt                  physics time step [ sec ]
!     i, j                 loop indices
!
!---------------------------------------------------------------------


!-------------------------------------------------------------------
!    save the model timestep as an integer. verify that donner_deep_freq
!    is an integral multiple of it.
!--------------------------------------------------------------------
      idt = int (dt)
      if (MOD (donner_deep_freq, idt) /= 0) then
        call error_mesg ('donner_deep_time_vary',  &
          'donner_deep timestep NOT an integral multiple of '//&
           'physics timestep', FATAL)
       endif

!--------------------------------------------------------------------
!    on the first call to this subroutine in the job, initialize 
!    cape_alarm to be one physics step prior to first convection calcul-
!    ation, so that the time tendency of cape may be available when 
!    convection is calculated. after the first call, the alarm
!    will be incremented by the donner calling frequency, so that it 
!    will always be called one step prior to the convection call.
!--------------------------------------------------------------------
      if (first_call) then
        cape_alarm = conv_alarm - idt
        first_call = .false.
      endif

!--------------------------------------------------------------------
!    decrement the time remaining before the cape and convection calcu-
!    lations on the first entry to this routine on a given timestep.
!--------------------------------------------------------------------
      if (js == 1) then
        cape_alarm  = cape_alarm - idt
        conv_alarm  = conv_alarm - idt
      endif

!--------------------------------------------------------------------
!    set flags to indicate whether the convection and/or cape calcul-
!    ations are to be done on this timestep. if this is the first call
!    to donner_deep (i.e., coldstart), convection cannot be calculated.
!    otherwise, it is or not dependent on whether the convection 
!    "alarm" has gone off. the calculation or non-calculation of cape
!    is controlled by the cape alarm.
!---------------------------------------------------------------------
      if (coldstart) then
        conv_calc_on_this_step = .false.
      else
        if (conv_alarm <= 0) then
          conv_calc_on_this_step = .true.
        else
          conv_calc_on_this_step = .false.
        endif
      endif
      if (cape_alarm <= 0) then   
        conv_calc_on_next_step    = .true.
      else
        conv_calc_on_next_step    = .false.
      endif

!-------------------------------------------------------------------
!    define the physics window size.
!-------------------------------------------------------------------
      isize = ie - is + 1
      jsize = je - js + 1

!-------------------------------------------------------------------
!    call initialize_local_variables to allocate and initialize the
!    elements of the donner_conv and donner_cape derived type 
!    variables.
!-------------------------------------------------------------------
      call initialize_local_variables (Don_conv, Don_cape)

!---------------------------------------------------------------------
!    calculate time integrated low-level displacement (Cu Closure E 
!    notes, 8 feb 98). omint_acc is the time integrated vertical dis-
!    placement; omint is the value to be used on this step, which 
!    is the time integrated value, unless the current displacement is 
!    downward. 
!---------------------------------------------------------------------
      do j=1,jsize       
        do i=1,isize        
          omint_acc(i+is-1,j+js-1) = omint_acc(i+is-1,j+js-1) +   &
                                omega(i,j, nlev)*dt
          omint_acc(i+is-1,j+js-1) = MIN (0.0, omint_acc(i+is-1,j+js-1))
          omint_acc(i+is-1,j+js-1) = MAX (omint_acc(i+is-1,j+js-1),  &
                                     -phalf(i,j,nlev+1))
          if (omega(i,j,nlev) > 0.)   then
            omint(i,j) = 0. 
          else
            omint(i,j) = omint_acc(i+is-1,j+js-1)
          endif
        end do
      end do

!--------------------------------------------------------------------
!    initialize column diagnostics control in this subdomain. 
!-------------------------------------------------------------------
      call donner_column_control (is, js, Time)

!----------------------------------------------------------------------
!    call donner_column_input_fields to print out input fields, 
!    location and control information for diagnostics columns.   
!----------------------------------------------------------------------
      if (in_diagnostics_window) then
        call donner_column_input_fields (dt, conv_calc_on_this_step, &
                                         temp, mixing_ratio, phalf, &
                                         omega) 
      endif

!---------------------------------------------------------------------
!    call donner_convection_driver to calculate the effects of deep 
!    convection, on timesteps where that calculation is desired.
!---------------------------------------------------------------------
      if (conv_calc_on_this_step) then
        if (present (qlin)) then
          if (present(tracers)) then
            call donner_convection_driver (is, ie, js, je, temp, &
                                           mixing_ratio, phalf, pfull,&
                                           dt, omint, land,  &
                                           Don_cape, Don_conv, &
                                           qrat, ahuco,   &
                                           tracers=tracers, &
                                           qlin=qlin, qiin=qiin, cf=cf)
          else
            call donner_convection_driver (is, ie, js, je, temp, &
                                           mixing_ratio, phalf, pfull,&
                                           dt, omint, land,  &
                                           Don_cape, Don_conv, &
                                           qrat, ahuco, &
                                           qlin=qlin, qiin=qiin, cf=cf)
          endif
        else
          if (present(tracers)) then
            call donner_convection_driver (is, ie, js, je, temp, &
                                           mixing_ratio, phalf, pfull,&
                                           dt, omint, land, &
                                           Don_cape, Don_conv, &
                                           qrat, ahuco,  &
                                           tracers=tracers)
          else
            call donner_convection_driver (is, ie, js, je, temp, &
                                           mixing_ratio, phalf, pfull,&
                                           dt, omint, land, &
                                           Don_cape, Don_conv, &
                                           qrat, ahuco)
          endif

!--------------------------------------------------------------------
!    save tracer tendencies for access on other steps when donner_deep
!    not activated every step.
!--------------------------------------------------------------------
        endif
        if (do_donner_tracer) then
            tracer_tends(is:ie,js:je,:,:) = Don_conv%qtceme(:,:,:,:)
        endif
      endif ! (conv_calc_on_this_step)

!----------------------------------------------------------------------
!    define tendencies from donner deep to be passed back to calling
!    routine. these values are supplied even on non-calculation steps.
!    when a prognostic cloud scheme is active, then deltas to cloud
!    water, cloud ice, cloud area and the convective mass flux are also
!    returned. when tracers are being transported by donner convection,
!    their time tendencies are returned.
!----------------------------------------------------------------------
      if (present (qlin)) then
        if (present (qtrceme)) then
          call define_output_fields (is, ie, js, je, dt, mixing_ratio, &
                                   ttnd, qtnd, precip,  &
                                   Don_conv, &
                                   mtot,  &
                                 delta_ql, delta_qi, delta_qa, &
                                 qtrceme=qtrceme)
        else
          call define_output_fields (is, ie, js, je, dt, mixing_ratio, &
                                   ttnd, qtnd, precip, &
                                   Don_conv,  mtot,  &
                                   delta_ql, delta_qi, delta_qa)
        endif
      else
        if (present (qtrceme)) then
          call define_output_fields (is, ie, js, je, dt, mixing_ratio, &
                                   ttnd, qtnd, precip, &
                                   Don_conv, qtrceme=qtrceme)
        else
          call define_output_fields (is, ie, js, je, dt, mixing_ratio, &
                                   ttnd, qtnd, precip, &
                                   Don_conv)
        endif
      endif

!---------------------------------------------------------------------
!    the following if loop is executed if the calculation of cape is
!    needed on this step so that a cape time tendency may be calculated
!    when donner_deep convection is calculated on the next step.
!---------------------------------------------------------------------
      if (conv_calc_on_next_step) then
        call cape_calculation_driver (is, ie, js, je, temp,   &
                                      mixing_ratio, ttnd, qtnd, &
                                      pfull, Don_cape)
      endif  ! (conv_calc_on_next_step   )
 
 !-------------------------------------------------------------------
 !  on convection calculation steps, print out desired column diagnos-
 !  tics, and also send the desired diagnostic data to the diagnostics 
 !  handler module so that netcdf output may be produced.
 !-------------------------------------------------------------------
      if (conv_calc_on_this_step) then
        if (in_diagnostics_window) then
          call donner_column_end_of_step (Don_conv, Don_cape)
        endif
        do k=1,nlev
!! value used for gravity makes some difference here, since module
!! not yet unified value for grav . using here value from mulsub_vect.
!         pmass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/GRAV   
          pmass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/9.80616
        end do
        call donner_deep_netcdf (is, ie, js, je, Time, Don_conv,  &
                                 Don_cape, pmass)
      endif  ! (conv_calc_on_this_step)

!--------------------------------------------------------------------
!   call deallocate_local_variables to deallocate space used by the
!   derived-type variables.
!--------------------------------------------------------------------
      call deallocate_local_variables (Don_conv, Don_cape)

!--------------------------------------------------------------------



end subroutine donner_deep



!####################################################################

subroutine donner_deep_avg (is, ie, js, je,   &
                            cell_cloud_frac_out,  cell_liquid_amt_out,&
                            cell_liquid_size_out, cell_ice_amt_out,   &
                            cell_ice_size_out,    meso_cloud_frac_out,&
                            meso_liquid_amt_out,  meso_liquid_size_out,&
                            meso_ice_amt_out,     meso_ice_size_out)

!------------------------------------------------------------------
!   this subroutine provides the cloud microphysical quantities assoc-
!   iated with donner_deep convection to the radiation package. these
!   include the cloud liquid and ice amounts, liquid and ice sizes,
!   and fractional coverage, for both the convective cell and mesoscale
!   components of the convective system.
!------------------------------------------------------------------

integer, intent(in)                    :: is, ie, js, je
real,    intent(out), dimension(:,:,:) :: cell_cloud_frac_out,        &
                                          cell_liquid_amt_out,   &
                                          cell_liquid_size_out,  &
                                          cell_ice_amt_out, &
                                          cell_ice_size_out,   &
                                          meso_cloud_frac_out,   &
                                          meso_liquid_amt_out, &
                                          meso_liquid_size_out, &
                                          meso_ice_amt_out,  &
                                          meso_ice_size_out

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!
!  intent(out) variables:
!
!     cell_cloud_frac_out  fractional coverage of convective cells in 
!                          grid box [ dimensionless ]
!     cell_liquid_amt_out  liquid water content of convective cells 
!                          [ kg(h2o)/kg(air) ]
!     cell_liquid_size_out assumed effective size of cell liquid drops
!                          [ micrometers ]
!     cell_ice_amt_out     ice water content of cells 
!                          [ kg(h2o)/kg(air) ]
!     cell_ice_size_out    generalized effective diameter for ice in
!                          convective cells [ micrometers ]
!     meso_cloud_frac_out  fractional area of mesoscale clouds in grid 
!                          box [ dimensionless ]
!     meso_liquid_amt_out  liquid water content in mesoscale clouds
!                          currently set to 0.0
!                          [ kg(h2o)/kg(air) ]
!     meso_liquid_size_out assumed effective size of mesoscale drops
!                          currently set to 0.0 [ micrometers ]
!     meso_ice_amt_out     ice water content of mesoscale elements 
!                          [ kg(h2o)/kg(air) ]
!     meso_ice_size_out    generalized ice effective size for anvil ice
!                          [ micrometers ]
!
!---------------------------------------------------------------------

   
!------------------------------------------------------------------
!   local variables

      real, dimension(size(cell_cloud_frac_out,1), &
                      size(cell_cloud_frac_out,2))   :: inv_nsum
      integer                                        :: num, k
   
!---------------------------------------------------------------------
!   local variables:
!
!       num           number of grid columns which have not been given
!                     values for the output variables by donner_deep_mod
!       k             do loop index
!       inv_sum       inverse of number of elements in the time averaged
!                     output fields
!
!--------------------------------------------------------------------- 

!---------------------------------------------------------------------
!   check to make sure dimensions of arguments match the module
!   variable dimensions.
!---------------------------------------------------------------------
      if (size(cell_cloud_frac_out,3) /= size(cell_cloud_frac,3)) &
        call error_mesg (  'donner_deep_mod',  &
                         'input argument has the wrong size',FATAL)

!---------------------------------------------------------------------
!    check to see that all columns have been given values for the module
!    variables that are going to be averaged. 
!----------------------------------------------------------------------
      num = count(nsum(is:ie,js:je) == 0)
!----------------------------------------------------------------------
!    if any columns have not been given values, stop execution with an 
!    error message. 
!----------------------------------------------------------------------
      if (num > 0) then
        call error_mesg ( 'donner_deep_mod', &
                         'nsum has some zero entries', FATAL)

!----------------------------------------------------------------------
!    if all columns have valid data, produce time averaged values of
!    the desired output fields.
!----------------------------------------------------------------------
      else
        inv_nsum(:,:) = 1./float(nsum(is:ie,js:je))
        do k=1,size(cell_cloud_frac_out,3)
          cell_cloud_frac_out(:,:,k) =   &
                              cell_cloud_frac(is:ie,js:je,k)*inv_nsum
          cell_liquid_amt_out(:,:,k) =   &
                              cell_liquid_amt(is:ie,js:je,k)*inv_nsum
          cell_liquid_size_out(:,:,k) =  &
                              cell_liquid_size(is:ie,js:je,k)*inv_nsum
          cell_ice_amt_out(:,:,k) =      &
                              cell_ice_amt(is:ie,js:je,k)*inv_nsum
          cell_ice_size_out(:,:,k) =     &
                              cell_ice_size(is:ie,js:je,k)*inv_nsum
          meso_cloud_frac_out(:,:,k) =   &
                               meso_cloud_frac(is:ie,js:je,k)*inv_nsum
          meso_liquid_amt_out(:,:,k) =   &
                               meso_liquid_amt(is:ie,js:je,k)*inv_nsum
          meso_liquid_size_out(:,:,k) =  &
                               meso_liquid_size(is:ie,js:je,k)*inv_nsum
          meso_ice_amt_out(:,:,k) =      &
                               meso_ice_amt(is:ie,js:je,k)*inv_nsum
          meso_ice_size_out(:,:,k) =     & 
                               meso_ice_size(is:ie,js:je,k)*inv_nsum
        end do
     
!----------------------------------------------------------------------
!    checks to eliminate unreasonable values and incompatibilities 
!    between output fields.
!----------------------------------------------------------------------
        where (cell_ice_size_out .ge. 13.3 .and. &
               cell_ice_size_out .le. 18.6)
          cell_ice_size_out   =  18.6
        end where
        where (cell_cloud_frac_out  > 0.0 .and.  &
               cell_liquid_amt_out == 0.0 .and.  &
               cell_ice_amt_out    == 0.0)
          cell_cloud_frac_out  = 0.0
        end where
        where (cell_cloud_frac_out == 0.0 .and.   &
               cell_liquid_amt_out  > 0.0)
          cell_liquid_amt_out  = 0.0
        end where
        where (cell_cloud_frac_out == 0.0 .and.   &
               cell_ice_amt_out     > 0.0)
          cell_ice_amt_out     = 0.0
        end where

        where (meso_ice_size_out .ge. 13.3 .and. &
               meso_ice_size_out .le. 18.6)
          meso_ice_size_out   =  18.6
        end where
        where (meso_cloud_frac_out  > 0.0 .and.  &
               meso_liquid_amt_out == 0.0 .and.  &
               meso_ice_amt_out    == 0.0)
          meso_cloud_frac_out  = 0.0
        end where
        where (meso_cloud_frac_out == 0.0 .and.   &
               meso_liquid_amt_out  > 0.0)
          meso_liquid_amt_out  = 0.0
        end where
        where (meso_cloud_frac_out == 0.0 .and.   &
               meso_ice_amt_out     > 0.0)
          meso_ice_amt_out     = 0.0
        end where
      endif
 
!----------------------------------------------------------------------
!    reset the variables just processed so that new sums may be begun 
!    when donner_deep is called again.
!----------------------------------------------------------------------
      nsum            (is:ie,js:je)   = 0
      cell_cloud_frac (is:ie,js:je,:) = 0.0
      cell_liquid_amt (is:ie,js:je,:) = 0.0
      cell_liquid_size(is:ie,js:je,:) = 0.0
      cell_ice_amt    (is:ie,js:je,:) = 0.0
      cell_ice_size   (is:ie,js:je,:) = 0.0
      meso_cloud_frac (is:ie,js:je,:) = 0.0
      meso_liquid_amt (is:ie,js:je,:) = 0.0
      meso_liquid_size(is:ie,js:je,:) = 0.0
      meso_ice_amt    (is:ie,js:je,:) = 0.0
      meso_ice_size   (is:ie,js:je,:) = 0.0
       
!----------------------------------------------------------------------



end subroutine donner_deep_avg




!##################################################################


subroutine donner_deep_end

!---------------------------------------------------------------------
!   this is the destructor for donner_deep_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variable

      integer :: unit          ! unit number for restart file
      integer :: n             ! do-loop index

!-------------------------------------------------------------------
!    open unit for restart file.
!-------------------------------------------------------------------
       unit = open_restart_file ('RESTART/donner_deep.res', 'write')

!-------------------------------------------------------------------
!    file writing is currently single-threaded. write out restart
!    version, time remaining until next call to donner_deep_mod and
!    the freaquency of calculating donner_deep convection.
!-------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe()) then
        write (unit) restart_versions(size(restart_versions(:)))
        write (unit) conv_alarm, donner_deep_freq
      endif

!-------------------------------------------------------------------
!    write out the donner_deep restart variables.
!    cemetf    - heating rate due to donner_deep
!    cememf    - moistening rate due to donner_deep
!    xcape_lag - cape value which will be used on next step in
!                calculation od dcape/dt
!-------------------------------------------------------------------
      call write_data (unit, cemetf)
      call write_data (unit, cememf)
      call write_data (unit, xcape_lag)
      
!----------------------------------------------------------------------
!    the following variables are needed for the current cape closure.
!----------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe()) then
        write (unit) model_levels_in_sfcbl
      endif
      call write_data (unit, tempbl)
      call write_data (unit, ratpbl)

!--------------------------------------------------------------------
!    the following variables are needed when a prognostic cloud scheme
!    is being used. they are always present in the restart file, having
!    been initialized to zero, if prognostic clouds are not active.
!--------------------------------------------------------------------
      call write_data (unit, mass_flux)
      call write_data (unit, delta_ql )
      call write_data (unit, delta_qi )
      call write_data (unit, delta_qa )

!----------------------------------------------------------------------
!    
!-------------------------------------------------------------------
!    write out more donner_deep restart variables.
!    qint_lag   - column integrated water vapor mixing ratio
!    omint_acc  - time-integrated low-level vertical displacement
!    tprea1     - precipitation due to donner_deep_mod
!----------------------------------------------------------------------
      call write_data (unit, qint_lag )
      call write_data (unit, omint_acc)
      call write_data (unit, tprea1)

!---------------------------------------------------------------------
!    the following variables contain cloud property information needed 
!    by the radiation package:
!---------------------------------------------------------------------
      call write_data(unit, cell_cloud_frac  )
      call write_data(unit, cell_liquid_amt  )
      call write_data(unit, cell_liquid_size )
      call write_data(unit, cell_ice_amt     )
      call write_data(unit, cell_ice_size    )
 
      call write_data(unit, meso_cloud_frac  )
      call write_data(unit, meso_liquid_amt  )
      call write_data(unit, meso_liquid_size )
      call write_data(unit, meso_ice_amt     )
      call write_data(unit, meso_ice_size    )
      call write_data(unit, nsum             )

!---------------------------------------------------------------------
!    write out the number of tracers that are being transported by
!    donner_deep_mod.
!---------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe()) then
        write (unit) ncont
      endif

!----------------------------------------------------------------------
!    if tracers are being transported, write out their names and current!    time tendencies.
!----------------------------------------------------------------------
      if (do_donner_tracer) then
        do n=1,ncont
          if (mpp_pe() == mpp_root_pe()) then
            write (unit) tracername(n)         
          endif
          call write_data(unit, tracer_tends(:,:,:,n))
        end do
      endif

!-------------------------------------------------------------------
!    close restart file unit.
!------------------------------------------------------------------
      call close_file (unit)

!-------------------------------------------------------------------
!    close any column diagnostics units which are open.
!------------------------------------------------------------------
      if (num_diag_pts > 0) then
        call close_column_diagnostics_units (col_diag_unit)
      endif

!---------------------------------------------------------------------
      module_is_initialized = .false.

 
end subroutine donner_deep_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      1. ROUTINES CALLED BY DONNER_DEEP_INIT
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
 

subroutine register_fields (Time, axes, tracers_in_donner)

!----------------------------------------------------------------------

type(time_type),               intent(in)           :: Time
integer,         dimension(4), intent(in)           :: axes
logical,         dimension(:), intent(in)         :: tracers_in_donner
!----------------------------------------------------------------------

       integer :: n, num_tracers, nn
       character(len=32) :: tracer_units


            id_cemetf_deep =  &
           register_diag_field (mod_name, 'cemetf_deep', axes(1:3),   &
           Time, 'heating rate, c + m ', 'K/s',   &
           missing_value=missing_value)
            id_ceefc_deep =  &
           register_diag_field (mod_name, 'ceefc_deep', axes(1:3),   &
           Time, 'cell entrpy flx cnvrgnc', 'K/s',   &
           missing_value=missing_value)
            id_cecon_deep =  &
           register_diag_field (mod_name, 'cecon_deep', axes(1:3),   &
           Time, 'cell cond/evap ', 'K/s',   &
           missing_value=missing_value)
            id_cemfc_deep =  &
           register_diag_field (mod_name, 'cemfc_deep', axes(1:3),   &
           Time, 'cell moist flx cnvgnc', 'kg(h2o)/kg/s',   &
           missing_value=missing_value)
            id_cememf_deep =  &
           register_diag_field (mod_name, 'cememf_deep', axes(1:3),   &
           Time, 'moistening rate, c + m ', 'kg(h2o)/kg/s',   &
           missing_value=missing_value)
            id_cememf_mod_deep =  &
           register_diag_field (mod_name, 'cememf_mod_deep', axes(1:3),&
           Time, 'mod cememf due to negative q ', 'kg(h2o)/kg/s',   &
           missing_value=missing_value)
            id_cual_deep =  &
           register_diag_field (mod_name, 'cual_deep', axes(1:3),   &
           Time, 'c + m cld frac ', 'percent',   &
           missing_value=missing_value)
            id_fre_deep =  &
           register_diag_field (mod_name, 'fre_deep', axes(1:3),   &
           Time, 'freezing ', 'K/sec',   &
           missing_value=missing_value)
            id_elt_deep =  &
           register_diag_field (mod_name, 'elt_deep', axes(1:3),   &
           Time, 'melting', 'K/sec',   &
           missing_value=missing_value)
            id_cmus_deep =  &
           register_diag_field (mod_name, 'cmus_deep', axes(1:3),   &
           Time, 'meso-up deposition', 'kg(h2o)/kg/sec)',   &
           missing_value=missing_value)
            id_ecds_deep =  &
           register_diag_field (mod_name, 'ecds_deep', axes(1:3),   &
           Time, 'convective dwndrft evap ', 'kg(h2o)/kg/sec', &
           missing_value=missing_value)
            id_eces_deep =  &
           register_diag_field (mod_name, 'eces_deep', axes(1:3),   &
           Time, 'convective updrft evap/subl ',    &
           'kg(h2o)/kg/sec',   missing_value=missing_value)
            id_emds_deep =  &
           register_diag_field (mod_name, 'emds_deep', axes(1:3),   &
           Time, 'meso-dwn subl ', 'kg(h2o)/kg/sec',   &
           missing_value=missing_value)
            id_emes_deep =  &
           register_diag_field (mod_name, 'emes_deep', axes(1:3),   &
           Time, 'meso-up subl ', 'kg(h2o)/kg/sec',   &
           missing_value=missing_value)
            id_qmes_deep =  &
           register_diag_field (mod_name, 'qmes_deep', axes(1:3),   &
           Time, 'meso moist flux conv', 'kg(h2o)/kg/sec',   &
           missing_value=missing_value)
            id_wmps_deep =  &
           register_diag_field (mod_name, 'wmps_deep', axes(1:3),   &
           Time, 'meso redistrib of vapor from cells',    &
           'kg(h2o)/kg/sec', missing_value=missing_value)
            id_wmms_deep =  &
           register_diag_field (mod_name, 'wmms_deep', axes(1:3),   &
           Time, 'meso depo of vapor from cells',    &
           'kg(h2o)/kg/sec',  missing_value=missing_value)
            id_tmes_deep =  &
           register_diag_field (mod_name, 'tmes_deep', axes(1:3),   &
           Time, 'meso entropy flux conv',  'K/sec',   &
           missing_value=missing_value)
            id_dmeml_deep =  &
           register_diag_field (mod_name, 'dmeml_deep', axes(1:3), &
           Time, 'mass flux meso dwndrfts', 'kg/((m**2) s)',   &
           missing_value=missing_value)
            id_uceml_deep =  &
           register_diag_field (mod_name, 'uceml_deep', axes(1:3), &
           Time, 'mass flux cell updrfts', 'kg/((m**2) s)',   &
           missing_value=missing_value)
            id_umeml_deep =  &
           register_diag_field (mod_name, 'umeml_deep', axes(1:3), &
           Time, 'mass flux meso updrfts', 'kg/((m**2) s)',   &
           missing_value=missing_value)
            id_xice_deep =  &
           register_diag_field (mod_name, 'xice_deep', axes(1:3),  &
           Time, 'meso ice mass mixing ratio ', 'kg(ice)/kg',   &
           missing_value=missing_value)

!---------------------------------------------------------------------
!    if tracers are being transported by donner_deep_mod, define 
!    appropriate diagnostics.
!---------------------------------------------------------------------
             num_tracers = size(tracers_in_donner(:))
             allocate (id_qtren1 (ncont))
             allocate (id_qtmes1 (ncont))
             allocate (id_wtp1   (ncont))
             allocate (id_qtceme (ncont))
             allocate (tracername(ncont))
             allocate (id_qtren1_col (ncont))
             allocate (id_qtmes1_col (ncont))
             allocate (id_wtp1_col   (ncont))
             allocate (id_qtceme_col (ncont))
             nn = 1
             do n=1,num_tracers
               if (tracers_in_donner(n)) then
                 call get_tracer_names (MODEL_ATMOS, n,  &
                                        name = tracername(nn), &
                                        units = tracer_units)

                 id_qtren1(nn) =  &
                       register_diag_field (mod_name,  &
                            trim(tracername(nn)) // '_qtren1',  &
                            axes(1:3), Time,  &
                            trim(tracername(nn)) // ' cell tendency ', &
                            trim(tracer_units)//'/s', &
                                       missing_value=missing_value)
                 id_qtmes1(nn) =  &
                       register_diag_field (mod_name,  &
                            trim(tracername(nn)) // '_qtmes1',  &
                            axes(1:3), Time,   &
                          trim(tracername(nn)) //' mesoscale tendency',&
                            trim(tracer_units)//'/s', &
                                       missing_value=missing_value)
                 id_wtp1(nn) =  &
                      register_diag_field (mod_name,  &
                          trim(tracername(nn)) // '_wtp1', axes(1:3), &
                          Time,  &
                          trim(tracername(nn)) //' mesoscale redist',&
                            trim(tracer_units)//'/s', &
                                     missing_value=missing_value)
                 id_qtceme(nn) =  &
                       register_diag_field (mod_name,  &
                            trim(tracername(nn)) // '_qtceme',  &
                            axes(1:3), Time,  &
                            trim(tracername(nn)) // ' total tendency ',&
                            trim(tracer_units)//'/s', &
                                       missing_value=missing_value)

                 id_qtren1_col(nn) =  &
                       register_diag_field (mod_name,  &
                            trim(tracername(nn)) // '_qtren1_col',  &
                            axes(1:2), Time,  &
                     'column integrated ' //trim(tracername(nn)) //  &
                                                   ' cell tendency ', &
                            trim(tracer_units)//'/s', &
                                       missing_value=missing_value)
                 id_qtmes1_col(nn) =  &
                       register_diag_field (mod_name,  &
                            trim(tracername(nn)) // '_qtmes1_col',  &
                            axes(1:2), Time,   &
               'column integrated ' //trim(tracername(nn)) //  &
                                             ' mesoscale tendency',&
                            trim(tracer_units)//'/s', &
                                       missing_value=missing_value)
                 id_wtp1_col(nn) =  &
                      register_diag_field (mod_name,  &
                          trim(tracername(nn)) // '_wtp1_col', axes(1:2), &
                          Time,  &
                  'column integrated '//trim(tracername(nn)) //  &
                                            ' mesoscale redist',&
                            trim(tracer_units)//'/s', &
                                     missing_value=missing_value)
                 id_qtceme_col(nn) =  &
                       register_diag_field (mod_name,  &
                            trim(tracername(nn)) // '_qtceme_col',  &
                            axes(1:2), Time,  &
                   'column integrated ' //trim(tracername(nn)) // &
                                                 ' total tendency ', &
                            trim(tracer_units)//'/s', &
                                       missing_value=missing_value)
                 nn = nn + 1
               endif
             end do

            id_dgeice_deep =  &
           register_diag_field (mod_name, 'dgeice_deep', axes(1:3), &
           Time, 'meso ice gen eff size ', 'micrometers',   &
           missing_value=missing_value)
            id_cuqi_deep =  &
           register_diag_field (mod_name, 'cuqi_deep', axes(1:3),  &
           Time, 'cell ice ', 'kg(H2O)/kg',   &
           missing_value=missing_value)
            id_cuql_deep =  &
           register_diag_field (mod_name, 'cuql_deep', axes(1:3),  &
           Time, 'cell liquid ', 'kg(H2O)/kg',   &
           missing_value=missing_value)
            id_dgeliq_deep =  &
           register_diag_field (mod_name, 'dgeliq_deep', axes(1:3), &
            Time, 'cell liq gen eff size ', 'micrometers',   &
            missing_value=missing_value)

            id_plcl_deep =  &
           register_diag_field (mod_name, 'plcl_deep', axes(1:2),   &
           Time, 'pressure at lcl ', 'Pa ',   &
           missing_value=missing_value)
            id_plfc_deep =  &
           register_diag_field (mod_name, 'plfc_deep', axes(1:2),   &
           Time, 'pressure at lfc ', 'Pa ',   &
           missing_value=missing_value)
            id_plzb_deep =  &
           register_diag_field (mod_name, 'plzb_deep', axes(1:2),   &
            Time, 'pressure at lzb ', 'Pa ',   &
           missing_value=missing_value)
            id_xcape_deep =  &
           register_diag_field (mod_name, 'xcape_deep', axes(1:2),  &
           Time, 'cape', 'J/kg',   &
           missing_value=missing_value)
            id_coin_deep =  &
           register_diag_field (mod_name, 'coin_deep', axes(1:2),   &
           Time, 'convective inhibition ', 'J/kg',   &
           missing_value=missing_value)
            id_dcape_deep =  &
           register_diag_field (mod_name, 'dcape_deep', axes(1:2), &
           Time, 'time tendency of cape ', 'J/kg/sec',   &
           missing_value=missing_value)
            id_qint_deep =  &
           register_diag_field (mod_name, 'qint_deep', axes(1:2),   &
           Time, 'column moisture ', 'kg(h2o)/m**2',   &
           missing_value=missing_value)
            id_a1_deep =  &
           register_diag_field (mod_name, 'a1_deep', axes(1:2),   &
           Time, 'fractional area of cu subensemble ', 'percent',   &
           missing_value=missing_value)
            id_amax_deep =  &
           register_diag_field (mod_name, 'amax_deep', axes(1:2),   &
           Time, 'fractional area of largest cu subensemble ',  &
           'percent',  missing_value=missing_value)
            id_amos_deep =  &
           register_diag_field (mod_name, 'amos_deep', axes(1:2),   &
           Time, 'uppr lmt on frac area from moisture', 'percent',   &
           missing_value=missing_value)
            id_tprea1_deep =  &
           register_diag_field (mod_name, 'tprea1_deep', axes(1:2), &
            Time, 'area wtd total precip ', 'mm/day',   &
           missing_value=missing_value)
            id_ampta1_deep =  &
           register_diag_field (mod_name, 'ampta1_deep', axes(1:2), &
           Time, 'meso cld frac', 'percent',   &
           missing_value=missing_value)
            id_omint_deep =  &
           register_diag_field (mod_name, 'omint_deep', axes(1:2), &
            Time, 'integrated low-lvl displ', 'Pa ',   &
           missing_value=missing_value)
            id_rcoa1_deep =  &
           register_diag_field (mod_name, 'rcoa1_deep', axes(1:2),  &
           Time, 'area wtd cnvctv precip ', 'mm/day',   &
           missing_value=missing_value)


end subroutine register_fields 

!####################################################################

subroutine read_restart_donner_deep

!---------------------------------------------------------------------
!   this subroutine reads a restart file previously written by this 
!   module.
!---------------------------------------------------------------------

      integer             :: next(2), dt(2), cal
      integer             :: next_conv, old_freq
      integer             :: unit, vers
      character(len=8)    :: chvers
!      type(time_type)   :: New_donner_deep_time
      integer                 :: secs_from_start, time_to_donner_deep, &
                                                 secs, days
      real, dimension(:,:,:), allocatable :: tempbl_old, ratpbl_old
      integer :: old_model_levels_in_sfcbl
      integer :: k
      integer  :: ncont_in, n, nn
      character(len=32) :: tracername_in
      logical, dimension(ncont) :: success

!-------------------------------------------------------------------- 
!   open the restart file.
!--------------------------------------------------------------------- 
      unit = open_restart_file ('INPUT/donner_deep.res', 'read')

!--------------------------------------------------------------------- 
!   read and check restart version number
!-------------------------------------------------------------------- 
      read (unit) vers
      if ( .not. any(vers == restart_versions) ) then
        write (chvers,'(i4)') vers
        call error_mesg ('read_restart_donner_deep', &
            'restart version '//chvers//' cannot be read '//&
            'by this module version', FATAL)
      endif

!--------------------------------------------------------------------
!   read time step and next occurrence from restart file. override 
!   provisional values, unless calendar type has changed.
!---------------------------------------------------------------------
      if (vers < 5) then
      read (unit) next, dt, cal
      else
      read (unit) next_conv, old_freq
      endif

!      Old_time_step = set_time (dt(1), dt(2))
      if (vers < 3) then
        if (cal == get_calendar_type() ) then
          conv_alarm     =  next(1)  
        else
  if (mpp_pe() == mpp_root_pe()) then
            call error_mesg ('read_restart_donner_deep',  &
        'current calendar not same as restart calendar, '//&
         'calling donner_deep on first timestep', NOTE)
! ADD 12-04: if calendar changes, do call now
            conv_alarm = 0
  endif
        endif
!     else
      else if (vers == 3 .or. vers == 4) then
!   next is the time remaining until the next timestep.
        conv_alarm     = next(1)
      else if (vers >= 5) then
        conv_alarm = next_conv
      endif

!--------------------------------------------------------------------
!  override previously calculated next time to call donner_deep_mod if 
!  timestep as input from namelist differs from that in restart file.
!--------------------------------------------------------------------
!         if (Donner_deep_timestep /= Old_time_step) then
        if (vers < 5) then
          if (donner_deep_freq     /= dt(1)        ) then
             conv_alarm     = conv_alarm     - dt(1) + donner_deep_freq
             if (conv_alarm     > 0         ) then
               if (mpp_pe() == mpp_root_pe()) then
                 call error_mesg ('donner_deep_init',   &
                     'donner_deep time step has changed, '//&
                      'next donner_deep time also changed', NOTE)
               end if
            endif
          endif
        else
          if (donner_deep_freq     /= old_freq     ) then
             conv_alarm     = conv_alarm  - old_freq + donner_deep_freq
             if (conv_alarm     > 0         ) then
               if (mpp_pe() == mpp_root_pe()) then
                 call error_mesg ('donner_deep_init',   &
                     'donner_deep time step has changed, '//&
                    'next donner_deep time also changed', NOTE)
               end if
            endif
          endif
          endif ! (vers < 5)

!---------------------------------------------------------------------
!   read the restart data fields.
!---------------------------------------------------------------------
      call read_data (unit, cemetf)
      call read_data (unit, cememf)
      call read_data (unit, xcape_lag)

!! RSH42902
      if ( vers >= 4) then
        if (vers >= 5) then
          read (unit) old_model_levels_in_sfcbl
          if (old_model_levels_in_sfcbl == model_levels_in_sfcbl) then
            call read_data (unit, tempbl)
            call read_data (unit, ratpbl)
        else
          allocate (tempbl_old(idf, jdf, old_model_levels_in_sfcbl))
           allocate (ratpbl_old(idf, jdf, old_model_levels_in_sfcbl))
           call read_data (unit, tempbl_old)
           call read_data (unit, ratpbl_old)
           if (old_model_levels_in_sfcbl > model_levels_in_sfcbl) then
             tempbl(:,:,1:model_levels_in_sfcbl) = &
                       tempbl_old(:,:,1:model_levels_in_sfcbl)
             ratpbl(:,:,1:model_levels_in_sfcbl) = &
                       ratpbl_old(:,:,1:model_levels_in_sfcbl)
           else
             tempbl(:,:,1:old_model_levels_in_sfcbl) = &
                    tempbl_old(:,:,1:old_model_levels_in_sfcbl)
             ratpbl(:,:,1:old_model_levels_in_sfcbl) = &
                    ratpbl_old(:,:,1:old_model_levels_in_sfcbl)
    tempbl(:,:,old_model_levels_in_sfcbl+1:model_levels_in_sfcbl) = 0.0
    ratpbl(:,:,old_model_levels_in_sfcbl+1:model_levels_in_sfcbl) = 0.0
             deallocate (ratpbl_old)
             deallocate (tempbl_old)
           endif 
          endif
          call read_data (unit, mass_flux)
          call read_data (unit, delta_ql )
          call read_data (unit, delta_qi )
          call read_data (unit, delta_qa )
        else ! (vers=5)
          allocate (tempbl_old(idf, jdf, nlev))
          allocate (ratpbl_old(idf, jdf, nlev))
          call read_data (unit, tempbl_old)
           call read_data (unit, ratpbl_old)
           do k=1, model_levels_in_sfcbl
             tempbl(:,:,k) =  tempbl_old(:,:,nlev-k+1)
             ratpbl(:,:,k) =  ratpbl_old(:,:,nlev-k+1)
           end do
           deallocate (ratpbl_old)
           deallocate (tempbl_old)
         endif
       else !(vers >=4)
!! supply initial values here if not on restart file
         tempbl = 0.0 
         ratpbl = 0.0
       endif
!! END RSH42902

       call read_data (unit, qint_lag )
       call read_data (unit, omint_acc)
      if (vers /= 1) then
        call read_data (unit, tprea1)
      else
        tprea1 = 0.0
      endif

!read donner cloud variables from restart file
      call read_data(unit, cell_cloud_frac  )
      call read_data(unit, cell_liquid_amt  )
      call read_data(unit, cell_liquid_size )
      call read_data(unit, cell_ice_amt     )
      call read_data(unit, cell_ice_size    )
 
      call read_data(unit, meso_cloud_frac  )
      call read_data(unit, meso_liquid_amt  )
      call read_data(unit, meso_liquid_size )
      call read_data(unit, meso_ice_amt     )
      call read_data(unit, meso_ice_size    )
      call read_data(unit, nsum             )

!------------------------------------------------------------------
!    if tracers are to be transported, see if tendencies are available
!    on the restart.
!------------------------------------------------------------------
      if (do_donner_tracer) then

!------------------------------------------------------------------
!    read the number of tracers whose tendencies are included in 
!    this file.
!-------------------------------------------------------------------
        if (vers >= 6) then
          success = .false.
          read (unit) ncont_in 

!--------------------------------------------------------------------
!    read each restart file tracer's name and see if it is to be 
!    transported in the current job.
!--------------------------------------------------------------------
          do n=1,ncont_in
            read (unit) tracername_in
            do nn=1,ncont

!--------------------------------------------------------------------
!    if the tracer is needed in the current job, read its data and
!    store it in the appropriate array. write a note indicating that 
!    the data has bben found and set a logical variable to also 
!    indicate such. exit this loop and process the next tracer present
!    in the restart file.
!--------------------------------------------------------------------
              if (trim(tracername_in) == trim(tracername(nn))) then
                call read_data(unit, tracer_tends(:,:,:,nn))
                if (mpp_pe() == mpp_root_pe() ) then
                  call error_mesg ('donner_deep_mod', &
                         'found tracer restart data for ' // &
                         trim(tracername(nn)), NOTE)
                endif
                success(nn) = .true.
                exit 

!---------------------------------------------------------------------
!    if the tracer in the restart file is not needed by the current
!    job, do a dummy read to get to the next record.
!---------------------------------------------------------------------
              else 
                if (nn == ncont) then
                  read (unit)
                endif
              endif
            end do
          end do

!---------------------------------------------------------------------
!    after having completely read the file, initialize the time ten-
!    dencies to 0.0 for any tracers whose tinme tendencies were not
!    found on the restart file and enter a message in the output file.
!---------------------------------------------------------------------
          do nn=1,ncont
            if (success(nn) ) then
            else
              call error_mesg ('donner_deep_mod', &
                  'did not find tracer restart data for ' //  &
                  trim(tracername(nn)) //  &
                  '; am initializing tendency to 0.0', NOTE)
              tracer_tends(:,:,:,nn) = 0.0
            endif   
          end do
        endif ! (vers >= 6)
      endif  ! (do_donner_tracer)

!-------------------------------------------------------------------- 
!   close the restart file.
!--------------------------------------------------------------------- 
      call close_file (unit)

!--------------------------------------------------------------------- 



end subroutine read_restart_donner_deep




!##################################################################

subroutine initialize_variables

!--------------------------------------------------------------------
!      this subroutine initializes the subdomain-dimensioned module 
!      variables when the donner_deep scheme is being coldstarted.
!--------------------------------------------------------------------
     
      cemetf       = 0.0
      cememf       = 0.0
      
      tracer_tends = 0.

      mass_flux = 0.
      delta_ql =0.
      delta_qi =0.
      delta_qa=0.

      xcape_lag     = 0.0
      qint_lag     = 0.0
      omint_acc     = 0.0
      tprea1    = 0.0
      tempbl       = 0.0
      ratpbl       = 0.0

!initialize cloud variables
      cell_cloud_frac  = 0.0
      cell_liquid_amt  = 0.0
      cell_liquid_size = 0.0
      cell_ice_amt     = 0.0
      cell_ice_size    = 0.0

      meso_cloud_frac  = 0.0
      meso_liquid_amt  = 0.0
      meso_liquid_size = 0.0
      meso_ice_amt     = 0.0
      meso_ice_size    = 0.0
      nsum             = 1

end subroutine initialize_variables


!#####################################################################

subroutine allocate_variables

!--------------------------------------------------------------------
!   this subroutine allocates space for the subdomain-dimensioned 
!   module variables.
!--------------------------------------------------------------------

     allocate ( cemetf         (idf, jdf, nlev ) )
     allocate ( cememf         (idf, jdf, nlev ) )

     allocate ( mass_flux      (idf, jdf, nlev ) )
     allocate ( delta_ql       (idf, jdf, nlev ) )
     allocate ( delta_qi       (idf, jdf, nlev ) )
     allocate ( delta_qa       (idf, jdf, nlev ) )

     allocate (tracer_tends    (idf, jdf, nlev, ncont) )

     allocate ( xcape_lag      (idf, jdf ) )
     allocate ( tempbl       (idf, jdf,model_levels_in_sfcbl) )
     allocate ( ratpbl       (idf, jdf,model_levels_in_sfcbl ) )
     allocate ( qint_lag       (idf, jdf ) )
     allocate ( omint_acc       (idf, jdf ) )
     allocate ( tprea1      (idf, jdf ) )

!variables which are passed to the radiation
     allocate (cell_cloud_frac (idf, jdf, nlev ) )
     allocate (cell_liquid_amt (idf, jdf, nlev ) )
     allocate (cell_liquid_size(idf, jdf, nlev ) )
     allocate (cell_ice_amt    (idf, jdf, nlev ) )
     allocate (cell_ice_size   (idf, jdf, nlev ) )

     allocate (meso_cloud_frac (idf, jdf, nlev ) )
     allocate (meso_liquid_amt (idf, jdf, nlev ) )
     allocate (meso_liquid_size(idf, jdf, nlev ) )
     allocate (meso_ice_amt    (idf, jdf, nlev ) )
     allocate (meso_ice_size   (idf, jdf, nlev ) )
     allocate (nsum            (idf, jdf ) )
!-------------------------------------------------------------------


end subroutine allocate_variables



!####################################################################




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      2. ROUTINES CALLED BY DONNER_DEEP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################

subroutine cape_calculation_driver (is, ie, js, je, temp,   &
                                    mixing_ratio, ttnd, qtnd, pfull, &
                                    Don_cape   )

integer,                intent(in)    ::  is, ie, js, je
real, dimension(:,:,:), intent(in)    ::  temp, mixing_ratio, ttnd,  &
                                          qtnd, pfull
type(donner_cape_type), intent(inout) ::  Don_cape

       real, dimension (size(temp,1),          &
                        size(temp,2), size(temp,3) )     ::         &
                              cape_input_temp, cape_input_vapor
       integer :: i,j,k

!----------------------------------------------------------------------
!    define the input temperature and moisture fields as those updated
!    by the current step tendencies. these should be very close to the
!    values that willbe seen on the next step.
!----------------------------------------------------------------------
        cape_input_temp = temp + ttnd
        cape_input_vapor = mixing_ratio + qtnd

!--------------------------------------------------------------------
!    call generate_cape_sounding to produce a high-resolution atmos-
!    pheric sounding to be used to evaluate cape.
!--------------------------------------------------------------------
         call generate_cape_sounding (cape_input_temp,    &
                                      cape_input_vapor, pfull, Don_cape)

!--------------------------------------------------------------------
!    call calculate_cape to evaluate cape in each model column.
!--------------------------------------------------------------------
         call calculate_cape (is, ie, js, je, Don_cape)

!---------------------------------------------------------------------
!   call integrate_vapor to produce column integral of water vapor.
!--------------------------------------------------------------------
         call integrate_vapor (Don_cape)

!--------------------------------------------------------------------
!   save the temperature and water vapor profiles in the boundary layer
!   that were used as input in this cape calculation and the values of 
!   cape and vertically-integrated water vapor that are output for use 
!   on the next timestep. the temperature and vapor profiles from the 
!   surface boundary layer (lowest model_layers_in_sfcbl levels) are 
!   used on the next step to reduce temporal noise in cape, as presc-
!   ribed in the catendb closure. the cape and water vapor output
!   from this calculation are needed to define the time tendency of 
!   these quantities, which are needed by the convection parameteriz-
!   ation.
!--------------------------------------------------------------------

       do k=1,model_levels_in_sfcbl
         do j=1,jsize    
           do i=1,isize       
             tempbl(i+is-1,j+js-1 ,k) = cape_input_temp  (i,j,nlev-k+1)
             ratpbl(i+is-1,j+js-1 ,k) = cape_input_vapor(i,j,nlev-k+1) 
           end do
         end do
       end do

       do j=1,jsize    
         do i=1,isize       
           qint_lag(i+is-1,j+js-1 ) = Don_cape%qint(i,j)
           xcape_lag(i+is-1,j+js-1 ) = Don_cape%xcape(i,j)
         end do
       end do
         
       call donner_column_cape_call (tempbl, ratpbl, ttnd, qtnd)

!---------------------------------------------------------------------
!    keep count of the number of points in the subdomain that have been 
!    initialized, so that initialization can be turned off once all 
!    points have been treated.
!---------------------------------------------------------------------
         pts_processed_cape = pts_processed_cape +   &
                               size(temp  ,1)*size(temp  ,2)
         if (pts_processed_cape == total_pts) then
           if (coldstart) coldstart = .false.
           pts_processed_cape = 0
           cape_alarm              = cape_alarm   + donner_deep_freq
         endif

end subroutine cape_calculation_driver 




!####################################################################

subroutine donner_convection_driver (is, ie, js, je, temp,    &
                                     mixing_ratio, phalf, pfull, dt,&
                                     omint, land, Don_cape, Don_conv, &
                                     qrat, ahuco, tracers,  &
                                     qlin, qiin, cf)

!---------------------------------------------------------------------
!   donner_convection_driver manages the calculation of the effects
!   of deep convection on atmospheric fields.
!---------------------------------------------------------------------

integer,                  intent(in)           :: is, ie, js, je
real, dimension(:,:,:),   intent(in)           :: temp, mixing_ratio, &
                                                  phalf, pfull 
real,                     intent(in)           :: dt
real, dimension(:,:),     intent(in)           :: omint, land
type(donner_cape_type),   intent(inout)        :: Don_cape
type(donner_conv_type),   intent(inout)        :: Don_conv
real, dimension(:,:,:),   intent(inout)        :: qrat, ahuco
real, dimension(:,:,:),   intent(in), optional :: qlin, qiin, cf
real, dimension(:,:,:,:), intent(in), optional :: tracers        

      integer :: i,j,k
       real, dimension (size(temp,1),          &
                        size(temp,2), size(temp,3) )     ::         &
                              dmeso_3d, xliq_3d

       real, dimension (size(temp,1), size(temp,2),  &
                                      size(temp,3)+1)  :: mhalf_3d
       real, dimension (size(temp,1),          &
                        size(temp,2), size(temp,3) )     ::         &
                              cape_input_temp, cape_input_vapor

!--------------------------------------------------------------------
!    define the temperature and mixing ratio fields to be used in de-
!    termining the cape sounding values. they are the model input values
!    except for the lowest two levels, where the values from the pre-
!    vious step are used (catendb closure) to prevent time tendencies
!    of cape being produced  that are a result of temporal noise con-
!    fined to the surface layer.
!--------------------------------------------------------------------
        do k=1, nlev-model_levels_in_sfcbl
          do j=1,jsize
            do i=1,isize
              cape_input_temp(i,j,k)  = temp  (i,j,k)
              cape_input_vapor(i,j,k) = mixing_ratio(i,j,k)
            end do
          end do
        end do

        do k=nlev-model_levels_in_sfcbl+1, nlev
          do j=1,jsize
            do i=1,isize
! if loop needed if starting from data where tempbl was not present
              if (tempbl(i+is-1,j+js-1,nlev-k+1) > 0.0) then  
               cape_input_temp(i,j,k) = tempbl(i+is-1, j+js-1,nlev-k+1)
               cape_input_vapor(i,j,k) = ratpbl(i+is-1, j+js-1,nlev-k+1)
              else
               cape_input_temp(i,j,k)  = temp  (i,j,k)
               cape_input_vapor(i,j,k) = mixing_ratio(i,j,k)
              endif
            end do
          end do
        end do

!---------------------------------------------------------------------
!    call generate_cape_sounding to produce a high-resolution atmos-
!    pheric sounding to be used to evaluate cape.
!---------------------------------------------------------------------
        call generate_cape_sounding (cape_input_temp, cape_input_vapor,&
                                     pfull, Don_cape )

!--------------------------------------------------------------------
!    call calculate_cape to evaluate cape in each model column.
!--------------------------------------------------------------------
        call calculate_cape (is, ie, js, je, Don_cape)              

!---------------------------------------------------------------------
!    calculate column-integrated water vapor in each column.
!--------------------------------------------------------------------
        call integrate_vapor (Don_cape)

!---------------------------------------------------------------------
!   call cupar to calculate normalized deep convective forcing
!---------------------------------------------------------------------
        if (present(tracers)) then
        call cupar (is, ie, js, je, dt, omint, phalf, Don_conv,  &
                    Don_cape, tracers)
        else
        call cupar (is, ie, js, je, dt, omint, phalf, Don_conv,  &
                    Don_cape)
        endif

!---------------------------------------------------------------------
!    call define_donner_anvil_ice to define the distribution of anvil
!    ice associated with deep convection.
!---------------------------------------------------------------------
        call define_donner_anvil_ice (is, js, pfull, temp, Don_conv)

!---------------------------------------------------------------------
!    call define_donner_mass_flux to define the distribution of mass
!    flux associated with deep convection.
!---------------------------------------------------------------------
        call define_donner_mass_flux (is,ie,js,je, pfull, phalf, xliq_3d,   &
                                      dmeso_3d, mhalf_3d,   &
                                      Don_conv)

!---------------------------------------------------------------------
!    call adjust_tiedtke_inputs to make needed adjustments based on
!    donner_deep_mod outputs.
!---------------------------------------------------------------------
      call adjust_tiedtke_inputs (is, ie, js, je, pfull, temp, &
                                  mixing_ratio, Don_conv, qrat, &
                                  ahuco) 
 
!---------------------------------------------------------------------
!    when strat_cloud is active, call strat_cloud_donner_tend to 
!    define increments to cloudice, cloudwater and cloud area associated
!    with deep convection. these are returned to moist_processes_mod,
!    where the pre-donner fields are incremented and then passed to 
!    strat_cloud_mod.
!    strat_cloud_donner_tend also removes vapor and temperature tendencies
!    corresponding to these increments from the Donner cumulus
!    thermal forcing and moisture forcing, which included
!    them as evaporatation and/or sublimation in mulsub.
!    assumptions used in strat_cloud_donner_tend to relate detrainment
!    to net mass fluxes differ from those in mulsub, so the 
!    increments here do not balance those in mulsub. the difference
!    remains as a phase change.
!---------------------------------------------------------------------
         if (present(qlin  )) then
          call strat_cloud_donner_tend (is, ie, js, je, dmeso_3d,   &
                                        xliq_3d, dt, Don_conv%xice, &
                                        mhalf_3d, phalf, qlin, qiin, cf)
         endif

!---------------------------------------------------------------------
!    call cell_liquid_size_comp to compute the cell liquid effective 
!    droplet diameter.
!---------------------------------------------------------------------
!       call cell_liquid_size_comp(pfull, temp, Don_conv, land = land)
        call cell_liquid_size_comp(pfull, temp, Don_conv,        land)

!---------------------------------------------------------------------
!    save variables needed by the radiation package.
!---------------------------------------------------------------------
        call donner_deep_sum (is, js, Don_conv)


!---------------------------------------------------------------------
!  determine if all points in this subdomain have 
!  been integrated on this timestep. if so, set the point counter to
!  zero and define the next time at which this module is to be executed.
!----------------------------------------------------------------------
           pts_processed_conv = pts_processed_conv + isize*jsize
           if (pts_processed_conv >= total_pts) then
             conv_alarm     = conv_alarm     + donner_deep_freq 
             pts_processed_conv = 0 
           endif



end subroutine donner_convection_driver

!####################################################################

subroutine adjust_tiedtke_inputs (is, ie, js, je, pfull, temp, &
                                    mixing_ratio, Don_conv, qrat, &
                                    ahuco) 

integer, intent(in) :: is, ie, js, je
real, dimension(:,:,:), intent(in) :: pfull, temp, mixing_ratio
type(donner_conv_type), intent(inout) :: Don_conv
real, dimension(:,:,:), intent(inout) :: qrat, ahuco
 
real ::   acell,pzm,rfun,qsat,esat
real, dimension(size(temp,3))  :: prinp, trf, qrf
integer :: i, j, k

!     calculate adjustments to u00 for Tiedtke stratiform cloud
!     parameterization. See "Tiedtke u00 adjustment" notes, 11/22/02.
!     Calculate ratio of large-scale specific humidity to humidity 
!     in environment of convective system.
!     

         do j=1, size(temp,2)
            do i=1, size(temp,1)
             pzm = Don_conv%pmd_v(i,j)-200.e02
             do k=1,size(temp,3)
               acell=0.
               prinp(k)=pfull(i,j,k)
               trf(k)=temp(i,j,k)
               qrf(k)=mixing_ratio(i,j,k)
               if (qrf(k) .lt. 0.) qrf(k)=0.
               if (Don_conv%cual(i,j,k) .ne. 0.) then
                  acell = Don_conv%cual(i,j,k)
               end if
               if ((prinp(k) .le. pzm) .and.   &
                   (prinp(k) .ge. Don_conv%pztm_v(i,j))) then
                 ahuco(i,j,k)=Don_conv%cual(i,j,k)
                 acell=Don_conv%cual(i,j,k) - Don_conv%ampta1(i,j)
               end if
               if (acell .lt. 0.) acell=0.
               rfun=0.
               if ((prinp(k) .ge. Don_conv%pztm_v(i,j)) .and. &
                   (prinp(k) .le. DOn_conv%pmd_v(i,j))) then
                 ahuco(i,j,k)=Don_Conv%cual(i,j,k) + &
                                    Don_conv%ampta1(i,j)
                 if (ahuco(i,j,k) .gt. 1.) ahuco(i,j,k)=1.
                 rfun=1.
               end if
               if ((prinp(k) .ge. Don_conv%pmd_v(i,j)) .and.  &
                   (prinp(k) .le. Don_conv%pb_v(i,j))) then
                 ahuco(i,j,k)=Don_conv%cual(i,j,k) + &
                                    Don_conv%ampta1(i,j)
                 if (ahuco(i,j,k) .gt. 1.) ahuco(i,j,k)=1.
                 rfun=1.-.3*(Don_conv%pmd_v(i,j)-prinp(k))/ &
                       (Don_conv%pmd_v(i,j)-Don_conv%pb_v(i,j))
               end if
               call lookup_es(trf(k),esat)
               qsat=rocp*esat/(prinp(k)+esat*(rocp-1.))
               if (prinp(k) .lt. Don_conv%pztm_v(i,j)) rfun=0.
               if (prinp(k) .le. Don_conv%pb_v(i,j)) then
                 qrat(i,j,k)=qrf(k)-acell*qsat-Don_conv%ampta1(i,j)* &
                 rfun*qsat
               if ((iequ(Don_conv%cual(i,j,k),1.) .eq. -10.) &
                    .and. (qrf(k) .ne. 0.)) then
                      if (qrat(i,j,k) .gt. 0. ) then
                            qrat(i,j,k)=  qrat(i,j,k)/      &
                            (1.-Don_conv%cual(i,j,k))
                           qrat(i,j,k)=qrf(k)/qrat(i,j,k)
                       else

!    qrat < 0 indicates insufficient large-scale moisture for
!    convective system if convective system is saturated or
!    downdraft has relative humidity r  set qrat to -10 as a
!    flag that this situation exists
!
                          qrat(i,j,k)=-10.
                      end if
                   else
                     qrat(i,j,k)=1.

                  end if
               end if
!cljd
!                if (Don_conv%ampta1(i,j) .ne. 0.) then
!                write(6,*) 'i,j,k,cual,ampt= ',i,j, &
!                k, Don_conv%cual(i,j,k),Don_conv%ampta1(i,j)
!                write (6,*) 'pb,prinp,trf, acell= ',Don_conv%pb_v(i,j),prinp(k), &
!                trf(k),acell
!                 write(6,*) 'rfun,qsat,qrat,qrf= ',rfun,qsat,qrat,qrf(k)
!                if (k .eq. size(temp,3)) stop 
!                end if
!cljd
                end do
             end do
         end do



end subroutine adjust_tiedtke_inputs 




!#####################################################################

subroutine define_output_fields (is, ie, js, je, dt, mixing_ratio, &
                                 ttnd, qtnd, precip,  &
                                 Don_conv, mtot, qltend, qitend,  &
                                 qatend, qtrceme)
integer, intent(in) :: is, ie, js, je
real, intent(in)  :: dt
real, dimension(:,:,:), intent(in) :: mixing_ratio
real, dimension(:,:,:), intent(out) :: ttnd, qtnd   
real, dimension(:,:), intent(out) :: precip         
real, dimension(:,:,:), intent(out), optional :: mtot, qltend,  &
                                                 qitend, qatend
real, dimension(:,:,:,:), intent(out),optional  :: qtrceme           
                                                 
type(donner_conv_type), intent(inout) :: Don_conv


     integer :: i, j, k
!----------------------------------------------------------------------
!    define the temperature and moisture tendencies from deep
!    convection. if the moisture tendency results in
!    the production of a negative value, reduce the tendency to produce
!    a zero value for moisture, and save the modified tendency for 
!    analysis purposes.
!----------------------------------------------------------------------
         do k=1,nlev
           do j=1,jsize      
             do i=1,isize       
               ttnd  (i,j,k) = cemetf(i+is-1,j+js-1 ,k)*dt
               qtnd  (i,j,k) = cememf(i+is-1,j+js-1 ,k)*dt
               if ((mixing_ratio(i,j,k) + qtnd(i,j,k)) .lt. 0.) then
                 if (mixing_ratio(i,j,k) .gt. 0.) then
                    qtnd  (i,j,k) = -mixing_ratio(i,j,k)/dt
                   Don_conv%cememf_mod(i,j,k) = -mixing_ratio(i,j,k)/dt
                 else 
                   Don_conv%cememf_mod(i,j,k) = 0.
                   qtnd(i,j,k) = 0.0
                 end if
               else
                 Don_conv%cememf_mod(i,j,k) = cememf(i+is-1,j+js-1 ,k)
               end if
             end do
           end do
         end do

         mtot(:,:,:) = mass_flux(is:ie, js:je,:)
         qltend(:,:,:) = delta_ql (is:ie, js:je,:)
         qitend(:,:,:) = delta_qi (is:ie, js:je,:)
         qatend(:,:,:) = delta_qa (is:ie, js:je,:)

!----------------------------------------------------------------------
!    define the rainfall and snowfall accrued on the current timestep
!    from deep convection. 
!    note: tprea1    [mm/day] * 1/86400 [day/sec] * 1/1000 [ m/mm] * 
!          1000 [kg(h2o)/m**3] * dt [sec] = kg/m2, as desired. 
!----------------------------------------------------------------------
      do j=1,jsize       
        do i=1,isize      
!         if (coldT(i,j)) then
!           snow(i,j) = tprea1(i,j+js-1 )*dt/86400.
!           rain(i,j) = 0.0
!         else
!           rain(i,j) = tprea1(i,j+js-1 )*dt/86400.
!           snow(i,j) = 0.0
!         endif
         precip (i,j) = tprea1(i+is-1,j+js-1 )*dt/86400.
        end do
      end do

!--------------------------------------------------------------------
!    if tracers are being transported by donner_deep_mod, retrieve their
!    time tendencies for return to the calling routine.
!--------------------------------------------------------------------
!       if (do_donner_tracer) then
        if (present (qtrceme)) then
          qtrceme(:,:,:,:) = tracer_tends(is:ie,js:je,:,:)
!       else
!         qtrceme(:,:,:,:) = 0.0                            
        endif

end subroutine define_output_fields 



!####################################################################


subroutine deallocate_local_variables (Don_conv, Don_cape)

type(donner_conv_type), intent(inout) :: Don_conv
type(donner_cape_type), intent(inout) :: Don_cape

      deallocate (Don_conv%ceefc            )
      deallocate (Don_conv%cecon            )
      deallocate (Don_conv%cemfc            )
      deallocate (Don_conv%cememf_mod      )
      deallocate (Don_conv%cual             )
      deallocate (Don_conv%fre              )
      deallocate (Don_conv%elt              )
      deallocate (Don_conv%cmus             )
      deallocate (Don_conv%ecds             )
      deallocate (Don_conv%eces          )
      deallocate (Don_conv%emds          )
      deallocate (Don_conv%emes          )
      deallocate (Don_conv%qmes             )
      deallocate (Don_conv%wmps          )
      deallocate (Don_conv%wmms          )
      deallocate (Don_conv%tmes          )
      deallocate (Don_conv%dmeml         )
      deallocate (Don_conv%uceml         )
      deallocate (Don_conv%umeml         )
      deallocate (Don_conv%xice          )
      deallocate (Don_conv%qtren1          )
      deallocate (Don_conv%qtceme          )
      deallocate (Don_conv%xgcm1          )
      deallocate (Don_conv%qtmes1          )
      deallocate (Don_conv%wtp1          )
      deallocate (Don_conv%dgeice        )
      deallocate (Don_conv%cuqi          )
      deallocate (Don_conv%cuql          )
      deallocate (Don_conv%cell_liquid_eff_diam  )
      deallocate (Don_conv%cell_ice_geneff_diam  )
      deallocate (Don_conv%dcape            )
      deallocate (Don_conv%a1            )
      deallocate (Don_conv%amax             )
      deallocate (Don_conv%amos          )
      deallocate (Don_conv%ampta1        )
      deallocate (Don_conv%rcoa1           )
      deallocate (Don_conv%contot         )
      deallocate (Don_conv%emdi_v           )
      deallocate (Don_conv%pb_v             )
      deallocate (Don_conv%pmd_v            )
      deallocate (Don_conv%pztm_v )

      deallocate (Don_cape%coin   )
      deallocate (Don_cape%plcl   )
      deallocate (Don_cape%plfc    )
      deallocate (Don_cape%plzb    )
      deallocate (Don_cape%xcape  )
      deallocate (Don_cape%parcel_r)
      deallocate (Don_cape%parcel_t   )
      deallocate (Don_cape%cape_p   )
      deallocate (Don_cape%env_r   )
      deallocate (Don_cape%env_t  )
      deallocate (Don_cape%model_p  )
      deallocate (Don_cape%model_r )
      deallocate (Don_cape%model_t )
      deallocate (Don_cape%qint    )   

end subroutine deallocate_local_variables 





!###################################################################

subroutine donner_deep_netcdf (is, ie,js, je,Time, Don_conv, Don_cape, &
                               pmass)

type(time_type), intent(in) :: Time
integer, intent(in) :: is, ie, js, je
type(donner_conv_type), intent(in) :: Don_conv
type(donner_cape_type), intent(in) :: Don_cape
real, dimension(:,:,:), intent(in)   :: pmass


     logical :: used
     integer :: k, n
     real, dimension(size(Don_conv%ceefc,1),   &
                     size(Don_conv%ceefc,2))  ::   tempdiag


      used = send_data (id_cemetf_deep, cemetf(is:ie,js:je,:), Time,   &
                        is, js, 1)
      used = send_data (id_ceefc_deep,   Don_conv%ceefc(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_cecon_deep,  Don_conv%cecon(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_cemfc_deep, Don_conv%cemfc(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_cememf_deep, cememf(is:ie,js:je,:),  &
                         Time, is, js, 1)
      used = send_data (id_cememf_mod_deep, Don_conv%cememf_mod(:,:,:),&
                        Time, is, js, 1)
       used = send_data (id_cual_deep, Don_conv%cual(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_fre_deep, Don_conv%fre   (:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_elt_deep, Don_conv%elt   (:,:,:), Time,&
                         is, js, 1)
       used = send_data (id_cmus_deep, Don_conv%cmus   (:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_ecds_deep, Don_conv%ecds(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_eces_deep, Don_conv%eces(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_emds_deep, Don_conv%emds(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_emes_deep, Don_conv%emes(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_qmes_deep, Don_conv%qmes(:,:,:), Time,&
                         is, js, 1)
       used = send_data (id_wmps_deep, Don_conv%wmps(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_wmms_deep, Don_conv%wmms(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_tmes_deep, Don_conv%tmes(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_dmeml_deep, Don_conv%dmeml(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_uceml_deep, Don_conv%uceml(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_umeml_deep, Don_conv%umeml(:,:,:), Time,&
                        is, js, 1)
       used = send_data (id_xice_deep, Don_conv%xice(:,:,:), Time,&
                         is, js, 1)
    do n=1,ncont  
       used = send_data (id_qtren1(n), Don_conv%qtren1(:,:,:,n), Time,     &
                         is, js, 1)
       used = send_data (id_qtmes1(n), Don_conv%qtmes1(:,:,:,n), Time,   &
                         is, js, 1)
       used = send_data (id_wtp1(n), Don_conv%wtp1(:,:,:,n), Time,&
                         is, js, 1)
       used = send_data (id_qtceme(n), Don_conv%qtceme(:,:,:,n), Time,     &
                         is, js, 1)


       tempdiag(:,:) = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtren1(:,:,k,n) * &
                         pmass(:,:,k)
       end do
       used = send_data (id_qtren1_col(n), tempdiag, Time, is, js)
       tempdiag(:,:) = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtmes1(:,:,k,n) * &
                         pmass(:,:,k)
       end do
       used = send_data (id_qtmes1_col(n), tempdiag, Time, is, js)
       tempdiag(:,:) = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%wtp1(:,:,k,n) *   &
                         pmass(:,:,k)
       end do
       used = send_data (id_wtp1_col(n), tempdiag, Time, is, js)
       tempdiag(:,:) = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtceme(:,:,k,n) * &
                         pmass(:,:,k)
       end do
       used = send_data (id_qtceme_col(n), tempdiag, Time, is, js)
    end do

       used = send_data (id_dgeice_deep, Don_conv%dgeice(:,:,:), Time,&
                         is, js, 1)
       used = send_data (id_cuqi_deep, Don_conv%cuqi(:,:,:), Time,&
                         is, js, 1)
       used = send_data (id_cuql_deep, Don_conv%cuql(:,:,:), Time,&
                         is, js, 1)
       used = send_data (id_dgeliq_deep,  &
                         Don_conv%cell_liquid_eff_diam(:,:,:), Time,&
                         is, js, 1)
             

        used = send_data (id_plcl_deep, Don_cape%plcl(:,:), Time,&
                          is, js)
        used = send_data (id_plfc_deep, Don_cape%plfc(:,:), Time,&
                        is, js)
        used = send_data (id_plzb_deep, Don_cape%plzb(:,:), Time,&
                        is, js)
        used = send_data (id_xcape_deep,    &
                          xcape_lag(is:ie,js:je),  &
                          Time, is, js)
        used = send_data (id_coin_deep, Don_cape%coin(:,:), Time,&
                        is, js)
         used = send_data (id_dcape_deep, Don_conv%dcape(:,:), Time,&
                        is, js)
         used = send_data (id_qint_deep,     &
                           qint_lag(is:ie,js:je),   &
                           Time, is, js)
         used = send_data (id_a1_deep, Don_conv%a1(:,:), Time,&
                           is, js)
         used = send_data (id_amax_deep, Don_conv%amax(:,:), Time,&
                           is, js)
         used = send_data (id_amos_deep, Don_conv%amos(:,:), Time,&
                           is, js)
         used = send_data (id_tprea1_deep,    &
                           tprea1(is:ie,js:je),  &
                           Time, is, js)
         used = send_data (id_ampta1_deep, Don_conv%ampta1(:,:), Time,&
                            is, js)
         used = send_data (id_omint_deep,    &
                           omint_acc(is:ie,js:je),  &
                           Time, is, js)
         used = send_data (id_rcoa1_deep, Don_conv%rcoa1(:,:), Time,&
                            is, js)

end subroutine donner_deep_netcdf



!####################################################################


subroutine initialize_local_variables (Don_conv, Don_cape)

type(donner_conv_type), intent(inout) :: Don_conv
type(donner_cape_type), intent(inout) :: Don_cape


!--------------------------------------------------------------------
!   this subroutine initializes the arrays of the donner_conv_type 
!   variable Don_conv.
!--------------------------------------------------------------------

      allocate (Don_conv%ceefc           (isize, jsize, nlev) )
      allocate (Don_conv%cecon           (isize, jsize, nlev) )
      allocate (Don_conv%cemfc           (isize, jsize, nlev) )
      allocate (Don_conv%cememf_mod      (isize, jsize, nlev) )
      allocate (Don_conv%cual            (isize, jsize, nlev) )
      allocate (Don_conv%fre             (isize, jsize, nlev) )
      allocate (Don_conv%elt             (isize, jsize, nlev) )
      allocate (Don_conv%cmus            (isize, jsize, nlev) )
      allocate (Don_conv%ecds         (isize, jsize, nlev) )
      allocate (Don_conv%eces         (isize, jsize, nlev) )
      allocate (Don_conv%emds         (isize, jsize, nlev) )
      allocate (Don_conv%emes         (isize, jsize, nlev) )
      allocate (Don_conv%qmes         (isize, jsize, nlev) )
      allocate (Don_conv%wmps         (isize, jsize, nlev) )
      allocate (Don_conv%wmms         (isize, jsize, nlev) )
      allocate (Don_conv%tmes         (isize, jsize, nlev) )
      allocate (Don_conv%dmeml        (isize, jsize, nlev) )
      allocate (Don_conv%uceml        (isize, jsize, nlev) )
      allocate (Don_conv%umeml        (isize, jsize, nlev) )
      allocate (Don_conv%xice         (isize, jsize, nlev) )
      allocate (Don_conv%qtren1       (isize, jsize, nlev,ncont) )
      allocate (Don_conv%qtceme       (isize, jsize, nlev,ncont) )
      allocate (Don_conv%xgcm1        (isize, jsize, nlev,ncont) )
      allocate (Don_conv%qtmes1       (isize, jsize, nlev,ncont) )
      allocate (Don_conv%wtp1         (isize, jsize, nlev,ncont) )
      allocate (Don_conv%dgeice       (isize, jsize, nlev) )
      allocate (Don_conv%cuqi         (isize, jsize, nlev) )
      allocate (Don_conv%cuql         (isize, jsize, nlev) )
      allocate (Don_conv%cell_liquid_eff_diam (isize, jsize, nlev) )
      allocate (Don_conv%cell_ice_geneff_diam (isize, jsize, nlev) )
      allocate (Don_conv%dcape           (isize, jsize) )
      allocate (Don_conv%a1           (isize, jsize) )
      allocate (Don_conv%amax         (isize, jsize) )
      allocate (Don_conv%amos         (isize, jsize) )
      allocate (Don_conv%ampta1       (isize, jsize) )
      allocate (Don_conv%rcoa1        (isize, jsize) )
      allocate (Don_conv%contot        (isize, jsize) )
      allocate (Don_conv%emdi_v          (isize, jsize) )
      allocate (Don_conv%pb_v            (isize, jsize) )
      allocate (Don_conv%pmd_v           (isize, jsize) )
      allocate (Don_conv%pztm_v          (isize, jsize) )


      Don_conv%cell_liquid_eff_diam = 0.0
      Don_conv%cell_ice_geneff_diam = 0.0


     Don_conv%ceefc      = 0.0
     Don_conv%cecon      = 0.0
     Don_conv%cemfc      = 0.0
     Don_conv%cememf_mod = 0.0
     Don_conv%cual       = 0.0
     Don_conv%fre        = 0.0
     Don_conv%elt        = 0.0
     Don_conv%cmus       = 0.0
     Don_conv%ecds    = 0.0
     Don_conv%eces    = 0.0
     Don_conv%emds    = 0.0
     Don_conv%emes    = 0.0
     Don_conv%qmes    = 0.0
     Don_conv%wmps    = 0.0
     Don_conv%wmms    = 0.0
     Don_conv%tmes    = 0.0
     Don_conv%dmeml   = 0.0
     Don_conv%uceml   = 0.0
     Don_conv%umeml   = 0.0
     Don_conv%xice    = 0.0
     Don_conv%qtren1    = 0.0
     Don_conv%qtceme    = 0.0
     Don_conv%xgcm1    = 0.0
     Don_conv%qtmes1    = 0.0
     Don_conv%wtp1    = 0.0
     Don_conv%dgeice  = 0.0
     Don_conv%cuqi    = 0.0
     Don_conv%cuql    = 0.0
     Don_conv%dcape      = 0.0
     Don_conv%a1      = 0.0
     Don_conv%amax    = 0.0
     Don_conv%amos    = 0.0
     Don_conv%ampta1  = 0.0
     Don_conv%rcoa1   = 0.0
     Don_conv%contot   = 0.0
     Don_conv%emdi_v     = 0.0
     Don_conv%pb_v       = 0.0
     Don_conv%pmd_v      = 0.0
     Don_conv%pztm_v     = 0.0

      allocate (Don_cape%coin       (isize, jsize) )
      allocate (Don_cape%plcl       (isize, jsize) )
      allocate (Don_cape%plfc       (isize, jsize) )
      allocate (Don_cape%plzb       (isize, jsize) )
      allocate (Don_cape%xcape      (isize, jsize) )
      allocate (Don_cape%qint       (isize, jsize) )
      allocate (Don_cape%parcel_r   (isize, jsize, ncap) )
      allocate (Don_cape%parcel_t      (isize, jsize, ncap) )
      allocate (Don_cape%cape_p      (isize, jsize, ncap) )
      allocate (Don_cape%env_r      (isize, jsize, ncap) )
      allocate (Don_cape%env_t      (isize, jsize, ncap) )
      allocate (Don_cape%model_p    (isize, jsize, nlev) )
      allocate (Don_cape%model_r    (isize, jsize, nlev) )
      allocate (Don_cape%model_t    (isize, jsize, nlev) )

                Don_cape%coin     =0. 
                Don_cape%plcl      = 0.0
                Don_cape%plfc      = 0.0
                Don_cape%plzb       = 0.0
                Don_cape%xcape   =  0.0 
                Don_cape%qint      = 0.0
                Don_cape%parcel_r   = 0.0
                Don_cape%parcel_t  = 0.0
                Don_cape%cape_p    = 0.0
                Don_cape%env_r     = 0.0
                Don_cape%env_t     = 0.0
                Don_cape%model_p  = 0.0
                Don_cape%model_r    = 0.0
                Don_cape%model_t   = 0.0

!-------------------------------------------------------------------


end subroutine initialize_local_variables



!####################################################################



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      3. LOWER LEVEL ROUTINES            
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!#####################################################################

subroutine define_donner_mass_flux (is, ie, js,je,  pfull, phalf,  &
                                    xliq_3d, dmeso_3d, mhalf_3d,&
                                    Don_conv)

integer, intent(in)  :: is, ie, js, je
real, dimension(:,:,:), intent(in) :: phalf, pfull                   
real, dimension(:,:,:), intent(inout) :: xliq_3d, dmeso_3d, mhalf_3d
type(donner_conv_type), intent(inout) :: Don_conv

            integer  :: ilon, jlat, k
    real :: dnna, dnnb
!      dnna       weight for averaging full-level mass fluxes to half levels
!      dnnb       weight for averaging full-level mass fluxes to half levels
!
!        Code to link Donner convection to Tiedtke/Rotstayn cloud
!        fraction/microphysics
!
!        Calculate dmeso_3d, the detrainment rate from the convective system
!        dmeso_3d has units of s**-1 and is (1/rho) dM/dz, where M is the
!        detraining mass flux of the convective system (kg/(m**2 s))
!
         do ilon=1,isize       
           do jlat=1,jsize      
             do k=1,nlev
               xliq_3d(ilon,jlat,k)=0.
               if (Don_conv%xice(ilon,jlat,k) .ge. 1.0e-10) then
                 dmeso_3d(ilon,jlat,k)=Don_conv%emes(ilon,jlat,k)/ &
                                       Don_conv%xice(ilon,jlat,k)
               else
                 dmeso_3d(ilon,jlat,k)=0.
               end if
             end do
             do k=1,nlev-1
               if ((Don_conv%uceml(ilon,jlat,k) .le. 1.0e-10)  .and.   &
                   (Don_conv%umeml(ilon,jlat,k) .le. 1.0e-10) .and. &
                   (Don_conv%dmeml(ilon,jlat,k) .le. 1.0e-10) ) then
                 mhalf_3d(ilon,jlat,k)=0.
                 mhalf_3d(ilon,jlat,k+1)=0.
                 cycle
               end if
               dnna=phalf(ilon,jlat,k+1)-pfull(ilon,jlat,k)
               dnnb=pfull(ilon,jlat,k+1)-phalf(ilon,jlat,k+1)
               mhalf_3d(ilon,jlat,k+1)=  &
                                 (dnnb*Don_conv%uceml(ilon,jlat,k) + &
                                  dnna*Don_conv%uceml(ilon,jlat,k+1)   &
                                 +dnnb*Don_conv%umeml(ilon,jlat,k)    &
                                 +dnna*Don_conv%umeml(ilon,jlat,k+1)  &
                                 +dnnb*Don_conv%dmeml(ilon,jlat,k)    &
                                 +dnna*Don_conv%dmeml(ilon,jlat,k+1))/ &
                                  (dnna+dnnb)
             end do
             mhalf_3d(ilon,jlat,nlev+1)=0.
             mhalf_3d(ilon,jlat,1)=0.
             do k=1,nlev
               if ((Don_conv%uceml(ilon,jlat,k) .le. 1.0e-10) .and.   &
                   (Don_conv%umeml(ilon,jlat,k) .le. 1.0e-10) .and.   &
                   (Don_conv%dmeml(ilon,jlat,k) .le. 1.0e-10) ) then
!                mtot(ilon,jlat,k)=0.
                 mass_flux(ilon+is-1,jlat+js-1,k)=0.
                 cycle
               end if
!              mtot(ilon,jlat,k) = Don_conv%uceml_3d(ilon,jlat,k) + &
               mass_flux(ilon+is-1,jlat+js-1,k) =   &
                                      Don_conv%uceml(ilon,jlat,k) + & 
                                      Don_conv%umeml(ilon,jlat,k) + & 
                                      Don_conv%dmeml(ilon,jlat,k) 
             end do 
           end do 
         end do 

end subroutine define_donner_mass_flux 



!######################################################################

subroutine define_donner_anvil_ice (is, js, pfull, temp, Don_conv)

integer, intent(in) :: is, js
real, dimension(:,:,:), intent(in) :: pfull, temp
type(donner_conv_type), intent(inout) :: Don_conv
            

!      ampu       fractional area of mesoscale circulation
!      emdi       vertical integral of  mesoscale-downdraft 
!                 sublimation (mm/d)
!      prinp      pressures at full levels (index 1 at model top) (Pa)
!      rmuf       mesoscale updraft mass flux  (kg/[(m**2) s])
!      trf        temperature at full GCM levels (K) 
!                 index 1 at model top
!      tprei      total convective-system precipitation (mm/d)
!      xice       anvil ice (kg(ice)/kg)
!      przm       pressure at anvil base (Pa), inferred from presence 
!                 of ice at full GCM levels
!      prztm      pressure at anvil top (Pa), inferred from presence 
!                 of ice at full GCM levels
!      dgeicer    generalized effective size of ice crystals in anvil 
!                 (micrometers) defined as in Fu (1996, J. Clim.)
  integer :: ilon, jlat, k, jgl, n, igl
  real :: ampu, emdi,         tprei, przm, prztm, prent, dgeicer
  real, dimension(size(pfull,3)) :: prinp, rmuf, trf, xice

    if (in_diagnostics_window) then
      do n=1,ncols_in_window
        if (Don_conv%ampta1 (i_dc(n),j_dc(n)) /= 0.0) then
          write (unit_dc(n), '(a, 4e20.12)') 'pre prean: ampta1, &
                            &emdi, contot, tprei', &
                            Don_conv%ampta1 (i_dc(n),j_dc(n)), &
                            Don_conv%emdi_v    (i_dc(n),j_dc(n)), &
                            Don_conv%contot  (i_dc(n),j_dc(n)), &
                            tprea1 (igl_dc(n), jgl_dc(n))
        endif
      end do
    endif

    do ilon=1,isize        
      do jlat=1,jsize       
        jgl=jlat+js-1
        igl=ilon+is-1
        ampu=Don_conv%ampta1(ilon,jlat)
        emdi=Don_conv%emdi_v(ilon,jlat)
!        contot=Don_conv%contot_v(ilon,jlat)
        tprei=tprea1(igl ,jgl)
        do k=1,nlev
          prinp(k)=pfull(ilon,jlat,k)
          rmuf(k)=Don_conv%umeml(ilon,jlat,k)
          trf(k)=temp  (ilon,jlat,k)
        end do
        przm=prinp(nlev)
        prztm=-10.

        call prean(ampu,Don_conv%contot(ilon,jlat),emdi,prinp,  &
                   rmuf,trf,tprei,xice)

        do k=1,nlev
          Don_conv%xice(ilon,jlat,k)=xice(k)
          if ((xice(k) .ge. 1.e-10) .and. (prztm .le. 0.)) then 
            prztm=prinp(k)
            cycle
          end if
        end do
        do k=1,nlev
          if ((xice(k) .lt. 1.e-11) .and. (prinp(k) .ge. prztm)) then
            przm=prinp(k-1)
            exit
          end if
        end do

        do k=1,nlev
          prent=prinp(k)
          dgeicer=0.
          if (xice(k) .ge. 1.e-10) call andge(prent,przm,prztm,&
                                              dgeicer)
          Don_conv%dgeice(ilon,jlat,k)=dgeicer
        end do
      end do
    end do

    if (in_diagnostics_window) then
      do n=1,ncols_in_window
        do k=1,nlev
          if (Don_conv%xice   (i_dc(n),j_dc(n),k) > 0.0) then
             write (unit_dc(n), '(a, i4, e10.3, e20.12)')  &
                               'post prean: pressure, xice', &
                                k, pfull(i_dc(n), j_dc(n),k),    &
                                Don_conv%xice(i_dc(n), j_dc(n), k)   
           endif
         end do
       end do
     endif

end subroutine define_donner_anvil_ice


!###################################################################

subroutine integrate_vapor (Cape)

type(donner_cape_type), intent(inout) :: Cape

     integer :: i, j, k, n
     real  :: sum

         do i=1,isize      
           do j=1,jsize      
             sum    =     Cape%env_r(i,j,1)*(    Cape%cape_p(i,j,1) -  &
                          Cape%cape_p(i,j,2))
             do k=2,ncap-1
               sum   = sum            +     Cape%env_r(i,j,k)*  &
                        0.5*(    Cape%cape_p(i,j,k-1) -  &
                                 Cape%cape_p(i,j,k+1))
             end do
             sum     = sum            +     Cape%env_r(i,j,ncap)*  &
                      (    Cape%cape_p(i,j,ncap-1) -  &
                           Cape%cape_p(i,j,ncap))
             Cape%qint(i,j) = sum           /gravm
           end do
         end do

!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
         if (in_diagnostics_window) then
           do n=1,ncols_in_window
            write (unit_dc(n), '(a, f20.10)')  &
                 'integrate_vapor: qint= ',Cape%qint(i_dc(n),j_dc(n))  
           end do
         endif
end subroutine integrate_vapor 



!######################################################################
!###################################################################

subroutine prean(ampu,contot,emdi,prf,rmuf,trf,tprei,xice)

! Calculates ice content for mesoscale circulation.
!
!      Leo Donner
!      GFDL
!      17 May 2001
!
!   Input arguments
!
real,  intent(in)  ::  ampu  ! fractional area of mesoscale
                             !   circulation
real,  intent(in)  ::  contot ! ratio of convective to total
                              !   precipitation
real,  intent(in)  ::  emdi  ! vertical integral of
                             !    mesoscale-downdraft sublimation  
                             !  (mm/d)
real, dimension(:), intent(in)  ::  prf ! pressure at full levels
                                        !   (Pa) index 1 at top of model
real, dimension(:), intent(in)  ::  rmuf ! mesoscale updraft mass flux
                                         !  (kg/[(m**2) s]) Index 1 
                                         !  at model top
real, dimension(:), intent(in)  ::  trf  ! temperature at full GCM 
                                         ! levels
                                         !   (K)  index 1 at model top
real,  intent(in)  ::  tprei ! total convective-system
                             !   precipitation (mm/d)
!
!   Input arguments as module variables:
!
!   nlev                       number of GCM levels
!   rair  (module variable)    gas constant for dry air (J/(kg K))
!
!   Output arguments
!
real, dimension(:), intent(out) :: xice ! anvil ice (kg(ice)/kg)
! index 1 at model top
!
!      local workspace
!
  integer :: k                    ! vertical index
  integer :: kou                  ! counter
  real    :: emdiw
 real    :: rho              ! anvil (ht. avg.) air density (kg/(m**3))
 real    :: rhol
  real    :: tprew
 real    :: xicet                ! anvil ice work variable (kg(ice)/kg)
!
!   Initialize
!
 xicet=0.
 do k=1,nlev
    xice(k)=0.
   end do
!
!    Return if no mesoscale anvil exists.
!
  if (ampu .le. 0.) return
   if (contot .ge. 1.) return
 if ((tprei .le. 0.) .and. (emdi .le. 0.)) return
!
!    Convert precipitation and sublimation integral to kg/[(m**2) s].
!
   tprew=tprei/86400.
   emdiw=emdi/86400.
!
!     Calculate average air density over vertical extent of anvil.
!
   rho=0.
   kou=0
   do k=1,nlev
    if (rmuf(k) .ne. 0.) then
         kou=kou+1
         rhol=prf(k)/(rair*trf(k))
         rho=rho+rhol
        cycle

 end if
 end do
 if (kou .eq. 0) return
   rho=rho/kou
!
!      Calculate mesoscale ice content by balancing fall at anvil base with
!      mesoscale precipitation and sublimation in mesoscale downdraft.
!
  xicet=(1.-contot)*tprew
   xicet=xicet+emdiw
   xicet=xicet/(3.29*ampu)
   xicet=xicet**.862
  xicet=xicet/rho
!
!      Assign anvil ice to layers with postive meso updraft mass flux
!
  do k=1,nlev
    if (rmuf(k) .le. 0.) xice(k)=0.
      if (rmuf(k) .gt. 0.) xice(k)=xicet
  end do
! ljd
!   do k=1,nlev
! if (xice(k) .gt. .01) then
! write(6,*) 'prean ampu,contot,emdi= ',ampu,contot,emdi
! write(6,*) 'rho,rhol,kou= ',rho,rhol,kou
! write(6,*) 'tprei,rair= ',tprei, rair
!    do kk=1,nlev
!    write(6,*) 'k,prf,xice= ',kk,prf(kk),xice(kk)
!    write(6,*) 'k,prf,trf= ',kk,prf(kk),trf(kk)
!    write(6,*) 'k,prf,rmuf= ',kk,prf(kk),rmuf(kk)
!    end do
!    stop
!    end if
! end do
! ljd


end subroutine prean

!###################################################################
subroutine andge(p,pzm,pztm,dgeicer)
!
!  test routine for calculating anvil d_ge
!
!----------------------------------------------------------------------
       implicit none
!      Input Arguments
!
      integer, parameter          :: nmclev=6    ! number of levels
                                                 ! in anvil
      real, intent(in)            :: p           ! pressure
      real, intent(in)            :: pzm         ! anvil base pressure
     real, intent(in)            :: pztm         ! anvil top pressure
     real, intent(out)           :: dgeicer      ! generalized effective
                                                 !  size at pressure p 
                                                 !  (micrometers)
                                                 ! p, pzm, and pztm must
                                                 !  have identical units
!
!      Local Variables
!
     real, dimension(nmclev)     :: dgeice  ! generalized effective
!  size of hexagonal ice crystals, defined as in Fu (1996, J. Clim.)
!  values from Table 2 of McFarquhar et al. (1999, JGR) are averaged over
!  all grid boxes for which D_ge is defined for all altitudes between
!  9.9 and 13.2 km. index 1 at bottom of anvil
        real, dimension(nmclev)     :: relht   ! distance from anvil
!  base, normalized by total anvil thickness. from Table 2 of McFarquhar
!  et al. (1999, JGR) for grid boxes with data between 9.9 and 13.2 km.
!  index 1 at anvil bottom
!
        integer                     :: k               ! vertical index
       real                        :: znor            ! normalized distance

!  from anvil base
      data dgeice/38.5,30.72,28.28,25.62,24.8,13.3/
      data relht/0.,.3,.45,.64,.76,1./
!      write(6,*) 'p,pzm,pztm= ',p,pzm,pztm
        if (pzm .lt. pztm) then
          write(6,*) 'error in angde pzm,pztm= ',pzm,pztm
          stop
        end if
       if (pzm .eq. pztm) then
         znor=.5
         go to 12
       end if
       znor=(pzm-p)/(pzm-pztm)
12   continue
!      write(6,*) 'znor= ',znor
        do k=2,nmclev
!       write(6,*) 'relhts= ',relht(k-1),relht(k)
          if ((znor .ge. relht(k-1)) .and. (znor .le. relht(k))) then
             dgeicer=dgeice(k-1)+( (znor-relht(k-1))*(dgeice(k)- &
               dgeice(k-1))/(relht(k)-relht(k-1)) )
!       write(6,*) 'dgeicer= ',dgeicer
        go to 11
       end if
      end do
11   continue
!      do k=1,nmclev
!      write(6,*) 'k,dgeice,relht= ',k,dgeice(k),relht(k)
!      end do

 end subroutine andge

!###################################################################

subroutine cupar (is, ie, js, je, dt, omint, phalf, Don_conv,   &
                  Don_cape, tracers)

!---------------------------------------------------------------------
integer,                intent(in)    :: is, ie, js, je
real,                   intent(in)    :: dt
real, dimension(:,:),   intent(in)    :: omint
real, dimension(:,:,:), intent(in)    :: phalf
type(donner_conv_type), intent(inout) :: Don_conv
type(donner_cape_type), intent(inout) :: Don_cape
real, dimension(:,:,:,:), intent(in), optional :: tracers
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      Cupar drives the parameterization for deep cumulus convection
!      based on Donner (1993, J. Atmos. Sci.).
!
!-----------------------------------------------------------------------
!     On Input:
!
!        cape     convective available potential energy (J/kg)
!        cin      convective inhibtion (J/kg)
!        cpd      specific heat of dry air at constant pressure (J/(kg K))
!        cpv      specific heat of water vapor [J/(kg K)]
!        dcape    local rate of CAPE change by all processes
!                 other than deep convection [J/(kg s)]
!        dqls     local rate of change in column-integrated vapor
!                 by all processes other than deep convection
!                 {kg(H2O)/[(m**2) s]}
!        epsilo   ratio of molecular weights of water vapor to dry air
!        gravm    gravity constant [m/(s**2)]
!        ilon     longitude index
!        jlat     latitude index
!        mcu      frequency (in time steps) of deep cumulus
!        omint    integrated low-level displacement (Pa)
!        cape_p   pressure at Cape.F resolution (Pa)
!                 Index 1 at bottom of model.
!        plfc     pressure at level of free convection (Pa)
!        plzb     pressure at level of zero buoyancy (Pa)
!        pr       pressure at Skyhi vertical resolution (Pa)
!                 Index 1 nearest ground  
!        q        large-scale vapor mixing ratio at Skyhi vertical resolution
!                 [kg(h2O)/kg]
!                 Index 1 nearest ground 
!        qlsd     column-integrated vapor divided by timestep for cumulus
!                 parameterization {kg(H2O)/[(m**2) s]}
!        r        large-scale vapor mixing ratio at Cape.F resolution
!                 [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rpc      parcel vapor mixing ratio from Cape.F [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rd       gas constant for dry air (J/(kg K))
!        rlat     latent heat of vaporization (J/kg)
!        rv       gas constant for water vapor (J/(kg K))
!        t        large-scale temperature at Skyhi vertical resolution (K)
!                 Index 1 nearest ground
!        tcape    large-scale temperature at Cape.F resolution (K)
!                 Index 1 at bottom of model.
!        tpc      parcel temperature from from Cape.F (K)
!                 Index 1 at bottom of model.
!
!     On Input as Parameters:
!
!        kmax     number of vertical levels at Skyhi resolution
!        kpar     number of cumulus sub-ensembles
!        ncap     number of vertical levels in Cape.F resolution
!

      real, dimension (size(Don_Cape%model_t,1),   &
                       size(Don_cape%model_t,2)) :: &
                                                 qlsd_v, dqls_v
      real, dimension(kpar) :: arat
      real, dimension (size(Don_Cape%model_t,1),   &
                       size(Don_cape%model_t,2), kpar    ) :: arat_v

      real, dimension(size(Don_cape%model_t,3)) :: cmus, ecds,  &
                                                   eces, emds, pr, q, &
                                     t, cual, emes, disa, disb, disc, &
                                     disd, dise, elt, fre, qmes,   &
                                      tmes, wmms, wmps, sig, dmeml, &
                                      uceml, umeml
      real, dimension (ncap) ::qli0, qli1, qr, qt, r, ri, rl, tc,   &
                               cape_p, tcape, tpc, rpc
      real, dimension (size(Don_cape%model_t,1),   &
                       size(Don_cape%model_t,2),   &
                       size(Don_cape%model_t,3) ) ::   &
                                sig_v, cmus_v, cual_v, ecds_v, &
                                 eces_v, emds_v, emes_v, disa_v, &
                                  disb_v, disc_v, disd_v, dise_v, &
                                   dmeml_v, elt_v, fre_v, qmes_v, &
                                 tmes_v, uceml_v, umeml_v, wmms_v, &
                                  wmps_v, qtest, a1_vk, cuq_v,cuql_v

      real, dimension (size(Don_cape%model_t,1),   &
                       size(Don_cape%model_t,2),   &
                       size(Don_cape%model_t,3),ncont ) ::   &
                                       xgcm_v ,qtmes_v, qtren_v, wtp_v

      real, dimension (size(Don_cape%model_t,1),  &
                       size(Don_cape%model_t,2), ncont  ) ::   &
                                        xba_v

      real, dimension (size(Don_cape%model_t,1),   &
                       size(Don_cape%model_t,2), ncap)  ::     &
                                   qli0_v, qli1_v, qr_v, qt_v, ri_v, &
                                   rl_v

      real, dimension (size(Don_cape%model_t,1),   &
                       size(Don_cape%model_t,2) ) ::     &
                                   ampt_v, sfcqf_v, sfcsf_v, amax_v, &
                                   cmui_v,  tpre_v, emei_v, &
                                    dcape_v, &
                                   amos_v, rco_v, tfint, a1_v, pdeet1, &
                                   pdeet2
     logical, dimension (size(Don_cape%model_t,1),  &
                         size(Don_cape%model_t,2)) :: exit_flag,    &
                                                       nine_flag
       
      logical, dimension (size(Don_cape%model_t,1),  &
                          size(Don_cape%model_t,2),  &
                          size(Don_cape%model_t,3) ) ::   &
                                                              exit_3d

      integer :: unit
      integer :: n, k
      integer  :: i, j, jlat, ilon, jgl, jinv
      integer  :: igl
      real  :: disbar

!--------------------------------------------------------------------
!
!     test set to true for graphics
!     must also dimension sig for graphics
!
      logical test
      test=.false.
!
!     arat(i) is the ratio at cloud base of the fractional area
!     of ensemble i to ensemble 1.
!
!     GATE:
!
      data arat/1.,.26,.35,.32,.3,.54,.66/

           do i=1,isize        
             do j=1,jsize     
               dqls_v(i,j) = (Don_cape%qint(i,j) -   &
                              qint_lag(i+is-1,j+js-1))/dt
               qlsd_v(i,j) = Don_cape%qint(i,j)/donner_deep_freq
               Don_conv%dcape(i,j) = (Don_cape%xcape(i,j) -    &
               xcape_lag(i+is-1,j+js-1 ))/dt
             end do
           end do

       dcape_v(:,:) = Don_conv%dcape(:,:)

     if (do_donner_tracer .and. (.not. present (tracers)) ) then
       call error_mesg ('donner_deep_mod', &
    ' are doing donner_tracers but tracer array not passed to cupar', &
                                                           FATAL) 
     endif

     if (do_donner_tracer) then
!-------------------------------------------------------------------
!    define the tracer values at lowest model level.
!-------------------------------------------------------------------
         xba_v(:,:,:) = tracers(:,:,nlev,:)

!--------------------------------------------------------------------
!    save the tracer profiles as elements of the don_conv_type 
!    variable. 
!--------------------------------------------------------------------
         Don_conv%xgcm1(:,:,:,:)  = tracers(:,:,:,:)

!-------------------------------------------------------------------
!   define inverted tracer profile (index 1 nearest ground) for use in
!   cloud routines.
!------------------------------------------------------------------
       do k=1,nlev
         xgcm_v(:,:,k,:)=tracers(:,:,nlev-k+1,:)
       end do

    else
        xba_v = 0.
        xgcm_v = 0.
        Don_conv%xgcm1 = 0.
    endif ! (do_donner_tracers)


       if (in_diagnostics_window) then
          do n=1,ncols_in_window
            write (unit_dc(n), '(a, 2i4)')  &
           'cupar: entering cupar with i_dc, j_dc:',  i_dc(n), j_dc(n)
           end do
        endif

      do jlat=1,je-js+1  
        do ilon=1,ie-is+1  

          nine_flag(ilon,jlat) = .false.

!---------------------------------------------------------------------
!  check that criteria for deep convection are satisfifed.
!---------------------------------------------------------------------
          pdeet1(ilon,jlat) = Don_cape%plfc(ilon,jlat) -  &
                              Don_cape%plzb(ilon,jlat)
          pdeet2(ilon,jlat) = Don_cape%plfc(ilon,jlat) -  &
                              Don_cape%model_p(ilon,jlat,1)
          if ((Don_Cape%xcape (ilon,jlat) <= 0.) .or.  &
               (dcape_v(ilon,jlat) <= 0.) .or. & 
               (pdeet1(ilon,jlat) < pdeep_cv )       .or. &
               (pdeet2(ilon,jlat) < omint(ilon,jlat)) .or.  &
              (Don_cape%coin(ilon,jlat) > cdeep_cv) ) then
            exit_flag(ilon,jlat) = .true.
          else
            exit_flag(ilon,jlat) = .false.
          endif
        end do
      end do

    if (in_diagnostics_window) then
      do n=1,ncols_in_window
        write (unit_dc(n), '(a, f20.12, f20.12, e20.12, l4)')   &
            'in cupar: omint,cape,dcape, exit_flag',    &
           omint(i_dc(n),j_dc(n)), Don_Cape%xcape(i_dc(n),j_dc(n)),   &
           dcape_v(i_dc(n),j_dc(n)), exit_flag(i_dc(n),j_dc(n))
       end do
     endif

      do jlat=1,je-js+1  
        do ilon=1,ie-is+1
          jgl = js+jlat-1
          igl = is+ilon-1 



!-------------------------------------------------------------------
!   initialize fields.
!--------------------------------------------------------------------
          do j=1,nlev
            cemetf(igl ,jgl ,j          )=0.
            cememf(igl ,jgl ,j          )=0.
            cuq_v(ilon,jlat,j) = 0.0
            cuql_v(ilon,jlat,j) = 0.0
            tprea1(igl ,jgl ) = 0.0
          end do
        end do
      end do

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          write (unit_dc(n), '(a, 2f12.4)')   &
                   'in cupar: cpi,cpv= ',cpi,cpv
          write (unit_dc(n), '(a, 2f12.6, f12.2)')  &
                   'in cupar: rocp,rair,latvap= ',rocp,rair, latvap
          write (unit_dc(n), '(a, f12.7)') 'in cupar: rvap= ',rvap
          write (unit_dc(n), '(a, 2f14.7, f19.10)')  &
                    'in cupar: cape,cin,plzb= ',  &
               Don_cape%xcape(i_dc(n),j_dc(n)), &
               Don_cape%coin(i_dc(n),j_dc(n)), &
               Don_cape%plzb(i_dc(n),j_dc(n))
          write (unit_dc(n), '(a, f19.10)') 'in cupar: plfc= ', &
               Don_cape%plfc(i_dc(n),j_dc(n))
          write (unit_dc(n), '(a, 2e20.12)') 'in cupar: dqls,qlsd= ', &
               dqls_v(i_dc(n),j_dc(n)),&
               qlsd_v(i_dc(n),j_dc(n))
           do k=1,nlev-kstart_diag+1
                write (unit_dc(n), '(a, i4, f19.10, f20.14, e20.12)') &
                                   'in cupar: k,pr,t,q= ',k,   &
                    Don_cape%model_p(i_dc(n),j_dc(n),k),   &
                    Don_cape%model_t(i_dc(n),j_dc(n),k),   &
                    Don_cape%model_r(i_dc(n),j_dc(n),k)
            end do
          end do
        endif

      
!---------------------------------------------------------------------
!  call routine to calculate normalized cumulus forcings.
!---------------------------------------------------------------------
      sfcqf_v(:,:      )=0.
      sfcsf_v(:,:      )=0.
       do k=1,kpar
         arat_v(:,:,k) = arat(k)
       end do



     call mulsub_vect (ampt_v,arat_v,Don_cape%plzb,   &
                       Don_cape%model_p, phalf, Don_cape%model_r,  &
                       sfcqf_v, sfcsf_v,Don_cape%model_t, &
                       xba_v, xgcm_v,   &
                       amax_v,cmui_v,cmus_v, cual_v,cuq_v, cuql_v, &
                       ecds_v,eces_v,emds_v,emei_v,emes_v,disa_v, &
                       disb_v, disc_v,disd_v,dise_v,dmeml_v, &
                       elt_v,fre_v,qmes_v, qtmes_v, qtren_v,  &
                       tmes_v,tpre_v,uceml_v,umeml_v,wmms_v, &
                       wmps_v,  wtp_v, exit_flag, Don_conv        )


!ljd
!         do jlat=jminp,jmaxp
!        do ilon=iminp,imaxp
!         if (Don_conv%pmd_v(ilon,jlat) .ge. 200.e02) then
!         write(6,*) 'ilon,jlat,pb,pmd,pztm= ',ilon,jlat,  &
!                 Don_conv%pb_v(ilon,jlat), &
!        Don_conv%pmd_v(ilon,jlat), Don_conv%pztm_v(ilon,jlat)
!        stop
!        end if
!           if (ampt_v(ilon,jlat) .ne. 0.) then
!           write(6,*) 'cupar 1 ilon,jlat,jgl= ',ilon,jlat,jgl
!             write(6,*) 'cupar 1 contot= ',Don_conv%contot_v(ilon,jlat)
!           end if
!        end do
!        end do
!ljd


        if (in_diagnostics_window) then
          do n=1,ncols_in_window
            if (.not. exit_flag(i_dc(n),j_dc(n))) then
              write  (unit_dc(n), '(a, 2e20.12)')   &
                     'in cupar:  ampt,tpre= ',  &
                   ampt_v(i_dc(n),j_dc(n)),  tpre_v(i_dc(n),j_dc(n))
              do k=1,nlev-kstart_diag+1    
                write (unit_dc(n), '(a, i4, e20.12)')  &
                      'in cupar: k,cual= ',k,cual_v(i_dc(n),j_dc(n),k)
              end do
            endif
          end do
        endif


      do jlat=1,je-js+1   
        do ilon=1,ie-is+1       
          if (.not. exit_flag(ilon,jlat)) then
            amos_v(ilon,jlat)=0.
            rco_v(ilon,jlat)=tpre_v(ilon,jlat)* &
                             Don_conv%contot(ilon,jlat)
            if (tpre_v(ilon,jlat) .eq. 0.) then
              a1_v(ilon,jlat)=0.
              amax_v(ilon,jlat)=0.
              nine_flag(ilon,jlat) = .true.
            end if
          endif
        end do
      end do




     if (in_diagnostics_window) then
       do n=1,ncols_in_window
         if (exit_flag(i_dc(n),j_dc(n)) ) then
           write (unit_dc(n), '(a)')   &
                'in cupar: deep convection not present in&
                   &test column. exit_flag is set to .true.'
          else
            write (unit_dc(n), '(a)')   &
                   'in cupar: deep convection is possible&
                     & -- exit_flag is .false.'
           endif
         end do
      endif

!---------------------------------------------------------------------
!  interpolate forcings to resolution of cuclo.
!--------------------------------------------------------------------
      call polat_vect (dise_v, Don_cape%model_p, Don_cape%cape_p, qr_v) 
      call polat_vect (disa_v, Don_cape%model_p, Don_cape%cape_p, qt_v) 

!--------------------------------------------------------------------
!   initialize variables.
!--------------------------------------------------------------------
      do k=1,ncap
        do j=1,je-js+1   
          do i=1,ie-is+1      
            qli0_v(i,j,k)=0.
            qli1_v(i,j,k)=0.
            rl_v(i,j,k)=0.
            ri_v(i,j,k)=0.
          end do
        end do
      end do

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          if (.not. exit_flag(i_dc(n),j_dc(n)) ) then
            do k=1,ncap
              write (unit_dc(n), '(a, i4, e20.12, f20.14)')  &
                           'in cupar: k,qr,qt= ',k,   &
                            qr_v(i_dc(n),j_dc(n),k), &
                            qt_v(i_dc(n),j_dc(n),k)
            end do
            do k=1,nlev
              write (unit_dc(n), '(a, i4, 2e20.12)')  &
                            'in cupar: k,dise,disa= ',k,   &
                            dise_v(i_dc(n),j_dc(n),k),   &
                            disa_v(i_dc(n),j_dc(n),k)
            end do
          endif
        end do
      endif

!--------------------------------------------------------------------
!   call routine to close deep-cumulus parameterization.
!--------------------------------------------------------------------
      call cuclo_vect(dcape_v, Don_cape%cape_p,Don_cape%plfc, &
                      Don_cape%plzb, qli0_v,qli1_v,qr_v, &
                      qt_v,Don_cape%env_r, ri_v,rl_v,   &
                      Don_cape%parcel_r, Don_Cape%env_t, &
                      Don_cape%parcel_t,a1_v, exit_flag, nine_flag, &
                      is, ie, js, je)


!--------------------------------------------------------------------
!   vertical integral of normalized moisture forcing:
!-------------------------------------------------------------------
      do j=1,je-js+1  
        do i=1,ie-is+1      
          tfint(i,j) = 0.0
        end do
      end do


      do k=2,nlev
        do j=1,je-js+1  
          do i=1,ie-is+1      
            if (.not. exit_flag(i,j) ) then
              disbar = 0.5*(dise_v(i,j,k-1)+dise_v(i,j,k))
              tfint(i,j) = tfint(i,j) - disbar*   &
                           (Don_cape%model_p(i,j,k-1) - &
                            Don_cape%model_p(i,j,k))
             endif
          end do
        end do
      end do

      do k=1,nlev
        exit_3d(:,:,k) = exit_flag(:,:)
      end do


      do jlat=1,je-js+1  
        do ilon=1,ie-is+1      
          if (.not. exit_flag(ilon,jlat) ) then

!---------------------------------------------------------------------
!       restrict fractional area of first ensemble. see 
!       "a Bounds 6/7/97" notes.
!---------------------------------------------------------------------
            a1_v(ilon,jlat) = min(amax_v(ilon,jlat), a1_v(ilon,jlat))

!---------------------------------------------------------------------
!      fractional area further restricted by moisture constraint.
!      See "Moisture Constraint," 8/8/97.
!---------------------------------------------------------------------
            if ((tfint(ilon,jlat) .eq. 0.) .or.    &
                 (tpre_v(ilon,jlat) .eq. 0.)) then
              a1_v(ilon,jlat) = 0.
              nine_flag(ilon,jlat) = .true.
            end if
          else
            a1_v(ilon,jlat) = 0.
          endif
        end do
      end do

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          if (.not. exit_flag(i_dc(n),j_dc(n)) .and.   &
              .not. nine_flag(i_dc(n),j_dc(n)))  then
             write (unit_dc(n), '(a, e20.12)')  &
                        'in cupar: tfint= ',tfint(i_dc(n),j_dc(n))
             write (unit_dc(n), '(a, e20.12)')  &
                        'in cupar: a1_v = ',a1_v (i_dc(n),j_dc(n))
          endif
        end do
      endif

      do jlat=1,je-js+1  
        do ilon=1,ie-is+1      
          if (.not. exit_flag(ilon,jlat) .and.   &
              .not. nine_flag(ilon,jlat))  then
            tfint(ilon,jlat)=tfint(ilon,jlat)/gravm
            amos_v(ilon,jlat) = (dqls_v(ilon,jlat) +     &
                                 qlsd_v(ilon,jlat))/tfint(ilon,jlat)
            if (a1_v(ilon,jlat) .gt. amos_v(ilon,jlat))     &
                a1_v(ilon,jlat) = max(amos_v(ilon,jlat),0.)

          endif
        end do
      end do

!-----------------------------------------------------------------------
!    constrain area so that mean specific humidity and specific
!    humidity are realizeable
 

    if (in_diagnostics_window) then
      do n=1,ncols_in_window
        if (.not. exit_flag(i_dc(n),j_dc(n)) .and.   &
            .not. nine_flag(i_dc(n),j_dc(n)))  then
          write (unit_dc(n), '(a, 3e20.12, f12.5)')  &
                          'in cupar: tfint,amos,a1,gravm= ',  &
                     tfint(i_dc(n),j_dc(n)), amos_v(i_dc(n),j_dc(n)), &
                     a1_v(i_dc(n),j_dc(n)), gravm 
        endif
      end do
    endif

    do k=1,nlev
      a1_vk(:,:,k) = a1_v(:,:)
    end do

!---------------------------------------------------------------------
!   if a1 allows negative mixing ratio, reset its value.
!--------------------------------------------------------------------
    do k=1,nlev
      do j=1,je-js+1  
        do i=1,ie-is+1      
          if (.not. exit_flag(i,j) .and.   &
              .not. nine_flag(i,j))  then
            qtest(i,j,k) = Don_cape%model_r(i,j,k) + a1_vk(i,j,k)* &
                           donner_deep_freq*dise_v(i,j,k)
            if (qtest(i,j,k) .lt. 0.) then
!             a1_vk(i,j,k) = -q_v(i,j,k)/   &
              a1_vk(i,j,k) = -Don_cape%model_r(i,j,k)/   &
                             (dise_v(i,j,k)*donner_deep_freq)
            endif
          endif
        end do
      end do
    end do

      do k=2,nlev
        do j=1,je-js+1   
          do i=1,ie-is+1      
            a1_v(i,j) = min(a1_vk(i,j,k), a1_v(i,j))
          end do
        end do
      end do

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          if (.not. exit_flag(i_dc(n),j_dc(n)) ) then
            write (unit_dc(n), '(a, e20.12)')   &
                   'in cupar: a1= ',a1_v (i_dc(n),j_dc(n))
          endif
        end do
      endif

!-----------------------------------------------------------------------
!      remove normalization from cumulus forcings.
!-----------------------------------------------------------------------
      do j=1,nlev
        do jlat=1,je-js+1  
          do ilon=1,ie-is+1  
            if (.not. exit_3d  (ilon,jlat,j) ) then


              disa_v(ilon,jlat,j)=disa_v(ilon,jlat,j)*a1_v(ilon,jlat)
              disb_v(ilon,jlat,j)=disb_v(ilon,jlat,j)*a1_v(ilon,jlat)
              disc_v(ilon,jlat,j)=disc_v(ilon,jlat,j)*a1_v(ilon,jlat)
              disd_v(ilon,jlat,j)=disd_v(ilon,jlat,j)*a1_v(ilon,jlat)
              dise_v(ilon,jlat,j)=dise_v(ilon,jlat,j)*a1_v(ilon,jlat)
              fre_v(ilon,jlat,j)=fre_v(ilon,jlat,j)*a1_v(ilon,jlat)
              elt_v(ilon,jlat,j)=elt_v(ilon,jlat,j)*a1_v(ilon,jlat)
              cmus_v(ilon,jlat,j)=cmus_v(ilon,jlat,j)*a1_v(ilon,jlat)
              ecds_v(ilon,jlat,j)=ecds_v(ilon,jlat,j)*a1_v(ilon,jlat)
              eces_v(ilon,jlat,j)=eces_v(ilon,jlat,j)*a1_v(ilon,jlat)
              emds_v(ilon,jlat,j)=emds_v(ilon,jlat,j)*a1_v(ilon,jlat)
              emes_v(ilon,jlat,j)=emes_v(ilon,jlat,j)*a1_v(ilon,jlat)
              wmms_v(ilon,jlat,j)=wmms_v(ilon,jlat,j)*a1_v(ilon,jlat)
              wmps_v(ilon,jlat,j)=wmps_v(ilon,jlat,j)*a1_v(ilon,jlat)
              tmes_v(ilon,jlat,j)=tmes_v(ilon,jlat,j)*a1_v(ilon,jlat)
              qmes_v(ilon,jlat,j)=qmes_v(ilon,jlat,j)*a1_v(ilon,jlat)
              cual_v(ilon,jlat,j)=cual_v(ilon,jlat,j)*a1_v(ilon,jlat)
              uceml_v(ilon,jlat,j)=uceml_v(ilon,jlat,j)*a1_v(ilon,jlat)
              umeml_v(ilon,jlat,j)=umeml_v(ilon,jlat,j)*a1_v(ilon,jlat)
              dmeml_v(ilon,jlat,j)=dmeml_v(ilon,jlat,j)*a1_v(ilon,jlat)
!! RSH ADDS:
              qtmes_v(ilon,jlat,j,:)=qtmes_v(ilon,jlat,j,:)*a1_v(ilon,jlat)
              qtren_v(ilon,jlat,j,:)=qtren_v(ilon,jlat,j,:)*a1_v(ilon,jlat)
              wtp_v(ilon,jlat,j,:) = wtp_v(ilon,jlat,j,:)*a1_v(ilon,jlat)
            endif
          end do
        end do
      end do


!do k=nlev-4,nlev
!           do i=1,ie-is+1  
!              if (in_diagnostics_window) then
!      do n=1,ncols_in_window
!                 if (.not. exit_flag(i,j_dc(n)) ) then
!                   if (disa_v(i,j_dc(n),k) .ne. 0.) then
!              write (unit_dc(n), '(a, 3i4, e20.12)') 'in cupar:k, jtest,i, disa= ',&
! k, jgl_dc(n), i,  disa_v(i,j_dc(n),k)
!             call error_mesg ('cupar', &
!      'disa found to be /= 0.0 in top 5 levels', FATAL)
!                   endif
!                   if (fre_v(i,j_dc(n),k) .lt. 0.) then
!              write (unit_dc(n), '(a, 3i4, 2e20.12)') 'in cupar:  k,jtest,i,fre,a1=',&
!       k,jgl_dc(n),i, fre_v(i,j_dc(n),k),  &
!a1_v(i,j_dc(n))
!     call error_mesg ('cupar',  &
!                   ' fre_v found to be < 0.0', FATAL)
!            end if
!                   if (cmus_v(i,j_dc(n),k) .lt. 0.) then
!              write (unit_dc(n), '(a, 3i4, e20.12)') 'in cupar: k,jtest,i,cmus= ', &
!k,jgl_dc(n),i, cmus_v(i,j_dc(n),k)
!     call error_mesg ('cupar',  &
!                 ' cmus_v found to be < 0.0', FATAL)
!                   end if    
!               endif
!     end do
!             endif
!           end do
!       end do


        do jlat=1,je-js+1  
          jgl = jlat+js-1 
          do ilon=1,ie-is+1  
            if (.not. exit_flag(ilon,jlat)) then
              igl = ilon+is-1
              do j=1,nlev
                jinv=nlev+1-j

                cemetf(igl ,jgl ,jinv          )=disa_v(ilon,jlat,j)
                Don_conv%ceefc(ilon,jlat,jinv )=disb_v(ilon,jlat,j)
                Don_conv%cecon(ilon,jlat,jinv   )=disc_v(ilon,jlat,j)
                Don_conv%cemfc(ilon,jlat,jinv  )=disd_v(ilon,jlat,j)
                cememf(igl ,jgl ,jinv   )=dise_v(ilon,jlat,j)
                Don_conv%cual(ilon,jlat,jinv  )=cual_v(ilon,jlat,j)
                Don_conv%fre   (ilon,jlat,jinv )=fre_v(ilon,jlat,j)
                Don_conv%elt   (ilon,jlat,jinv) = elt_v(ilon,jlat,j)
                Don_conv%cmus   (ilon,jlat,jinv )=cmus_v(ilon,jlat,j)
                Don_conv%ecds(ilon,jlat,jinv )=ecds_v(ilon,jlat,j)
                Don_conv%eces(ilon,jlat,jinv )=eces_v(ilon,jlat,j)
                Don_conv%emds(ilon,jlat,jinv )=emds_v(ilon,jlat,j)
                Don_conv%emes(ilon,jlat,jinv  )=emes_v(ilon,jlat,j)
                Don_conv%qmes(ilon,jlat,jinv )=qmes_v(ilon,jlat,j)
                Don_conv%wmps(ilon,jlat,jinv  )=wmps_v(ilon,jlat,j)
                Don_conv%wmms(ilon,jlat,jinv )=wmms_v(ilon,jlat,j)
                Don_conv%tmes(ilon,jlat,jinv )=tmes_v(ilon,jlat,j)
                Don_conv%dmeml(ilon,jlat,jinv)=dmeml_v(ilon,jlat,j)
                Don_conv%uceml(ilon,jlat,jinv )=uceml_v(ilon,jlat,j)
                Don_conv%umeml(ilon,jlat,jinv  )=umeml_v(ilon,jlat,j)
                Don_conv%cuqi(ilon,jlat,jinv )  = cuq_v(ilon,jlat,j)
                Don_conv%cuql(ilon,jlat,jinv )  = cuql_v(ilon,jlat,j)
                Don_conv%qtren1(ilon,jlat,jinv,:) = qtren_v(ilon,jlat,j,:)
                Don_conv%qtmes1(ilon,jlat,jinv,:) = qtmes_v(ilon,jlat,j,:)
                Don_conv%qtceme(ilon,jlat,jinv,: )  =   &
                    qtmes_v(ilon,jlat,j,:) +    qtren_v(ilon,jlat,j,:)
                Don_conv%wtp1(ilon,jlat,jinv,: )  = wtp_v(ilon,jlat,j,:)
              end do
            endif
          end do
        end do

!      unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' donner_deep,in cupar'
!        write (unit,*) ' Don_conv%cuql'
!        write (unit,*) Don_conv%cuql
!      call close_file (unit)
!     if (mpp_pe() == 2) then
!print *, 'cemetf in cupar', (cemetf(33,3,k),k=1,nlev)
!     endif



      do jlat=1,je-js+1     
        do ilon=1,ie-is+1      
        jgl=jlat+js-1
        igl=ilon+is-1
          if (.not. exit_flag(ilon,jlat) ) then
            Don_conv%a1(ilon,jlat) = a1_v(ilon,jlat)
            Don_conv%amax(ilon,jlat) = amax_v(ilon,jlat)
            Don_conv%amos(ilon,jlat         )=amos_v(ilon,jlat)
            tprea1(igl ,jgl      )=tpre_v(ilon,jlat)*a1_v(ilon,jlat)
            Don_conv%ampta1(ilon,jlat)=ampt_v(ilon,jlat)*a1_v(ilon,jlat)
            Don_conv%rcoa1(ilon,jlat)=rco_v(ilon,jlat)*a1_v(ilon,jlat)
            Don_conv%emdi_v(ilon,jlat)=  &
                              Don_conv%emdi_v(ilon,jlat)*a1_v(ilon,jlat)
          endif
        end do
      end do

!ljd
!         do jlat=jminp,jmaxp
!        do ilon=iminp,imaxp
!           if (ampt_v(ilon,jlat) .ne. 0.) then
!           write(6,*) 'cupar 3 ilon,jlat,jgl= ',ilon,jlat,jgl
!             write(6,*) 'cupar 3 ampt,tpre= ',Don_conv%ampta1(ilon,jlat), &
!              tprea1_3d(ilon,jgl)
!             do k=1,nlev
!             write(6,*) 'k,cual= ',k,Don_conv%cual(ilon,jlat,k)
!             end do
!           end if
!        end do
!        end do
!ljd
!ljd
!      do jlat=jminp,jmaxp
!      do ilon=iminp,imaxp
!              if (Don_conv%ampta1(ilon,jlat) .ne. 0.) then
!             write(6,*) ' Cupar 4 i,j,jgl= ',ilon,jlat,jgl
!             write(6,*) 'emdi= ',Don_conv%emdi_v(ilon,jlat)
!             do k=1,nlev
!               write(6,*) 'k,umem= ',k,Don_conv%umeml(ilon,jlat,k)
!             end do
!             end if
!             end do
!             end do
!ljd

      if ( in_diagnostics_window) then
        do n=1,ncols_in_window
          do i=1,ie-is+1     
            if (.not. exit_flag(i,j_dc(n)) ) then
              if (Don_conv%rcoa1(i,j_dc(n) ) .gt. 1.) then
                 call error_mesg ('cupar',  &
             ' the following columns contain deep convection', NOTE)
                write (unit_dc(n), '(a, 2i4, e20.12)')  &
                              'in cupar: i,jtest,meso a= ',  &
                        i, jgl_dc(n),Don_conv%rcoa1(i,j_dc(n))
                 write (unit_dc(n), '(a, 2e20.12)') &
                    'in cupar: amax,ampt= ',Don_conv%amax(i,j_dc(n)), &
                     Don_conv%ampta1(i,j_dc(n))
              end if
            endif
          end do
        end do
      endif
!  end do

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          if (.not. exit_flag(i_dc(n),j_dc(n)) ) then
            if (tprea1(igl_dc(n),jgl_dc(n) ) .ne. 0.) then
              do k=1,nlev
                if ((Don_cape%model_p(i_dc(n),j_dc(n),k) .gt. 100.e02) .and.   &
                    (Don_cape%model_p(i_dc(n),j_dc(n),k) .lt. 500.e02)) then    
                  if (disa_v(i_dc(n),j_dc(n),k) .ne. 0.) then
                    write (unit_dc(n), '(a, 3i4, f20.14)')    &
                                'in cupar: j_dc,i_dc,k,t= ',  &
                              j_dc(n), i_dc(n),k,  &
                                 Don_cape%model_t(i_dc(n),j_dc(n),k)
                    write (unit_dc(n), '(a, e20.12, i4, 2e20.12)')&
                           'in cupar: tprea1,k,pr,cemetf= ',  &
                            tprea1(igl_dc(n),jgl_dc(n)), k,    &
                            Don_cape%model_p(i_dc(n),j_dc(n),k),   &
                            cemetf(igl_dc(n),jgl_dc(n) ,k )
                  endif
                endif
              end do
            endif
          endif
        end do
      endif

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          write (unit_dc(n), '(a, 2e20.12)') 'in cupar: contot,tpre=', &
              Don_conv%contot(i_dc(n),j_dc(n)), tpre_v(i_dc(n),j_dc(n))
          write (unit_dc(n), '(a, 2e20.12)') 'in cupar: a1,ampt =',  &
                         Don_conv%a1 (i_dc(n),j_dc(n)), &
                         Don_conv%ampta1(i_dc(n),j_dc(n))
        end do
      endif

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          write (unit_dc(n), '(a, e20.12)')  'in cupar: amax= ', &
                         Don_conv%amax(i_dc(n),j_dc(n))
          do k=1,nlev
            write (unit_dc(n), '(a, i4, f19.10, 3e20.12)')  &
                           'in cupar: k,pr,uceml,dmeml,umeml= ',  &
                         k,Don_cape%model_p(i_dc(n),j_dc(n),k),    &
                          Don_conv%uceml(i_dc(n),j_dc(n),k), &
                           Don_conv%dmeml(i_dc(n),j_dc(n),k),  &
                           Don_conv%umeml(i_dc(n),j_dc(n),k)
          end do
        end do
      end if

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          do k=kstart_diag,nlev
            write (unit_dc(n), '(a, i4, e20.12)')  &
                              'in donner_deep: k,cuql', &
                              k,Don_conv%cuql (i_dc(n),j_dc(n),k)
          end do
        end do
      endif

!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          do k=kstart_diag,nlev
            if (abs(cemetf(igl_dc(n),jgl_dc(n),k     )) .gt. .002) then
              write (unit_dc(n), '(a, i4, e20.12)')  &
                             'in donner_deep: k, cemetf= ',k,   &
                              cemetf(igl_dc(n),jgl_dc(n),k)
              write (unit_dc(n), '(a, i4, e20.12)')  &
                               'in donner_deep: k, cual= ',k,    &
                               Don_conv%cual(i_dc(n),j_dc(n),k )
              write (unit_dc(n), '(a, i4, e20.12)')  &
                               'in donner_deep: k, xcape= ',k,    &
                               xcape_lag(igl_dc(n),jgl_dc(n)) 
              write (unit_dc(n), '(a, i4, e20.12)')   &
                                'in donner_deep: k, dcape = ',k,    &
                                Don_conv%dcape(i_dc(n),j_dc(n))
              write (unit_dc(n), '(a, i4, e20.12)')  &
                                'in donner_deep: k,a1    = ',k,    &
                                 Don_conv%a1 (i_dc(n),j_dc(n))
              write (unit_dc(n), '(a, i4, e20.12)')   &
                                'in donner_deep: k, amax  = ',k,   &
                                 Don_conv%amax(i_dc(n),j_dc(n))
            end if
          end do
        end do
      endif



end subroutine cupar



subroutine generate_cape_sounding (temp, mixing_ratio, press, Don_cape)      

real, dimension(:,:,:), intent(in) :: temp, mixing_ratio , press
type(donner_cape_type), intent(inout) :: Don_cape

!-------------------------------------------------------------------
!!    precu_vect defines pressure, temperature and moisture values on
!     generate_cape_sounding defines pressure, temperature and moisture values on
!     an enhanced vertical grid that will be used in the deep cumulus
!     convection parameterization (Donner, 1993, JAS).
!-------------------------------------------------------------------
!     On Input:
!
!        ilon          longitude index
!        jlat          latitude index
!        nlev          Cape.F resolution (parameter)
!        temperature   temperautre at GCM full levels (K)
!                      For third array index, 1 nearest top.
!
!     On Output:
!
!        cape_p         pressure at Cape.F resolution  (Pa)
!                      Index 1 nearest surface.
!        rcape         vapor mixing ratio at Cape.F resolution (kg(H2O)/kg)
!                      Index 1 nearest surface.
!        rini          vapor mixing ratio at Skyhi resolution (kg(H2O)/kg)
!                      Index 1 nearest surface.
!        tcape         temperautre at Cape.F resolution (K)
!        tini          temperaute at Skyhi resolution (K)
!                      Index 1 nearest surface.
!        prini         pressure at Skyhi resolution (Pa)
!                      Index 1 nearest surface.
!!
!

   real, dimension (size(temp,1), size(temp,2) ) :: dp

integer    :: n, k
      integer :: j, i
!     integer :: isize, jsize

!     isize = size (temperature  ,1)
!     jsize = size (temperature  ,2)

       if ( (conv_calc_on_this_step)  ) then
         if (in_diagnostics_window) then
           do n=1,ncols_in_window
             do k=1,nlev
               write (unit_dc(n), '(a, i4, f20.14, e20.12, e15.5)') &
                               'calculate_cape input profiles:&
                            & k,temp, mixing ratio, pressure',  k,  &
                             temp(i_dc(n),j_dc(n),k), &
                            mixing_ratio(i_dc(n),j_dc(n),k), &
                            press(i_dc(n), j_dc(n), k)
             end do
           end do
         end if
       end if
!--------------------------------------------------------------------
!   create arrays of moisture, temperature (K) and pressure with 
!   index 1 being closest to the surface. require that the moisture 
!   value be  non-negative.
!--------------------------------------------------------------------
      do k=1,nlev
        do j=1,jsize       
          do i=1,isize      
            Don_cape%model_r (i,j,nlev+1-k) =  &
                                   amax1 (mixing_ratio (i,j,k), 0.0e00)
            Don_cape%model_t (i,j,nlev+1-k) = temp  (i,j,k)
            Don_cape%model_p(i,j,nlev+1-k) = press(i,j,k)
          end do
        end do
      end do

!-------------------------------------------------------------------
!   define the vertical resolution of the convection parameterization
!   grid. define the top level pressure in that grid to be zero.
!   interpolate to define the pressure levels of that grid.
!-------------------------------------------------------------------
      do j=1,jsize        
        do i=1,isize      
          dp(i,j) = (Don_cape%model_p(i,j,1) -   &
                     Don_cape%model_p(i,j,nlev))/(ncap-1)
          Don_cape%cape_p(i,j,ncap) = 0.
        end do
      end do
      do k=1,ncap-1
        do j=1,jsize         
          do i=1,isize       
            Don_cape%cape_p(i,j,k) = Don_cape%model_p(i,j,1) - &
                                     (k-1)*dp(i,j)
          end do
        end do
      end do

!--------------------------------------------------------------------
!   call polat_vect to define values of temperature and moisture on
!   the enhanced convection grid by interpolating from the model grid
!   values. insure that the moisture field is positive-definite.
!--------------------------------------------------------------------

       call polat_vect (Don_cape%model_t, Don_cape%model_p,   &
                        Don_cape%cape_p, Don_cape%env_t           )
       call polat_vect (Don_cape%model_r, Don_cape%model_p,  &
                        Don_cape%cape_p, Don_cape%env_r           )
       do k=1,ncap
         do j=1,jsize        
           do i=1,isize        
               Don_cape%env_r(i,j,k) = MAX(Don_cape%env_r(i,j,k), 0.0)
           end do
         end do
       end do

!---------------------------------------------------------------------
!   if debugging is activated, print out the input pressure, moisture 
!   and temperature fields in the  desired column.
!---------------------------------------------------------------------
      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          do k=1,nlev-kstart_diag+1
            write(unit_dc(n), '(a, i4, f19.10, f20.14, e20.12)') &
                       'in precu: k, model_p,model_t,model_r = ', k,   &
                         Don_cape%model_p(i_dc(n),j_dc(n),k),   &
                         Don_cape%model_t(i_dc(n),j_dc(n),k), &
                         Don_cape%model_r(i_dc(n),j_dc(n),k)
          end do
        end do
      endif

!--------------------------------------------------------------------

end subroutine generate_cape_sounding


!###################################################################

subroutine polat_vect (xv, pv, p, x, exit_flag)

!------------------------------------------------------------------
!    polat_vect interpolates the field xv on a  pressure grid pv to
!    a field x on a pressure grid p.
!------------------------------------------------------------------
 
real, dimension(:,:,:), intent(in)   :: xv, pv, p
real, dimension(:,:,:), intent(out)  :: x
logical, dimension(:,:), intent(in), optional ::exit_flag

!     ON INPUT:
!
!         XV(N)   DATA AT RESOLUTION N
!         PV(N)   PRESSURE AT N LEVELS
!
!     ON OUTPUT:
!
!         X       DATA AT PRESSURE P
!
!       PARAMETER(N=40      ,NM1=N-1)
!       DIMENSION XV(iminp:imaxp,jminp:jmaxp,N),    &
! PV(iminp:imaxp,jminp:jmaxp,N)
!       dimension x (iminp:imaxp, jminp:jmaxp,ncap)
!       dimension p (iminp:imaxp, jminp:jmaxp,ncap)

      real :: xs
!integer, dimension(iminp:imaxp, jminp:jmaxp)  :: kstart
      integer      :: kstart, kkstart
      real    :: xt1, xt2
      real, dimension(size(pv,3)) :: pv_1d, xv_1d, dxdpv
      real, dimension(size(p,3)) :: p_1d
      real, dimension(size(x,1), size(x,2), size(x,3)) :: x2

      integer :: isize, jsize
      integer :: n, j, i, k, kk, nc

       n = size(pv,3)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        isize = size(x,1)
        jsize = size(x,2)

        do j=1,jsize       
         do i=1,isize        
           if (present(exit_flag) .and. exit_flag(i,j)) cycle
           pv_1d(:) = pv(i,j,:)
           xv_1d(:) = xv(i,j,:)
           p_1d(:)  = p(i,j,:)
           dxdpv(1:n-1) = (xv_1d(2:n)-xv_1d(1:n-1))/(pv_1d(2:n) - &
                                                     pv_1d(1:n-1))
           kkstart = 1
           do k=1,ncap
            IF (p_1d(    k) .GE. PV_1d(    1)) THEN
               xt2      =dxdpv(1)*(P_1d(k)-PV_1d(1))+XV_1d(1)
            else IF (P_1d(k) .LE. PV_1d(n)) THEN
                Xt2     =dxdpv(n-1) *(P_1d(k)-PV_1d(n))+XV_1d(n)
            else
              DO  kk=kkstart     ,N-1
                 IF ((P_1d(k) .GE. PV_1d(kk+1)) ) then
                   xt2   =dxdpv(kk) *(P_1d(k)-PV_1d(kk+1))+XV_1d(kk+1)
                    kkstart      = kk
                   exit
                 ENDIF
              end do
            endif
            x(i,j,k) = xt2
           end do
         end do
        end do

        if (in_diagnostics_window) then
          do nc=1,ncols_in_window
            do k=1,ncap
              write(unit_dc(nc), '(a, i4, f19.10, f20.14)')  &
                    'in polat: k,p,x=', k, p(i_dc(nc),j_dc(nc),k),    &
                     x(i_dc(nc),j_dc(nc),k)
            end do
          end do
        endif




end subroutine polat_vect









subroutine calculate_cape (is, ie, js, je, Cape, exit_flag, perturbed)
 
!----------------------------------------------------------------------
!    calculate_cape calculates convective available potential energy for a 
!    cloud whose temperature follows a saturated adiabat.
!---------------------------------------------------------------------

integer :: is, ie, js, je
type(donner_cape_type), intent(inout) :: Cape      
logical, dimension(:,:), intent(in), optional  :: exit_flag
logical, intent(in), optional  :: perturbed

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      p_v     pressure (Pa)
!              Index 1 refers to level nearest earth's surface.
!      r_v     mixing ratio (kg(H2O)/kg)
!              Index 1 refers to level nearest earth's surface.
!      t_v     temperature (K)
!              Index 1 refers to level nearest earth's surface.
!
!   intent(out) variables:
!
!     cin_v    convective inhibition (J/kg)
!              energy required to lift parcel from level istart to
!              level of free convection
!     plcl_v   pressure at lifting condensation level (Pa)
!     plfc_v   pressure at level of free convection (Pa)
!              height of plfc .le. height of plcl
!              If parcel becomes buoyant below plcl, cin can be .lt. 0
!     plzb_v   pressure at level of zero buoyancy (Pa)
!     rpc_v    parcel mixing ratio
!              Set to environment below istart.
!              Index 1 refers to level nearest earth's surface.
!     tpc_v    parcel temperature (K)
!              Set to environment below istart.
!              Index 1 refers to level nearest earth's surface.
!     xcape_v  convective available potential energy (J/kg)
!              energy released as parcel moves from level of free
!              convection to level of zero buoyancy
!
!     Calculated:
!   
!     tot      xcape+cin (J/kg)
!--------------------------------------------------------------------
!

   integer, dimension(size(Cape%env_t,1),   &
                      size(Cape%env_t,2) ) :: klcl, klzb, klfc, ieqv      

   logical, dimension(size(Cape%env_t,1),   &
                      size(Cape%env_t,2) ) :: capepos, cape_exit,&
                                              lcl_found, lzb_found,  &
                                              skip_search, &
                                              cin_done, cape_skip

      real, dimension(size(Cape%env_t,1),   &
                      size(Cape%env_t,2) ) :: tot, ro, tc_1, tp, &
                                  q_ve, cp_v, rs_v, es_v, tlcl, rlcl, &
                                  dt_v, dtdp_v, tc_v, qe_v, qs_v, &
                                  tve_v, tv_v, tvc_v,        delt_v

     integer :: s_found, z_found, cin_counter, error_flag
     real    :: q_ve_s, cp_v_s,               tp_s, es_v_s, rs_v_s, &
                dtdp_v_s, dt_v_s, rlcl_s, tlcl_s, plcl_v_s, p_v_s, &
                p_v_sm, p_v_sp
     integer :: ieqv_s, klcl_s

     real :: tc_v_s, qe_v_s, tve_v_s, qs_v_s, tv_v_s
     real :: tvc_v_s, delt_v_s
     real, dimension(ncap)  :: tpc_v_s, rpc_v_s
     real, dimension(size(Cape%cape_p,3)) :: p_v_k, r_v_k, t_v_k
     integer :: klzb_s, klfc_s
     logical :: capepos_s, cape_exit_s
     real  :: plfc_v_s, plzb_v_s, cin_v_s, xcape_v_s, tot_s
     integer :: j, i
     real    :: rc, pb, fact1, fact2, fact3, dtdp, rbc, rbe, qc, qe
!--------------------------------------------------------------------
!     Stop CAPE calculation when pressure falls to pstop or
!     parcel temperature falls below tmin.
!     istart-index of level whose mixing ratio is conserved as a parcel 
!            leaves it undergoing dry adiabatic ascent
!--------------------------------------------------------------------
      real ::   pstop=4.0e03
      real ::   tmin=154.
      integer :: istart=1
      integer :: isize, jsize
      integer  :: n, k


!--------------------------------------------------------------------
    if (in_diagnostics_window) then
      do n=1,ncols_in_window
        write(unit_dc(n), '(a, f14.5)') 'in cape: cpi= ', cpi
        write(unit_dc(n), '(a, f14.5)') 'in cape: cpv= ', cpv
        write(unit_dc(n), '(a, f14.5)') 'in cape: rocp= ', rocp
        write(unit_dc(n), '(a, f14.5)') 'in cape: rair= ', rair
        write(unit_dc(n), '(a, f14.5)') 'in cape: cpi= ', latvap
        write(unit_dc(n), '(a, f14.5)') 'in cape: rvap= ', rvap
        do k=1,ncap
          write (unit_dc(n), '(a, i4, f19.10, f20.14, e20.12)')   &
                           'press, temp, vapor in cape: k, p,t,r = ',  &
                            k, Cape%cape_p  (i_dc(n),j_dc(n),k), &
                            Cape%env_t(i_dc(n), j_dc(n),k),   &
                             Cape%env_r(i_dc(n),j_dc(n),k)  
        end do
      end do
    endif

      isize = size(Cape%env_t,1)
      jsize = size(Cape%env_t,2)

!-------------------------------------------------------------------
!     loop over columns in physics window.
!------------------------------------------------------------------
      do j=1,jsize       
        do i=1,isize        

!------------------------------------------------------------------
!  if exit_flag is .true., it is already known that the cape 
!  calculation is not needed in this column. cycle to next column.
!-------------------------------------------------------------------
          if (present (exit_flag) .and. exit_flag(i,j)) cycle

!------------------------------------------------------------------
!  if lowest level temp is lower than tmin, set output fields in this
!  column to defaults and cycle to next column.
!------------------------------------------------------------------
!         if (t_v(i,j,istart) .lt. tmin)  then
          if (Cape%env_t(i,j,istart) .lt. tmin)  then
!RETURN
            do k=1,ncap
!             rpc_v(i,j,k) = r_v(i,j,k)
              Cape%parcel_r(i,j,k) = Cape%env_r(i,j,k)
!             tpc_v(i,j,k) = t_v(i,j,k)
              Cape%parcel_t(i,j,k) = Cape%env_t(i,j,k)
            end do 
            Cape%plfc(i,j)=pstop       
            Cape%plzb(i,j)=pstop      
            Cape%plcl(i,j)=pstop      
            Cape%coin(i,j) = 0.       
!           Don_cape%coin_vect(i,j) = 0.       
            Cape%xcape(i,j) = 0.             
            cycle
          endif
!-------------------------------------------------------------------
!     initialize variables.
!------------------------------------------------------------------
          plfc_v_s   =pstop
          plzb_v_s   =pstop
          plcl_v_s   =pstop

          klfc_s = ncap - 1
          klcl_s = ncap - 1

          cin_v_s = 0.
          xcape_v_s  = 0.
          tot_s = 0.

          capepos_s = .false.
          cape_exit_s = .false.
!--------------------------------------------------------------------
!  define parcel departure point values. convert mixing ratio to 
!  specific humidity.
!--------------------------------------------------------------------
!         tpc_v_s(istart) = t_v(i,j,istart)
          tpc_v_s(istart) = Cape%env_t(i,j,istart)
!  rpc_v_s(istart) = r_v(i,j,istart)
  rpc_v_s(istart) = Cape%env_r(i,j,istart)
          tp_s   =tpc_v_s(istart)
          q_ve_s    = rpc_v_s(istart)/( 1.+rpc_v_s(istart) )
          cp_v_s    = cpi*(1.+((cpv/cpi)-1.)*q_ve_s   )

!-------------------------------------------------------------------
!  define 1-d arrays of pressure, temperature and mixing ratio in the 
!  column.
!-------------------------------------------------------------------
!         p_v_k(:) = p_v(i,j,:)
          p_v_k(:) = Cape%cape_p(i,j,:)
!         t_v_k(:) = t_v(i,j,:)
          t_v_k(:) = Cape%env_t(i,j,:)
!         r_v_k(:) = r_v(i,j,:)
          r_v_k(:) = Cape%env_r(i,j,:)

!--------------------------------------------------------------------
!  move the parcel upwards to find the lcl in the column.
!--------------------------------------------------------------------
          do k=istart,ncap  ! k loop to find lcl
!           p_v_k_s = p_v_k(k)
!--------------------------------------------------------------------
!  if the temperature and pressure are still within limits, continue 
!  parcel movement. determine saturation mixing ratio for parcels at 
!  this level.
!---------------------------------------------------------------------
!           if (tp_s >= tmin .and. p_v_k(k) >= pstop) then
            if (tp_s >= tmin .and. p_v_k(k) >= pdeep_cv) then
!           if (tp_s >= tmin .and. p_v_k_s  >= pstop) then
!             call escomp (tp_s, es_v_s)
              call lookup_es (tp_s, es_v_s)
!!! THIS FORMULA IS INCORRECT -- no rocp*es in denom -- mixing ratio
!!!  is in use here.
              rs_v_s  = rocp*es_v_s /(p_v_k(k)  + (rocp - 1.)*es_v_s)
!      denom = p_v_k_s + (rocp-1.0)*es_v_s
!             rs_v_s  = rocp*es_v_s /(p_v_k_s   + (rocp - 1.)*es_v_s)
!             rs_v_s  = rocp*es_v_s /denom                            

!--------------------------------------------------------------------
!  check if the parcel is now saturated.
!---------------------------------------------------------------------
              call ieq_x (rs_v_s, rpc_v_s(istart), ieqv_s)

!--------------------------------------------------------------------
!   if saturation is exact or if parcel is super-saturated at its 
!   starting level, save pressure, temp, mixing ratio and cloud base 
!   level. exit vertical loop.
!--------------------------------------------------------------------
              if ( (ieqv_s ==  0) .or. &
                   (ieqv_s < 0 .and. k == istart) ) then
                plcl_v_s    = p_v_k(k)  
!               plcl_v_s    = p_v_k_s   
                rlcl_s    = r_v_k(  k)
                tlcl_s    = t_v_k(  k)
                klcl_s    = k
                exit
              endif

!--------------------------------------------------------------------
!   if parcel is super-saturated, define cloud-base pressure, temp 
!   and mixing ratio as the average value between the current level
!   and the next lower level, and this level as the cloud base level.
!   exit the column.
!--------------------------------------------------------------------
              if (ieqv_s < 0) then
                plcl_v_s  = (p_v_k(k)     +p_v_k(k-1)  )/2.
!               plcl_v_s  = (p_v_k_s      +p_v_k(k-1)  )/2.
                tlcl_s    = (t_v_k(  k)+t_v_k(  k-1))/2.
                rlcl_s    = (r_v_k(  k)+r_v_k(  k-1))/2.
                klcl_s    = k
                exit

!---------------------------------------------------------------------
!    if the parcel remains unsaturated at this level and the top of the 
!    model has not been reached, move parcel along dry adiabat to next 
!    pressure level. define temperature at this level; verify that it 
!    is warmer than tmin. save parcel temperature and mixing ratio at 
!    this next higher level.
!---------------------------------------------------------------------
              else  ! (ieqv_s    .gt. 0) 
                if (k .lt. ncap) then
                  dtdp_v_s = rair*tp_s/cp_v_s   
                  dt_v_s   = dtdp_v_s*alog( p_v_k(k+1)/p_v_k(k) )
!                 dt_v_s   = dtdp_v_s*alog( p_v_k(k+1)/p_v_k_s  )
                  tp_s   = tp_s  + dt_v_s    
                  if (tp_s    .lt. tmin)  then
                    cape_exit_s    = .true.
!RETURN
                    tpc_v_s(k+1:ncap) = t_v_k(  k+1:ncap)
                    rpc_v_s(k+1:ncap) = r_v_k(  k+1:ncap)
                    exit
                  else  
                    tpc_v_s(  k+1)=tp_s   
                    rpc_v_s(  k+1)=rpc_v_s(istart)
                  endif

!-------------------------------------------------------------------
!    if have reached top of model, set flag to stop integration in this
!    column.
!-------------------------------------------------------------------
                else
                  cape_exit_s    = .true.
                endif
              endif ! (ieqv_s    .gt. 0)

!--------------------------------------------------------------------
!    if either parcel temperature or pressure is below cutoff values,
!    set remainder of parcel sounding to the environment and stop 
!    searching in this column.
!--------------------------------------------------------------------
           else
!RETURN
             tpc_v_s(k+1:ncap) = t_v_k(  k+1:ncap)
             rpc_v_s(k+1:ncap) = r_v_k(  k+1:ncap)
             cape_exit_s = .true.
           endif
         end do   ! k loop to find lcl

!---------------------------------------------------------------------
!   if lcl not found, stop calculations in this column.
!---------------------------------------------------------------------
         if (cape_exit_s) go to 12

!--------------------------------------------------------------------
!   if in debug mode, print out info on lcl in debug column.
!--------------------------------------------------------------------
           if (in_diagnostics_window ) then
             do n=1,ncols_in_window
               if (i == i_dc(n) .and. j == j_dc(n)) then
                 write (unit_dc(n), '(a, f19.10,i4, f20.14, e20.12)')  &
                           'in cape: plcl,klcl,tlcl,rlcl= ',   &
                             plcl_v_s, klcl_s, tlcl_s, rlcl_s       
                 write (unit_dc(n), '(a, f19.10)')   &
                              'in cape: p(klcl)= ',  p_v_k(klcl_s)
               endif
             end do
           endif

!-------------------------------------------------------------------
!   calculate temperature along saturated adiabat, starting at p(klcl)
!   and a temperature tp to find the level of free convection and
!   the level of zero buoyancy. 
!--------------------------------------------------------------------
         do k=klcl_s,ncap-1

!--------------------------------------------------------------------
!   search for the lfc only up to a pressure pstop.
!--------------------------------------------------------------------
           if ( (p_v_k(k+1)   >= pstop) .or.   &
                (p_v_k(k+1)   >= pdeep_cv .and. .not. capepos_s) ) then

!--------------------------------------------------------------------
!   define saturation vapor pressure for the parcel.
!--------------------------------------------------------------------
!            call escomp (tp_s, es_v_s)
             call lookup_es (tp_s, es_v_s)

!--------------------------------------------------------------------
!   define the environmental and parcel virtual temperature and specific
!   humidity.
!--------------------------------------------------------------------
!     denom = 1. + r_v_k(k)
             qe_v_s = r_v_k(k)/(1.+r_v_k(  k))
!            qe_v_s = r_v_k(k)/denom             
             tve_v_s = t_v_k(k)*(1.+.61*qe_v_s   )
             rs_v_s = rocp*es_v_s   /(p_v_k(k)  +   &
!            denom  = p_v_k_s + (rocp-1.)*es_v_s
!            rs_v_s = rocp*es_v_s   /(p_v_k_s   +   &
                         (rocp  -1.)*es_v_s   )
!            rs_v_s = rocp*es_v_s   /denom             
             qs_v_s = rs_v_s/(1.+rs_v_s   )
             tv_v_s = tp_s*(1.+.61*qs_v_s   )

!--------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
             call ieq_x (tv_v_s, tve_v_s, ieqv_s)

!---------------------------------------------------------------------
!   integrate parcel upward, finding level of free convection and 
!   level of zero buoyancy.
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!   determine if the level of free convection has been reached. 
!-------------------------------------------------------------------
             if ((ieqv_s >= 0) .and. (.not. capepos_s   )) then
               capepos_s    = .true.
               plfc_v_s    = p_v_k(k)  
!              plfc_v_s    = p_v_k_s   
               klfc_s    = k
             end if

!-------------------------------------------------------------------
!   determine if the level of zero buoyancy has been reached.  if so,
!   set flag so that calculation will be ended in this column.
!-------------------------------------------------------------------
             if ((ieqv_s    .lt. 0) .and. (capepos_s   )) then
               klzb_s    = k
               plzb_v_s    = (p_v_k(k)  +p_v_k(k-1)  )/2.
!              plzb_v_s    = (p_v_k_s   +p_v_k(k-1)  )/2.
               tpc_v_s(k+1:ncap) = t_v_k(  k+1:ncap)
               rpc_v_s(k+1:ncap) = r_v_k(  k+1:ncap)
               exit
!-----------------------------------------------------------------
!   if not, continue moving parcel up pseudo-adiabat to next cape-
!   calculation pressure level. define new parcel temperature and
!   mixing ratio at this level; if temperature is colder than allowed,
!   end integration.
!-------------------------------------------------------------------
             else   !  (cape is pos, parcel warmer than env)
               rc = (1.-qs_v_s   )*rair+qs_v_s   *rvap
               pb = 0.5*(p_v_k(k)     +p_v_k(k+1))
!              pb = 0.5*(p_v_k_s      +p_v_k(k+1))
              fact1 = rair/cpi
              fact2 = tv_v_s   +(latvap*qs_v_s   /rc)
              fact1 = fact1*fact2
              fact3 = rocp  *(latvap**2)*es_v_s   /    &
                      (cpi*pb*rvap*(tv_v_s   **2))
              fact3 = 1.+fact3
              dtdp = fact1/fact3
              tp_s    = tp_s   +dtdp*     &
                        alog(p_v_k(k+1)  /p_v_k(k)  )
!                       alog(p_v_k(k+1)  /p_v_k_s   )
              if (tp_s    .lt. tmin        )  then
                cape_exit_s    = .true.
!RETURN
                tpc_v_s(k+1:ncap) = t_v_k(  k+1:ncap)
                rpc_v_s(k+1:ncap) = r_v_k(  k+1:ncap)
                exit  ! exit k loop
              else
                tpc_v_s(  k+1) = tp_s   
                rpc_v_s(  k+1) = rs_v_s   
              endif
            endif   !  (ieq < 0, capepos)

!--------------------------------------------------------------------
!    if pressure has gone below the minimum at which deep convection 
!    is allowed, set flag to end calculation in this column.
!--------------------------------------------------------------------
            if (in_diagnostics_window ) then
              do n=1,ncols_in_window
                if (i == i_dc(n) .and. j == j_dc(n)) then
                  write (unit_dc(n), '(a, i4, 2f20.14)')  &
                                  'in cape: k,tv,tve= ',k,tv_v_s,  &
                                   tve_v_s            
                  write (unit_dc(n), '(a, i4, 2f19.10)')  &
                                  'in cape: klzb,plzb,p(klzb)= ',  &
                                   klzb_s            ,   &
                                  plzb_v_s            ,   &
                                 p_v_k(             klzb_s            )
                  write (unit_dc(n), '(a, 3f17.10)')  &
                            'in cape: fact1,fact2,rc= ',fact1,fact2,rc
                  write (unit_dc(n), '(a, 2f17.10)') &
                           'in cape: fact1,fact3= ',fact1,fact3
                  write (unit_dc(n), '(a, f17.10)')  &
                            'in cape: dtdp= ',dtdp
                  write (unit_dc(n), '(a,  2f20.14)') &
                             'in cape: tc,t= ',tp_s      ,  &
                               t_v_k(           k+1)
                  write (unit_dc(n), '(a, f19.10, 2e20.12)') &
                          'in cape: p,r,rs= ',p_v_k(           k+1),   &
                              r_v_k(           k+1),  rs_v_s            
                endif
              end do
            endif
          else
!RETURN
            tpc_v_s(k+1:ncap) = t_v_k(  k+1:ncap)
            rpc_v_s(k+1:ncap) = r_v_k(  k+1:ncap)
            exit
          endif
        end do   ! k loop
        if (cape_exit_s   ) go to 12

!--------------------------------------------------------------------
!   if this was a call to calculate perturbd profile, bypass cape and
!   cin calculation, since only the tpc and rpc profiles are needed.
!--------------------------------------------------------------------
!     if (present (exit_flag)) go to 12
      if (present (perturbed)) go to 13

!-------------------------------------------------------------------
!   define flag to indicate where to calculate convective inhibition.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!   calculate convective inhibition.
!--------------------------------------------------------------------
      do k=istart,ncap-1
!        p_v_k_s = p_v_k(k)

!              if (in_diagnostics_window .and. &
!           i == itest .and. j == j_dc) then
!              else
!             if ( p_v_k(k+1)   <= pstop              ) then
!               cape_exit_s    = .true.
!exit
!             endif
!             endif
!------------------------------------------------------------------
!   determine if sounding fails to produce a level of free convection.
!   if so, set flag to avoid cape calculation. If desired, print out
!   columns where lcl exists, but no lfc.
!------------------------------------------------------------------
!             if ( p_v_k(k+1)   <= pstop              ) then
              if ( p_v_k(k+1)   <= pdeep_cv           ) then
                cape_exit_s    = .true.
!               if (debug) then
!                 print *, 'cape = 0 (NO LFC): i, jrow, cin= ',   &
!   i, j+js-1 , cin_v_s   
!endif
                exit   ! exit k loop
              end if

!--------------------------------------------------------------------
!    define the specific humidity and virtual temperature of the
!    parcel and environment.
!--------------------------------------------------------------------
      rbc =     (rpc_v_s(  k) + rpc_v_s(  k+1))/2.
      rbe =     (r_v_k(  k) + r_v_k(  k+1))/2.
      qc = rbc/(1. + rbc)
      qe = rbe/(1. + rbe)
      tvc_v_s    =     (tpc_v_s(  k) + tpc_v_s(k+1))/2.
      tve_v_s    =     (t_v_k(  k) + t_v_k(  k+1))/2.
      tvc_v_s    = tvc_v_s *  (1.+.61*qc)
      tve_v_s    = tve_v_s   *(1.+.61*qe)

      if (in_diagnostics_window ) then
        do n=1,ncols_in_window
          if (i == i_dc(n) .and. j == j_dc(n)) then
            write (unit_dc(n), '(a, i4, 2f20.14)')  &
                            'in cape: k,tvc,tve= ', k,  &
                         tvc_v_s,  tve_v_s            
          endif
         end do
       endif

!---------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
        call ieq_x(tvc_v_s   , tve_v_s   , ieqv_s   )

!---------------------------------------------------------------------
!   add the contribution to cin from this pressure layer.
!---------------------------------------------------------------------
      if ((ieqv_s    .lt. 0) .or.      &
          (p_v_k(k)   .gt. plfc_v_s   ))  then
!         (p_v_k_s    .gt. plfc_v_s   ))  then
        delt_v_s    = rair*(tvc_v_s   -tve_v_s   )*   &
                      alog(p_v_k(k) /p_v_k(k+1))
!                     alog(p_v_k_s  /p_v_k(k+1))
        cin_v_s    = cin_v_s    - delt_v_s   
      else

!------------------------------------------------------------------
!   determine if sounding fails to produce a level of free convection.
!   if so, set flag to avoid cape calculation. If desired, print out
!   columns where lcl exists, but no lfc.
!------------------------------------------------------------------
!               if (p_v_k(k) < plfc_v_s .or.  &
!                   p_v_k(k) <= pstop) then
!          cape_exit_s    = .true.
!               if (debug) then
!                 print *, 'cape = 0 (NO LFC): i, jrow, cin= ',   &
!   i, j+js-1 , cin_v_s   
!endif
        exit
!endif
      end if


    end do  ! k loop

      if (column_diagnostics_desired .and. cape_exit_s) then
        print *, 'cape = 0 (NO LFC): i, jrow, cin= ',   &
                                i, j+js-1 , cin_v_s   
      endif

      if (cape_exit_s    ) go to 12



!-------------------------------------------------------------------
!   if desired, print out lfc k index and pressure.
!-------------------------------------------------------------------
      do n=1,ncols_in_window
        if (i == i_dc(n) .and. j == j_dc(n)) then
           write (unit_dc(n), '(a, i4, f19.10)')  &
                         'in cape: klfc, p(klfc)= ', klfc_s       ,  &
                         p_v_k(             klfc_s            )
        endif
      end do

!--------------------------------------------------------------------
!  calculate convective available potential energy.
!--------------------------------------------------------------------
      do k=klfc_s,ncap-1
!        p_v_k_s = p_v_k(k)

!--------------------------------------------------------------------
!  define flag to indicate which columns are actively computing cape.
!-------------------------------------------------------------------
        if ( p_v_k(k+1) .gt. plzb_v_s    ) then

!--------------------------------------------------------------------
!  define virtual temperature and specific humidity of parcel and 
!  environment.
!-------------------------------------------------------------------
          rbc = (rpc_v_s(  k)+rpc_v_s(  k+1))/2.
          rbe = (r_v_k(  k)+r_v_k(  k+1))/2.
          qc = rbc/(1.+rbc)
          qe = rbe/(1.+rbe)
          tvc_v_s    = (tpc_v_s(  k)+tpc_v_s(  k+1))/2.
          tve_v_s    = (t_v_k(  k)+t_v_k(  k+1))/2.
          tvc_v_s    = tvc_v_s   *(1.+.61*qc)
          tve_v_s    = tve_v_s   *(1.+.61*qe)

!--------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
          call ieq_x(tvc_v_s   , tve_v_s   , ieqv_s              )

!---------------------------------------------------------------------
!   add the contribution to column cape from this pressure layer.
!---------------------------------------------------------------------
          if (ieqv_s    .ge. 0) then
            delt_v_s    = rair*(tvc_v_s   -tve_v_s   )*    &
                          alog(p_v_k(k)/p_v_k(k+1))
!                         alog(p_v_k_s /p_v_k(k+1))
            xcape_v_s   =xcape_v_s   +delt_v_s   
          end if

!---------------------------------------------------------------------
!   print out cape and cape contribution from this level.
!---------------------------------------------------------------------
          if (in_diagnostics_window) then
            do n=1,ncols_in_window
              if ( i == i_dc(n) .and. j == j_dc(n)) then
                if (ieqv_s             .ge. 0) then
                  write (unit_dc(n), '(a,i4, 2f12.8)')  &
                                  'in cape: k,delt,xcape= ',k,   &
                                    delt_v_s            ,    &
                                     xcape_v_s            
                endif
              endif
            end do
          endif
        else
          exit
        endif
      end do  ! end of k loop

!--------------------------------------------------------------------
!  print out diagnostics (cape, cin, tot), if desired.
!--------------------------------------------------------------------
      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          if ( i == i_dc(n) .and. j == j_dc(n)) then
            if (cin_v_s    /= 0.0 .or. &
                xcape_v_s    /= 0.0) then
              tot_s             = xcape_v_s             -   &
                                  cin_v_s             
              write (unit_dc(n), '(a, f12.6, a)')  &
                           'in cape: cin= ',cin_v_s            ,' J/kg'
              write (unit_dc(n), '(a, f12.6, a)')  &
                           'in cape: xcape= ',xcape_v_s,  &
                            ' J/kg'
              write (unit_dc(n), '(a, f12.6, a)')  &
                            'in cape: tot= ',tot_s            ,' J/kg'
            endif
          endif
        end do
      endif

!--------------------------------------------------------------------
!  check for error in cape calculation. stop execution if present.
!--------------------------------------------------------------------
      if (xcape_v_s    .lt. 0.) then            
        call error_mesg ( 'calculate_cape',  &
                  ' xcape error -- value < 0.0 ', FATAL)
      endif


!--------------------------------------------------------------------
!  fill appropriate (i,j) indices in output variables.
!--------------------------------------------------------------------
   12 continue
      Cape%plfc(i,j)=plfc_v_s
      Cape%plzb(i,j)=plzb_v_s
      Cape%plcl(i,j)=plcl_v_s
      Cape%coin(i,j) = cin_v_s
!     Don_cape%coin_vect(i,j) = cin_v_s
      Cape%xcape(i,j) = xcape_v_s
  13  continue
!RETURN
!       if (.not. present(exit_flag) .or. present(perturbed)) then
!         if ( .not. cape_exit_s) then
      do k=1,ncap
        Cape%parcel_r(i,j,k) = rpc_v_s(k)
        Cape%parcel_t(i,j,k) = tpc_v_s(k)
      end do 
!          endif
!       endif

    end do ! i loop
    end do ! j loop


      if ( (conv_calc_on_this_step) .and.   &
            (.not. present(perturbed)) ) then

!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
         if (in_diagnostics_window) then
           do n=1,ncols_in_window
             write (unit_dc(n), '(a, 2f19.10)')  &
                               'in donner_deep: plfc,plzb= ',  &
                                Cape%plfc(i_dc(n),j_dc(n)),  & 
                                Cape%plzb(i_dc(n),j_dc(n))
             write (unit_dc(n), '(a, 3f19.10)')  &
                           'in donner_deep: plcl,coin,xcape= ',   &
                            Cape%plcl(i_dc(n),j_dc(n)),  &
                            Cape%coin(i_dc(n),j_dc(n)),   &
                            Cape%xcape(i_dc(n),j_dc(n))
!23 changed j to j_dc:
             do k=1,ncap      
               write (unit_dc(n), '(a, i4, f19.10)') &
                                 'in donner_deep: k,cape_p= ',k,  &
                                    Cape%cape_p(i_dc(n),j_dc(n),k)
               write (unit_dc(n), '(a, i4, 2f20.14)')  &
                                'in donner_deep: k,tcape,tpca= ',k,   &
                                 Cape%env_t(i_dc(n),j_dc(n),k),   &
                                 Cape%parcel_t(i_dc(n),j_dc(n),k)
               write (unit_dc(n), '(a, i4, 2e20.12)')  &
                                'in donner_deep: k,rcape,rpca= ',k,   &
                                 Cape%env_r(i_dc(n),j_dc(n),k),    &
                                 Cape%parcel_r(i_dc(n),j_dc(n),k)
               if (Cape%cape_p(i_dc(n),j_dc(n),k) <     &
                   Cape%plzb(i_dc(n),j_dc(n)))  exit
             end do
           end do
         endif

     endif

!---------------------------------------------------------------------






end subroutine calculate_cape







!#####################################################################

function iequ(x,y)

real,   intent(in)   :: x
real,   intent(in)   :: y
!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!

      integer :: iequ
      real :: eps, epsm, d

      iequ=0
!     eps=1.e-14
      eps=1.e-13
      epsm=-eps
      d=x-y
      if (d .gt. eps) iequ=10
      if (d .lt. epsm) iequ=-10
!      return



end function iequ

subroutine  ieq_y (x, y, ieq)

!-------------------------------------------------------------------
real, dimension(:,:),   intent(in)   :: x
real,                   intent(in)   :: y
integer,dimension(:,:), intent(out)  :: ieq
!-------------------------------------------------------------------
!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!
!      iequ=0
!!     eps=1.e-14
     real, dimension(size(x,1), size(x,2)) :: d
     real :: eps, epsm
     integer :: i, j

      eps=1.e-13
      epsm=-eps

      d(:,:) = x(:,:) - y

      do j=1,size(x,2)
      do i=1,size(x,1)

        if (d(i,j) .gt. eps) then
          ieq(i,j) = 10
        else if (d (i,j).lt. epsm) then
          ieq(i,j) = -10
        else 
          ieq(i,j) = 0
        endif
      end do
      end do



!     if (d .gt. eps) iequ=10
!     if (d .lt. epsm) iequ=-10
!      return



end subroutine ieq_y



subroutine  ieq_z (x, y, ieq, flag)

!-------------------------------------------------------------------
real, dimension(:,:),   intent(in)   :: x, y
integer,dimension(:,:), intent(out)  :: ieq
logical, dimension(:,:), intent(in)  :: flag
!-------------------------------------------------------------------
!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!
!      iequ=0
!!     eps=1.e-14
     real, dimension(size(x,1), size(x,2)) :: d
     real :: eps, epsm
     integer :: i, j

      eps=1.e-13
      epsm=-eps
      

!      d(:,:) = x(:,:) - y(:,:)

      do j=1,size(x,2)
      do i=1,size(x,1)

        if (.not. flag(i,j)) then
          d(i,j) = x(i,j) - y(i,j)
          if (d(i,j) .gt. eps) then
            ieq(i,j) = 10
          else if (d (i,j).lt. epsm) then
            ieq(i,j) = -10
          else 
            ieq(i,j) = 0
          endif
        endif
      end do
      end do



end subroutine ieq_z






subroutine  ieq_x (x, y, ieq)

real,                   intent(in)   :: x, y
integer,                intent(out)  :: ieq

!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!
       real :: eps, epsm, d
8      integer :: i, j
!      iequ=0
!!     eps=1.e-14

      eps=1.e-13
      epsm=-eps
      d=x-y
      if (d .gt. eps) then
        ieq = 10
      else if (d .lt. epsm) then
        ieq = -10
      else 
        ieq = 0
      endif
!     if (d .gt. eps) iequ=10
!     if (d .lt. epsm) iequ=-10
!      return



end subroutine ieq_x



!----------------------------------------------------------------------

subroutine cuclo_vect (dcape_v,         pcape_v, plfc_v, plzb_v,  &
                       qli0_v, qli1_v, qr_v, qt_v,  &
                        r_v, ri_v, rl_v, rpc_v, tcape_v, tpc_v, a1_v,  &
                        exit_flag, nine_flag, is, ie, js, je)

!---------------------------------------------------------------------
integer,  intent(in)   :: is, ie, js, je
real, dimension(:,:  ), intent(in) :: dcape_v, plfc_v, plzb_v
logical, dimension(:,:), intent(in) ::  exit_flag, nine_flag
real, dimension(:,:,:), intent(in) :: pcape_v, qli0_v, qli1_v,    &
                                      qr_v, qt_v, &
                                      r_v, ri_v, rl_v, rpc_v, tcape_v, &
                                      tpc_v
real, dimension(:,:), intent(out)  :: a1_v
!---------------------------------------------------------------------

!
!     Calculates a_1(p_b) for closing cumulus parameterization.
!     See LJD notes, "Cu Closure D," 6/11/97
!
!
!     On Input:
!
!     Parameters:
! 
!        nc        number of levels at cloud-model resolution
!
!     Variables:
!
!        cpd       specific heat of dry air at constant pressure [J/(kg K)]
!        cpv       specific heat of water vapor
!        dcape     rate of change of convective available potential
!                  energy due to large-scale processes [J/(kg s)]
!        epsilo    ratio of molecular weights of water vapor to air
!        p         pressure (Pa)
!                  level 1 nearest ground.
!        plfc      pressure at level of free convection (Pa)
!        plzb      pressure at level of zero buoyancy (Pa)
!        qli0      normalized component of cumulus condensate forcing
!                  [kg(H2O)/(kg sec)]
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 4.
!        qli1      un-normalized component of condensate forcing
!                  [kg(h2O)/(kg sec)]
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 4.
!        qr        normalized cumulus moisture forcing [kg(H2O)/(kg sec)]
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 1.
!        qt        normalized cumulus thermal forcing (K/sec)
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 1.
!        r         large-scale water-vapor mixing ratio (kg(H2O)/kg)
!                  level 1 nearest ground.
!        rd        gas constant for dry air [J/(kg K)]
!        ri        large-scale ice mixing ratio (kg(H2O)/kg)
!                  level 1 nearest ground.
!        rl        large-scale liquid mixing ratio (kg(H2O)/kg)
!                  level 1 nearest ground.
!        rlat      latent heat of vaporization (J/kg)
!        rpc       parcel vapor mixing ratio from Cape.F [kg(H2O)/kg]
!                  Index 1 at bottom of model.
!        rv        gas constant for water vapor (J/(kg K))
!        t         large-scale temperature (K)
!                  level 1 nearest ground.
!        tpc       parcel temperature from Cape.F (K)
!                  Index 1 at bottom of model.
!   
!     On Output:
! 
!         a1       fractional area of index-1 cu subensemble
!        

      real, dimension (ncap) :: p, qli0, qli1, qr, qt, r, ri, rl, t, &
                                rt, tpc, tpca, ta, ra, tden, tdena, &
                                rpc, rpca, dtpdta

      real, dimension (size(tcape_v,1), size(tcape_v,2), ncap) :: &
                                         ra_v, ta_v, rpca_v, tpca_v

   logical, dimension (size(tcape_v,1), size(tcape_v,2), ncap) :: &
                                                              exit_3d

   logical, dimension (size(tcape_v,1), size(tcape_v,2)      ) :: &
                                                               exit_2d

      real, dimension (size(tcape_v,1), size(tcape_v,2) ) :: &
                     cin_v2, plcl_v2, plfc_v2, plzb_v2, xcape_v2

      logical  :: debug_ijt
      logical  :: perturbed
      integer :: isize, jsize
      type(donner_cape_type) :: Cape_pert
      integer :: i, j, k, n     
      real  :: dcape, plfc, plzb, tdens, tdensa, ri1, ri2, rild, rile, &
               rilf, ri2b, sum2, rilak, rilbk, rilck, rilakm, rilbkm, &
               rilckm, rila, rilb, rilc, ri2ak, ri2akm, ri2a, sum1



       perturbed = .true.

      isize = size(tcape_v,1)
      jsize = size(tcape_v,2)

      allocate (Cape_pert%coin       (isize, jsize) )
      allocate (Cape_pert%plcl       (isize, jsize) )
      allocate (Cape_pert%plfc       (isize, jsize) )
      allocate (Cape_pert%plzb       (isize, jsize) )
      allocate (Cape_pert%xcape      (isize, jsize) )
      allocate (Cape_pert%parcel_r   (isize, jsize, ncap) )
      allocate (Cape_pert%parcel_t       (isize, jsize, ncap) )
      allocate (Cape_pert%cape_p      (isize, jsize, ncap) )
      allocate (Cape_pert%env_r      (isize, jsize, ncap) )
      allocate (Cape_pert%env_t      (isize, jsize, ncap) )
!--------------------------------------------------------------------
!  define 3d execution flag.
!---------------------------------------------------------------------
      do j=1,jsize        
        do i=1,isize        
            if ( (.not. exit_flag(i,j))   .and. &
                 ( .not. nine_flag(i,j)) ) then
               exit_2d(i,j) = .false.
           else
              exit_2d(i,j) = .true.
          endif
         end do
       end do

      do k=1,ncap
        do j=1,jsize    
          do i=1,isize         
            if ( (.not. exit_flag(i,j))   .and. &
                 ( .not. nine_flag(i,j)) ) then
              exit_3d(i,j,k) = .false.
            else
              exit_3d(i,j,k) = .true.
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!  initialize perturbed parcel profile and perturbed parcel 
!  environmental profiles where calculation is to be done.
!--------------------------------------------------------------------
      do k=1,ncap
        do j=1,jsize       
          do i=1,isize      
            if ( .not. exit_3d(i,j,k) ) then

!--------------------------------------------------------------------
!  initialize perturbed parcel profiles to actual parcel profiles.
!---------------------------------------------------------------------
              Cape_pert%parcel_t(i,j,k) = tcape_v(i,j,k)
              Cape_pert%parcel_r(i,j,k) = r_v(i,j,k)

!-------------------------------------------------------------------
!  define environmental profiles for perturbed parcel as the actual 
!  parcel soundings.
!-------------------------------------------------------------------
!             ra_v(i,j,k) = r_v(i,j,k)
              Cape_pert%env_r(i,j,k) = r_v(i,j,k)
!             ta_v(i,j,k) = tcape_v(i,j,k)
              Cape_pert%env_t(i,j,k) = tcape_v(i,j,k)

!--------------------------------------------------------------------
!  if calculation not needed in column, set temperature to ~0.0 to
!  prevent calculation in calculate_cape subroutine.
!--------------------------------------------------------------------
            else
!! CHANGE in connecting to fez
!             ta_v(i,j,k) = 0.007
!             ta_v(i,j,k) = 150.0     
              Cape_pert%env_t(i,j,k) = 150.0     
!             ra_v(i,j,k) = 0.0          
              Cape_pert%env_r(i,j,k) = 0.0          
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!   perturb surface mixing ratio and temperature to calculate 
!   derivative of parcel density temperature w.r.t. surface
!   large-scale density temperature.
!--------------------------------------------------------------------
      do j=1,jsize       
        do i=1,isize      
          if ( (.not. exit_flag(i,j) )    .and. &
               (.not. nine_flag(i,j) ) )  then
!           ra_v(i,j,1) = ra_v(i,j,1) - 0.01*ra_v(i,j,1)
!           ra_v(i,j,1) = max(ra_v(i,j,1), 0.0)
            Cape_pert%env_r(i,j,1) = Cape_pert%env_r(i,j,1) -   &
                                     0.01*Cape_pert%env_r(i,j,1)
            Cape_pert%env_r(i,j,1) = max(Cape_pert%env_r(i,j,1), 0.0)
!            if (ra_v(i,j,1) .lt. 0.) ra_v(i,j,1) = 0.
!           ta_v(i,j,1) = tcape_v(i,j,1) - 1.0
            Cape_pert%env_t(i,j,1) = tcape_v(i,j,1) - 1.0
          endif
        end do
      end do

!--------------------------------------------------------------------
!   call calculate_cape to calculate t and r profiles for perturbed parcel.
!   only Cape_pert%parcel_t and Cape_pert%parcel_r are needed as output from this call.
!--------------------------------------------------------------------

      Cape_pert%cape_p = pcape_v
!     Cape_pert%rcape = ra_v      
!     Cape_pert%tcape = ta_v        
!      Cape_pert%coin =  cin_v2
!      Cape_pert%plcl  = plcl_v2
!     Cape_pert%plfc =  plfc_v2
!     Cape_pert%plzb =  plzb_v2
!     Cape_pert%parcel_r =   rpca_v  
!     Cape_pert%parcel_t = tpca_v
!      Cape_pert%xcape = xcape_v2 

!     call cape_vect (js, Cape_pert, pcape_v, ra_v, ta_v, cin_v2, plcl_v2, plfc_v2,  &
      call calculate_cape (is, ie, js, je,  Cape_pert,   &
                           exit_flag=exit_2d, perturbed=perturbed)

!---------------------------------------------------------------------

      do j=1,jsize   
      do i=1,isize       

       if (.not. exit_3d(i,j,1)) then
         debug_ijt = .false.
         if (in_diagnostics_window ) then
           do n=1,ncols_in_window
             if (j == j_dc(n) .and. i == i_dc(n)) then
               debug_ijt = .true.
               exit
             endif
           end do
         endif

         dcape = dcape_v(i,j)
         plfc = plfc_v(i,j)
         plzb = plzb_v(i,j)
         do k=1,ncap

!          rpca(k) = rpca_v(i,j,k)
           rpca(k) = Cape_pert%parcel_r(i,j,k)
           tpca(k) = Cape_pert%parcel_t(i,j,k)

!          p(k) = pcape_v(i,j,k)
           p(k) = Cape_pert%cape_p(i,j,k)
           qli0(k) = qli0_v(i,j,k)
            qli1(k) = qli1_v(i,j,k)
           qr(k) = qr_v(i,j,k)
           qt(k) = qt_v(i,j,k)
           r(k) = r_v(i,j,k)
            t(k) = tcape_v(i,j,k)
           ri(k) = ri_v(i,j,k)
           rl(k) = rl_v(i,j,k)
           rpc(k) = rpc_v(i,j,k)
           tpc(k) = tpc_v(i,j,k)
         end do
         if ( (.not. exit_flag(i,j) ).and. &
              ( .not. nine_flag(i,j) )) then

           if (debug_ijt) then
             write (unit_dc(n), '(a, 3e20.12)')   &
                        'in cuclo: cpi,cpv,dcape=', cpi, cpv, dcape
             write (unit_dc(n), '(a, 4e20.12)')   &
                       'in cuclo: rocp, rair, latvap,rvap=', rocp, &
                          rair, latvap, rvap
             do k=1,ncap
               write (unit_dc(n), '(a, i4, f19.10, 3e20.12, f20.14)') &
                      'in cuclo: k,p,qr,qt,r,t  =', k, p(k), qr(k), &
                        qt(k), r(k), t(k)
             end do
             do k=1,ncap
               write (unit_dc(n), '(a, i4, f19.10, f20.14, e20.12)') &
                        'in cuclo: k,p,tpc, rpc   =', k, p(k), tpc(k),&
                           rpc(k)            
             end do
             do k=1,ncap
               write (unit_dc(n), '(a, i4, f19.10, 4e20.12)')   &
                       'in cuclo: k,p,qli0,qli1,ri,rl =', k, p(k),   &
                          qli0(k), qli1(k), ri(k), rl(k)
             end do
             write (unit_dc(n), '(a, 2e19.10)')   &
                          'in cuclo: plfc,plzb= ',plfc,plzb
           endif
!---------------------------------------------------------------------
!     calculate total-water mixing ratio. 
!---------------------------------------------------------------------
           do k=1,ncap
             rt(k)=r(k)+ri(k)+rl(k)
           end do



!
!     Calculate density temperatures. No liquid water in cape calculation.
!
           do k=1,ncap
             tden(k)=tpc(k)*(1.+(rpc(k)/rocp  )) 
             tdena(k)=tpca(k)*(1.+(rpca(k)/rocp  ))
           end do
           tdens=t(1)*(1.+(r(1)/rocp  ))
!          tdensa=ta_v(i,j,1)*(1.+(ra_v(i,j,1)/rocp  ))
!          tdensa=ta_v(i,j,1)*(1.+(Cape_pert%rcape(i,j,1)/rocp  ))
           tdensa=Cape_pert%env_t(i,j,1)*  &
                  (1.+(Cape_pert%env_r(i,j,1)/rocp  ))

!
!     Evaluate derivative of parcel density temperature w.r.t. surface
!     large-scale density temperature.
!

            do k=1,ncap
              dtpdta(k)=(tdena(k)-tden(k))/(tdensa-tdens)


            end do
            do k=1,ncap

              if (debug_ijt) then
                write (unit_dc(n), '(a, i4, 2f20.14, e20.12)')  &
                        'in cuclo: k,tden(k),tdena(k),dtpdta(k)= ',   &
                          k,tden(k), tdena(k),dtpdta(k)
              endif
              if (debug_ijt) then
                write (unit_dc(n), '(a, i4, f20.14, e20.12)')  &
                        'in cuclo: k,tpca,rpca= ', k, tpca(k),rpca(k)
              endif

!             write(6,*) 'k,tpca,rpca= ',k,tpca(k),rpca(k)

            end do
!
!     Calculate I1 and I2 integrals from p. 5 of "Cu Closure D" notes.
!
      ri1=0.
      ri2=0.
      k=1
      rild=qt(k)*(rocp  +r(k))/(rocp  *(1.+rt(k)))
      rile=t(k)*(1.+rl(k)+ri(k)-rocp  )*qr(k)
      rile=rile/(rocp  *((1.+rt(k))**2))
      rilf=-t(k)*(rocp  +r(k))*qli0(k)
      rilf=rilf/(rocp  *((1.+rt(k))**2))
      ri2b=t(k)*(rocp  +r(k))/(rocp  *((1.+rt(k))**2))
      ri2b=ri2b*qli1(k)
      sum2=rild+rile+rilf
      sum2=0.


      do k=2,ncap
        if (p(k) .eq. 0.) go to 3
        rilak=-qt(k)*(rocp  +r(k))/(rocp  *(1.+rt(k)))
        rilbk=-t(k)*(1.+rl(k)+ri(k)-rocp  )*qr(k)
        rilbk=rilbk/(rocp  *((1.+rt(k))**2))
        rilck=t(k)*(rocp  +r(k))*qli0(k)
        rilck=rilck/(rocp  *((1.+rt(k))**2))
        rilakm=-qt(k-1)*(rocp  +r(k-1))/(rocp  *(1.+rt(k-1)))
         rilbkm=-t(k-1)*(1.+rl(k-1)+ri(k-1)-rocp  )*qr(k-1)
        rilbkm=rilbkm/(rocp  *((1.+rt(k-1))**2))
        rilckm=t(k-1)*(rocp  +r(k-1))*qli0(k-1)
        rilckm=rilckm/(rocp  *((1.+rt(k-1))**2))
        rila=.5*(rilak+rilakm)
        rilb=.5*(rilbk+rilbkm)
        rilc=.5*(rilck+rilckm)
        ri2ak=t(k)*(rocp  +r(k))/(rocp  *((1.+rt(k))**2))
        ri2ak=ri2ak*qli1(k)
        ri2akm=t(k-1)*(rocp  +r(k-1))/(rocp  *((1.+rt(k-1))**2))
        ri2akm=ri2akm*qli1(k-1)
        ri2a=.5*(ri2ak+ri2akm)
        sum1=rila+rilb+rilc
        ri1=ri1+(alog(p(k-1)/p(k)))*(sum1+dtpdta(k)*sum2)
        ri2=ri2+(alog(p(k-1)/p(k)))*(ri2a-dtpdta(k)*ri2b)
        if (debug_ijt) then
          write(unit_dc(n), '(a, i4, e20.12)')   &
                        'in cuclo: k,dtpdta(k)= ',k,dtpdta(k)
          write (unit_dc(n),   '(a, 3e20.12)')  &
                        'in cuclo: rila,rilb,rilc= ', rila,rilb,rilc
          write (unit_dc(n), '(a, 2e20.12)')  &
                         'in cuclo: ri1,ri2= ',ri1,ri2
          write (unit_dc(n), '(a, 2e20.12)')  &
                         'in cuclo: sum1,sum2= ',sum1,sum2
        endif
      end do

      if (debug_ijt) then
        write (unit_dc(n), '(a, 3e20.12)')  &
                      'in cuclo: rild,rile,rilf= ', rild, rile, rilf
      endif

 3    continue
      if (ri1 .ge. 0) then
        a1_v(i,j) = 0.
        cycle 
      end if
      ri1=rair*ri1
      ri2=rair*ri2
 2    continue
      a1_v(i,j)=-(ri2+dcape)/ri1


     endif


     endif ! exit_3d

     end do
     end do


      deallocate (Cape_pert%coin   )
      deallocate (Cape_pert%plcl   )
      deallocate (Cape_pert%plfc   )
      deallocate (Cape_pert%plzb   )
      deallocate (Cape_pert%xcape  )
      deallocate (Cape_pert%parcel_r    )
      deallocate (Cape_pert%parcel_t   )
      deallocate (Cape_pert%cape_p )
      deallocate (Cape_pert%env_r  )
      deallocate (Cape_pert%env_t  )

end  subroutine cuclo_vect



!#####################################################################

subroutine mulsub_vect (ampt_v, arat_v, plzb_v, pr_v, phalf,  q_v,  &
                        sfcqf_v, sfcsf_v, t_v,  xba_v, xgcm_v, amax_v, &
                        cmui_v, cmus_v, &
                        cual_v,cuq_v, cuql_v, ecds_v, eces_v, emds_v, &
                        emei_v, emes_v, disa_v, disb_v, disc_v,    &
                        disd_v, dise_v, dmeml_v, elt_v, fre_v, qmes_v,&
                        qtmes_v, qtren_v,    &
                        tmes_v, tpre_v, uceml_v, umeml_v, wmms_v, &
                        wmps_v, wtp_v, exit_flag, Don_conv        )

!---------------------------------------------------------------------
real, dimension(:,:,:),intent(in)   :: arat_v
real, dimension(:,:), intent(in)    :: sfcqf_v, sfcsf_v, plzb_v
real, dimension(:,:), intent(out)   :: ampt_v, amax_v,   &
                                       cmui_v, emei_v, tpre_v

real, dimension(:,:,:), intent(in)  ::  pr_v, q_v, t_v
real, dimension(:,:,:), intent(in)  :: phalf
real, dimension(:,:,:), intent(in)  ::  xba_v
real, dimension(:,:,:,:), intent(in)  ::  xgcm_v
type(donner_conv_type), intent(inout) :: Don_conv
real, dimension(:,:,:), intent(out) ::  cmus_v, cual_v, ecds_v, &
                       eces_v, emds_v, emes_v, disa_v, disb_v, disc_v, &
                       disd_v, dise_v, dmeml_v, elt_v, fre_v, qmes_v, &
                      tmes_v, uceml_v, umeml_v, wmms_v, wmps_v, cuq_v, &
                       cuql_v
real, dimension(:,:,:,:), intent(out) :: qtmes_v,qtren_v,wtp_v
logical, dimension(:,:), intent(in) :: exit_flag
!--------------------------------------------------------------------




!
!     Calculates thermal and moisture forcing by an ensemble of cumulus
!     elements, following Donner (1993, JAS). See also LJD notes, "Cu 
!     Closure A," 2/97. Normalized forcings by a_1(p_b).
!     L. Donner  GFDL 27 Apr 97
!
!     On Input:
!
!     arat(kpar)       a_i(p_b)/a_1(p_b) for each of kpar subensembles
!                      Subsensemble kpar should have least entrainment.
!     plzb             level of zero buoyancy (Pa)
!     pr(nlev)         pressure at GCM resolution (Pa)
!                      index 1 at ground.
!     phalf            pressure on half levels (Pa)
!                      index 1 at top of model
!     q(nlev)          mixing ratio at GCM resolution (kg(H2O)/kg)          
!                      index 1 at ground.
!     sfcqf            surface sensible heat flux (W/((m**2))
!                      set to zero if treated by large-scale eddy-fluxes
!     sfcsf            surface moisture flux (kg(H2O)/((m**2) sec))
!                      set to zero if treated by large-scale moisture
!                      fluxes
!     t(nlev)          temperature at GCM resolution (K)
!                      index 1 at ground.
!     xba_v            tracer concentrations at cloud base
!                      (units of xgcm)
!                      (lon,lat,tracer index)
!     xgcm_v           tracer concentrations at GCM resolution
!                      (lon,lat,vert,tracer index)
!                      Vertical index 1 at ground.
!
!     On Input as Parameters:
!
!     kpar             number of cumulus subensembles
!     nlev             number of levels at coarse resolution
!     ncm              number of levels at cloud-model resolution
!
!     On Output:
!     
!     amax             maximum value for a_1(p_b)
!                      See "a Bounds 6/7/97" notes
!     ampt             mesoscale cloud fraction, normalized by a(1,p_b)
!     contot           ratio of convective to total precipitation
!     cmui             normalized vertical integral of mesoscale-updraft
!                      deposition (kg(H2O)/((m**2) sec)
!     cmus(nlev)       normalized mesoscale-updraft deposition
!                      (kg(H2O)/kg/sec)
!     cual(nlev)       cloud fraction, cells+meso, normalized by a(1,p_b)
!     cuq(nlev)        ice content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     cuqll(nlev)      liquid content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     ecds(nlev)       normalized convective downdraft evaporation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     eces(nlev)       normalzed convective-updraft evporation/sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emds(nlev)       normalized mesoscale-downdraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emei             normalized vertical integral of mesoscale-updraft
!                      sublimation (kg(h2O)/((m**2) sec)
!     emes(nlev)       normalized mesoscale-updraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disa(nlev)       normalized thermal forcing, cells+meso (K/sec)
!                      (excludes convergence of surface heat flux)
!                      index 1 at ground. Cumulus thermal forcing defined
!                      as in Fig. 3 of Donner (1993, JAS).
!     disb(nlev)       normalized cell entropy-flux convergence (K/sec)
!                      (excludes convergence of surface flux)
!                      index 1 at ground. Entropy-flux convergence divided
!                      by (p0/p)**(rd/cp).
!     disc(nlev)       normalized cell condensation/deposition
!                      (K/sec)
!                      index 1 at ground.
!     disd(nlev)       normalized cell moisture-flux convergence
!                      (excludes convergence of surface moisture flux)
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dise(nlev)       normalized moisture forcing, cells+meso (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dmeml(nlev)      mass flux in mesoscale downdraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!     elt(nlev)        normalized melting (K/sec)
!                      index 1 at ground.
!     fre(nlev)        normalized freezing (K/sec)
!                      index 1 at ground.
!     pb               pressure at base of cumulus updrafts (Pa)
!     pmd              pressure at top of mesoscale downdraft (Pa)
!     pztm             pressure at top of mesoscale updraft (Pa)
!     qmes(nlev)       normalized mesoscale moisture-flux convergence
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     qtmes(nlev,ncont)  tracer tendency due to mesoscale tracer-flux
!                        convergence (kg/kg/s) (normalized by a(1,p_b))
!                        index 1 at ground 
!     qtren_v          normalized tracer tendency due to cells...
!                      (lon,lat,vert,tracer index)
!                      Vertical index increases as height increases.
!     sfcq(nlev)       boundary-layer mixing-ratio tendency due to surface
!                      moisture flux (kg(H2O)/kg/sec)
!     sfch(nlev)       boundary-layer heating due to surface heat flux
!                      (K/sec)
!     tmes(nlev)       normalized mesoscale entropy-flux convergence
!                      (K/sec)
!                      Entropy-flux convergence is mesoscale component
!                      of second term in expression for cumulus thermal
!                      forcing in Fig. 3 of Donner (1993, JAS).
!                      index 1 at ground.
!     tpre             total normalized precipitation (mm/day)
!     uceml(nlev)      normalized mass fluxes in cell updrafts
!                      (kg/((m**2)*s) 
!     umeml(nlev)      mass flux in mesoscale updraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!                      index 1 at ground.
!     wmms(nlev)       normalized mesoscale deposition of water vapor from
!                      cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wmps(nlev)       normalized mesoscale redistribution of water vapor
!                      from cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wtp_v            tracer redistributed by mesoscale processes
!                      (kg/kg/s) (normalized by a(1,p_b))
!                      vertical index increases with increasing height
!                      (lon,lat,vert,tracer index)
!--------------------------------------------------------------------

      logical, dimension (size(t_v,1), size(t_v,2)) :: exit_mulsub
      real   , dimension (size(t_v,1), size(t_v,2)) :: tb_v, pb_v, &
                                                       qb_v, pt_v, &
                                                       rr_v, ps_v, &
                                                        precip_v, &
                                                       conint_v, &
                                                         dint_v
      real   , dimension (size(t_v,1), size(t_v,2), size(t_v,3) ) :: &
                                  sig_v         
      real   , dimension (size(t_v,1), size(t_v,2), ncap ) :: &
                           tcc_v, wv_v, rcl_v, te_v, qe_v, dpf_v, &
                           dfr_v, flux_v, qlw_v
    real   , dimension (size(t_v,1), size(t_v,2), ncap, ncont ) :: &   
                            xclo_v,xtrae_v
      real, dimension (ncap) :: rlsm, emsm, efchr, emfhr, rlhr, ctfhr, &
                                cmfhr, qllw, rsc, tcc, wv, te, qe, &
                                rcl, cuah, disp, dis, dpf, qlw, dfr, &
                                alp, flux, ucemh, cuql, cuqli, etfhrv
      real, dimension (ncap,ncont)  :: etfhr,etsm,xclo,xtrae
      real, dimension (size(xgcm_v,4))  :: stbl,sumetf
      real, dimension (size(t_v,3))   :: t, q, pr, h1, q1, sig, efc, &
                                         em, rlh, ctf, cmf, fre_vk, &
                                          elt_vk, frea, elta, cual_vk, &
                                         cuml, evap, disf, disg, dish, &
                                         disl, dism, disn, diso, cmu, &
                                         ecd, ece, emd, eme, wmm, wmp, &
                                         thlr, qlr, enctf, encmf, enev,&
                                         fres, elts, tmes_vk, qmes_vk, &
                                         sfcq, sfch, dmeml_vk,   &
                                          uceml_vk, umeml_vk, disg_sv, &
                                         cuq_vk, cuqll_vk, qtrv

      real, dimension (size(t_v,3),ncont) :: xgcm, wtp, qtmes_vk, &
                                             qtren_vk, wtp_vk, qtr
      real, dimension (size(t_v,3)+1) :: phr                   
      real, dimension (kpar) :: cuto, preto, pbma, ptma

    logical :: test2, constab, thetl, lmeso, debug_ijt, in_debug_column
      logical :: lcons
      real    ::         latsub,                       emdi

     integer :: isize, jsize
     integer :: kcont
     real :: contotxx
     integer :: i, j, k, n, kou, ncc, nccm, it
     integer :: ith, itl, jk, ndia
     real :: cappa, dfre, al, aptsum, cutot, catot, ca, apt, conint, &
             dint, cfracl, cfraci, cfract, conpre, cpre, cu, convrat, &
             dints, aal, aalm, cmuxxx, disga, dp, rr, gravit, tmel, &
             tfre, ps, pb, tb, qb, pretot, plzb, pmel, pt, precip, &
             pdeet, esc, qcc, pbmu, ptmu, rc, p, sumehf, sumhlr, &
             sumefc, sumthet, summel, pl, dpp, ph, esh, esl, targ, &
             rh, rl, pit, rgh, rgl, tveh, tvch, dpdzh, ehfh, tvel, &
             tvcl, dpdzl, ehfl, ptt, ehf, tve, tvc, dtv, dpdz, exf, &
             pi, qlhr, emfh, emfl, thetf, emff, emf, sbl, p1, dmela, &
             ssbl, emexxx, pmd, pztm, sumwmp, psmx, esumb, esumc, sumf,&
             summ, sumqme, sumg, sumn, sumelt, sumfre, summes, esum, &
             sumev, esuma, sumemf, sumlhr, es
     real :: etfh,etfl,etf
!ljdtest
     real :: qtrsum,qtmesum
!ljdtest



!---------------------------------------------------------------------
!      call tic ('mulsub', '1')

        isize = size (t_v,1)
        jsize = size (t_v,2)

      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          write(unit_dc(n), '(a, 2i4)') 'in mulsub: i_dc,j_dc= ',&
                                i_dc(n), j_dc(n)
        end do
      endif

!---------------------------------------------------------------------
!   define constants. 
!---------------------------------------------------------------------
       DP=-1000.
      RR=1000.

      GRAVIT=9.80616
      CAPPA=RDGAS   /CP_AIR
      LATSUB=latvap+LATICE
      tmel=273.15
      tfre=258.
      dfre= 10. 

       do j=1,jsize       
       do i=1,isize        

      tpre_v(i,j)=0.
      ampt_v(i,j)   = 0.
      amax_v(i,j)   = 0.
      Don_conv%contot(i,j) = 1.
      cmui_v(i,j)   = 0.
      emei_v(i,j)   = 0.
      Don_conv%pmd_v(i,j)=0.
      Don_conv%pztm_v(i,j)=0.
! add so that value is defined when exit_flag is .true.
      pb_v(i,j) = 0.
!
      DO k=1,nlev
      qtren_v(i,j,k,:)=0.
      qtmes_v(i,j,k,:)=0.
      wtp_v(i,j,k,:)=0.
      cual_v(i,j,k)=0.
      uceml_v(i,j,k)=0.
      umeml_v(i,j,k)=0.
      dmeml_v(i,j,k)=0.
      cmus_v(i,j,k)=0.
      emds_v(i,j,k)=0.
      emes_v(i,j,k)=0.
      wmms_v(i,j,k)=0.
      wmps_v(i,j,k)=0.
      fre_v(i,j,k)=0.
      elt_v(i,j,k)=0.
      tmes_v(i,j,k)=0.
      disd_v(i,j,k)=0.
      disa_v(i,j,k)=0.
      disb_v(i,j,k)=0.
      disc_v(i,j,k)=0.
      dise_v(i,j,k)=0.
      qmes_v(i,j,k)=0.
      ecds_v(i,j,k)=0.
      eces_v(i,j,k)=0.
     end do
     if (.not. exit_flag(i,j)) then
       in_debug_column = .false.
       debug_ijt = .false.
       if (in_diagnostics_window ) then
         do n=1,ncols_in_window
           if ( j == j_dc(n) .and. i == i_dc(n)) then
             in_debug_column = .true.
             debug_ijt = .true.
             exit
           endif
         end do
       endif



       exit_mulsub(i,j) = .false.

! lmeso = .true. for mesoscale, .false. for no mesoscale
!  initialize to .true.
!        lmeso = .true.

       pr(:)   = pr_v(i,j,:)
       q (:)   =  q_v(i,j,:)
       t (:)   =  t_v(i,j,:)
      ps_v(i,j) = phalf(i,j,nlev+1)
      
!---------------------------------------------------------------------
!   define constants. 
!---------------------------------------------------------------------

      PS=ps_v(i,j)
      CONSTAB=.FALSE.
!!! NOTE THIS IS NEEDED TO FIX BUG PRESENT IN LCL
      pb = 0.0
!
      DO 11 k=1,nlev
      SIG(k)=PR(k)/PS
 11   CONTINUE


      if (debug_ijt) then
        do k=1,nlev-kstart_diag+1
          write (unit_dc(n), '(a, i4, f20.14, e20.12, f19.10)')&
                    'in mulsub: k,T,Q,P= ',k,T(k),Q(k),PR(k)
        end do
      endif








      CALL LCL(TB,CONSTAB,PB,QB,PS,SIG,T,Q)

      tb_v(i,j) = tb
      pb_v(i,j) = pb
      qb_v(i,j) = qb

      if (debug_ijt) then
!       print *, 'DONNER_DEEP/mulsub: tb,pb,qb= ',TB,PB,QB
        write (unit_dc(n), '(a, f20.14, f19.10, e20.12)') &
                                'in mulsub: tb,pb,qb= ',TB,PB,QB
        if (constab) then
!         print *, 'DONNER_DEEP/mulsub: constab true'
          write (unit_dc(n), '(a)') 'in mulsub: constab true'
        endif
      endif

      if (constab) then
        exit_mulsub(i,j) = .true.
      endif 
 

     endif
     end do
     end do

!      call toc ('mulsub', '1')
!      call tic ('mulsub', '2')


!     if (mpp_pe() == 3) then
!       print *, 'in section2'
!     endif


       do j=1,jsize         
       do i=1,isize        

       if (.not. exit_flag(i,j)) then

          debug_ijt = .false.
          if (in_diagnostics_window ) then
            do n=1,ncols_in_window
              if( i == i_dc(n) .and. j == j_dc(n)) then
                debug_ijt =.true.
                exit
              endif
            end do
          endif

! lmeso = .true. for mesoscale, .false. for no mesoscale
!  initialize to .true.
          lmeso = .true.

          pr(:)   = pr_v(i,j,:)
          q (:)   =  q_v(i,j,:)
          t (:)   =  t_v(i,j,:)
          ps_v(i,j) = phalf(i,j,nlev+1)
          ps = ps_v(i,j)
      
!---------------------------------------------------------------------
!   define constants. 
!
!---------------------------------------------------------------------

       AL=.5
      THETL=.FALSE.
      TEST2=.false.
      PRETOT=0.
      APTSUM=0.
      CUTOT=0.
      CATOT=0.
      ca=0.
      APT=0.
!
      DO  k=1,nlev
      sfcq(k)=0.
      sfch(k)=0.
      cuml(k)=0.
      DISG(k)=0.
      DISN(k)=0.
      FRES(k)=0.
      frea(k)=0.
      elta(k)=0.
      ELTS(k)=0.
      WMM(k)=0.
      wmp(k)=0.
      ecd(k)=0.
      ece(k)=0.
      cmu(k)=0.
      emd(k)=0.
      eme(k)=0.
      ENEV(k)=0.
      ENCTF(k)=0.
      ENCMF(k)=0.
      wtp(k,:)=0.
      qtmes_vk(k,:)=0.
      qtren_vk(k,:)=0.
      SIG(k)=PR(k)/PS
      end do

      do k=1,nlev+1
        PHR(k)=phalf(i,j,nlev+2-k)
      end do

      DO k=1,ncap
      rcl(k)=0.
      EMSM(k)=0.
      RLSM(k)=0.
      cuah(k)=0.
      cuql(k)=0.
      cuqli(k)=0.
      ucemh(k)=0.
      flux(k)=0.
      etsm(k,:)=0.
      end do

       if ( .not. exit_mulsub(i,j)) then

      do k=1,ncap
        alp(k)=0.
      end do

      tb = tb_v(i,j)
      qb = qb_v(i,j)
      pb = pb_v(i,j)
      plzb=plzb_v(i,j)



      DO 31 KOU=1,KPAR
!     if (mpp_pe() == 3) then
!       print *, 'in      31 loop, kou = ', kou
!     endif
      DO 578 k=1,ncap
      do kcont=1,ncont
        etfhr(k,kcont)=0.
      end do
      EFCHR(k)=0.
      EMFHR(k)=0.
      RLHR(k)=0.
      CTFHR(k)=0.
      CMFHR(k)=0.
      QLLW(k)=0.
 578  RSC(k)=0.


      rr_v(i,j) = rr
      sig_v(i,j,:) = sig(:)

!     if (mpp_pe() == 3) then
!       print *, 'before cloudm'
!     endif
      call cloudm_vect(i, j, tb_v, pb_v, tcc_v, pt_v, wv_v,   &
               rr_v,rcl_v, te_v,qe_v,ps_v,t_v,q_v, sig_v,      cappa, &
               GRAVIT,LATICE,KOU,dfre,TFRE, xba_v,xgcm_v, precip_v,&
               conint_v,dpf_v,  dfr_v,dint_v,flux_v,qlw_v, xclo_v, &
               xtrae_v, debug_ijt, n)

       pmel = pb_v(i,j)
       
       xclo(:,:) = xclo_v(i,j,:,:)
       xtrae(:,:) = xtrae_v(i,j,:,:)
       tcc(:) = tcc_v(i,j,:)
       pt = pt_v(i,j)
       wv(:) = wv_v(i,j,:)
       rcl (:) = rcl_v(i,j,:) 
       te(:) = te_v(i,j,:)
       qe(:) = qe_v(i,j,:)
       precip = precip_v(i,j)
       conint = conint_v(i,j)
       dpf(:) = dpf_v(i,j,:)
       dfr(:) = dfr_v(i,j,:)
       dint = dint_v(i,j)
       flux(:) = flux_v(i,j,:)
       qlw(:) = qlw_v(i,j,:)
!     if (mpp_pe() == 3) then
!       print *, 'after  cloudm'
!     endif
!
!     accumulate normalized convective cloud fraction and
!     mass flux in cell updrafts
!
      do k=1,ncap
         cuah(k)=cuah(k)+arat_v(i,j,kou)*(rcl(k)/rr)**2
         cfracl=1.
         cfraci=0.
         cfract=tfre-dfre
         if (tcc(k) .lt. cfract) then
           cfracl=0.
           cfraci=1.
         end if
         if ((tcc(k) .ge. cfract) .and. (tcc(k) .le. tfre)) then
           cfraci=(tfre-tcc(k))/dfre
           cfracl=(tcc(k)-tfre+dfre)/dfre
         end if
         cuql(k)=cuql(k)+arat_v(i,j,kou)*cfraci*qlw(k)*(rcl(k)/rr)**2
         cuqli(k)=cuqli(k)+arat_v(i,j,kou)*cfracl*qlw(k)*(rcl(k)/rr)**2
         ucemh(k)=ucemh(k)+arat_v(i,j,kou)*flux(k)/(rr**2)
         if ((kou .eq. kpar) .and. (cuah(k) .gt. 0.)) then
           cuql(k)=cuql(k)/cuah(k)
           cuqli(k)=cuqli(k)/cuah(k)
         end if

         if (debug_ijt) then
           write (unit_dc(n), '(a, 2i4, e20.12)')  &
                    'in mulsub: kou,k,ucemh= ',kou,k,ucemh(k)
           do kcont=1,ncont
             write (unit_dc(n), '(a, 3i4, 2e20.12)')  &
                   'in muls:kou,k,kcont,xclo,xtrae=' &
                    ,kou,k,kcont,xclo(k,kcont),xtrae(k,kcont)
           end do
         endif

!        write(6,*) 'kou,k,ucemh= ',kou,k,ucemh(k)

       end do
!
!     If cloud thickness less than pdeep_mc, de-activate mesoscale
!     circulation.
!
      pdeet=pb-pt
      if (pdeet .lt. pdeep_mc) lmeso=.false.
!
!
      call lookup_es(tb, esc)
      qcc=epsilo*esc/(pb-esc)
!
      conint=conint/(rr**2)
      precip=precip/(rr**2)
      dint=dint/(rr**2)

      if (debug_ijt) then
        write (unit_dc(n), '(a, 3e20.12)')   &
                 'in mulsub: conint,precip,dint= ',conint,precip,dint
      endif


      do k=1,ncap

        if (debug_ijt) then
          write(unit_dc(n), '(a, i4, 3e20.12)')  &
                    'in mulsub: k,dfr,dpr,rcl= ',k,dfr(k),dpf(k),rcl(k)
          write (unit_dc(n), '(a, e20.12)')  &
                     'in mulsub: cuah(k)= ',cuah(k)
        endif


        dfr(k)=dfr(k)/(rr**2)
        dpf(k)=dpf(k)/(rr**2)

        if (debug_ijt) then
          write (unit_dc(n), '(a, i4, 3e20.12)')  &
                  'in mulsub: k,dfr,dpr,rcl= ',k,dfr(k),dpf(k),rcl(k)
        endif


      end do
!
      IF (KOU .EQ. 1) THEN
         PBMU=PB
         PTMU=PT
      END IF
!
!
      IF (PB .EQ. PT) THEN

       if (debug_ijt) then
         write (unit_dc(n), '(a, i4, 2f19.10)')  &
                           'in mulsub: kou,pb,pt= ',kou,pb,pt
       endif


       go to 165
      END IF
      IF (PRECIP .EQ. 0.) THEN

        if (debug_ijt) then
          write (unit_dc(n), '(a)') &
                       'in mulsub: PRECIP=0 AFTER CLOUD MODEL'
        endif


        go to 165
      END IF
      CONPRE=CONINT
      CPRE=PRECIP
      CU=CONPRE*86400.
      RC=CPRE*86400.
      PRETOT=PRETOT+RC*arat_v(i,j,kou)
      CUTOT=CUTOT+CU*arat_v(i,j,kou)
!
      if (debug_ijt) then
        write (unit_dc(n), '(a, 2e20.12, a)')  &
                    'in mulsub: CONPRE, CPRE= ', conpre, cpre, & 
                    ' KG/(M**2)/SEC'
        write (unit_dc(n), '(a, e20.12, a)')   &
                     'in mulsub: CONDENSATION PRE= ',CU,' MM/DAY'
        write (unit_dc(n), '(a, e20.12, a)')  &
                     'in mulsub: CLOUD MODEL PRE= ',RC,' MM/DAY'
      endif
 4    CONTINUE
 5    CONTINUE
      DISP(1)=PB
      DO 22 k=1,ncap-1
      P=PB+k*DP
      DISP(k+1)=P
      NCC=k
      IF (P .LT. PT) GO TO 21
  22  CONTINUE
 21   CONTINUE
      NCCM=NCC-1
!
!     Calculate quantities required for realizability check on
!     cloud fraction. See a bounds notes (7/6/97).
!
      do k=1,ncap
        alp(k)=alp(k)+((arat_v(i,j,kou)*(rcl(k)**2))/rcl(1)**2)
      end do
!
!     CALCULATE CUMULUS THERMAL FORCING AND MOISTURE FORCING AT
!     CLOUD-MODEL RESOLUTION
!

      if (debug_ijt) then
         write (unit_dc(n), '(a, 2f19.10)')  'in mulsub: PB,PT= ',PB,PT
      endif

      SUMEMF=0.
      SUMLHR=0.
      SUMEFC=0.
      SUMTHET=0.
      sumetf(:)=0.
      summel=0.
      DO 511 IT=1,ncap-1
      IF (IT.EQ.1) THEN
         PL=PB
         DPP=DP/2.
      END IF
      IF (IT.GT. 1) THEN
         PL=PB+(IT-2)*DP
         DPP=DP
      END IF
      PH=PB+(IT  )*DP
      IF (PH .GE. PT) ITH=IT+1
      IF (PH .LT. PT) THEN
         ITH=IT
         PH=PT
         P=PT
      END IF
!     if (mpp_pe() == 3) then
!print *, 'pl, ph, p, pb, pt = ', pl, ph, p, pb, pt
!     endif
      IF (PL .LE. PT) GO TO 502
      IF (IT.EQ. 1) ITL=1
      IF (IT.NE. 1) ITL =IT-1

      if (debug_ijt) then
        write (unit_dc(n), '(a, 3i4)')  &
                             'in mulsub: IT,ITL,ITH= ',IT,ITL,ITH
        write (unit_dc(n), '(a, 3e20.14)')  &
                         'in mulsub: TCC = ',TCC(IT),TCC(ITL),TCC(ITH)
      endif


!     CALL ESTABL(ESH,TCC(ITH))
!     CALL ESTABL(ESL,TCC(ITL))
!     CALL escomp(TCC(ITH), ESH)
!     CALL escomp(TCC(ITL), ESL)
      CALL lookup_es(TCC(ITH), ESH)
      CALL lookup_es(TCC(ITL), ESL)
      targ=tcc(it)
!     call establ(es,targ)   
!     call escomp(targ, es)   
      call lookup_es(targ, es)   
      rh=epsilo*esh/(ph+(epsilo-1.)*esh)
      rl=epsilo*esl/(pl+(epsilo-1.)*esl)
      pit=pb+(it-1)*dp
      rsc(it)=epsilo*es/(pit+(epsilo-1.)*es)
      RGH=287.05*(1.+.608*RH)
      RGL=287.05*(1.+.608*RL)

      if (debug_ijt) then
        write (unit_dc(n), '(a, 3e20.12)')  &
                              'in mulsub: QE= ',QE(IT),QE(ITL),QE(ITH)
        write (unit_dc(n), '(a, 3e20.12)')   &
                              'in mulsub: TE= ',TE(IT),TE(ITL),TE(ITH)
      endif


      TVEH=TE(ITH)*(1.+.61*QE(ITH))
      TVCH=TCC(ITH)*(1.+.61*RH)
! COMPUTE DP/DZ (DPDZH,DPDZL) AS EQUATION (B3) IN DONNER (1987),
! J. ATM. SCI.,??,PP ???-????.
      DPDZH=-GRAVIT*PH*(AL *TVEH+TVCH)/(RDGAS   *(1.+AL )*TVEH*TVCH)
! COMPUTE VERTICAL EDDY TRANSPORT OF T (EHFH,EHFL) USING EQUATION (3).
      EHFH=((rcl(ith)/rr)**2)*WV(ITH)*DPDZH*(TCC(ITH)-TE(ITH))
      TVEL=TE(ITL)*(1.+.61*QE(ITL))
      TVCL=TCC(ITL)*(1.+.61*RL)
      DPDZL=-GRAVIT*PL*(AL *TVEL+TVCL)/(RDGAS   *(1.+AL )*TVEL*TVCL)
      EHFL=((rcl(itl)/rr)**2)*WV(ITL)*DPDZL*(TCC(ITL)-TE(ITL))
! COMPUTE THE EDDY FLUX CONVERGENCE (EHF) BY FLUX DIFFERENCING
! ACROSS THE LAYER.
      IF (IT .EQ. ITH) THEN
         EFCHR(IT+1)=EHFH/DP
         SUMEFC=SUMEFC+EFCHR(IT+1)*DP/2.
         PTT=PT+DP
         SUMTHET=SUMTHET+EFCHR(IT+1)*( (1.0E05/PTT)**CAPPA )*DP/2.

        if (debug_ijt) then 
!          print *, 'DONNER_DEEP/mulsub: EFCHR(IT+1)= ',EFCHR(IT+1)
           write (unit_dc(n), '(a, e20.12)')  &
                      'in mulsub: EFCHR(IT+1)= ',EFCHR(IT+1)
         endif


         EHFH=0.
      END IF
      EHF=(EHFL-EHFH)/(2.*DPP)
      TVE=TE(IT)*(1.+.609*QE(IT))
      TVC=TCC(IT)*(1.+.609*RSC(IT))
      IF (ITH .NE. IT) P=(PL+PH)/2.
      dtv=tvc-tve
      DPDZ=-GRAVIT*P*(AL *TVE +TVC )/(RDGAS   *(1.+AL )*TVE*TVC)
! COMPUTE THE CUMULUS VERTICAL-FLUX CONVERGENCE OF ENTROPY (EXF)
! AS GIVEN IN UN-NUMBERED EQUATION ON P. 2163.
      EXF=RDGAS   *WV(IT)*DPDZ* (TCC(IT)-TE(IT))*((rcl(it)/rr)**2)/   &
          (CP_AIR   *p )
      EHF=EHF+EXF
      PI= (1.0E05/P)**CAPPA
      EFCHR(IT)       =EHF
      SUMTHET=SUMTHET+EFCHR(IT)* ( (1.0E05/P)**CAPPA )*DPP
      RLHR(IT)=-DPF(IT)
      if (tcc(it) .ge. tfre) then
        convrat=latvap/CP_AIR   
      end if
      if (tcc(it) .lt. tfre) then
        convrat=LATSUB/CP_AIR   
      end if
!     if (mpp_pe() == 3) then
!print *, 'it, ith, tcc(it), tcc(ith), tmel', it, ith, tcc(it), &
!  tcc(ith), tmel
!     endif
      if ((tcc(it) .ge. tmel) .and. (tcc(ith) .le. tmel)) pmel=p
!
!     NO MESOSCALE-NEXT LINE
!
      if (.not. lmeso)    &
            qlw(it)=-dpf(it)*(1.-(rc/cu))
      IF (RLHR(IT) .LT. 0.)  THEN

      if (debug_ijt) then
        write (unit_dc(n), '(a)') 'in mulsub: RLHR .LT. 0.'
      endif


        go to 165
      END IF
      CTFHR(IT)=RLHR(IT)*convrat +EHF
      QLHR  =DPF(IT)
      IF (QLHR   .GT. 0.) QLHR  =0.
! COMPUTE EDDY FLUX CONVERGENCE OF MOISTURE AS IN EQUATION (3).
      EMFH=((rcl(ith)/rr)**2)*WV(ITH)*DPDZH*(RH-QE(ITH))
      EMFL=((rcl(itl)/rr)**2)*WV(ITL)*DPDZL*(RL-QE(ITL))
      IF (IT .EQ. 1) THEN
         THETF=EHFL
         THETF=THETF*( (1.0E05/PL)**CAPPA )
         EMFF=EMFL

         if (debug_ijt) then
           write (unit_dc(n), '(a, 2e20.12)')  &
                             'in mulsub: THETF,EMFF= ',THETF,EMFF
         endif


      END IF
      IF (IT .EQ. ITH) THEN
         EMFHR(IT+1)=EMFH/DP

         if (debug_ijt) then
           write (unit_dc(n), '(a, e20.12)')  &
                          'in mulsub: EMFHRIT=ITH ',EMFHR(IT+1)
         endif


         SUMEMF=SUMEMF+EMFHR(IT+1)*DP/2.
         EMFH=0.
      END IF
      EMF=(EMFL-EMFH)/(2.*DPP)
      SUMEMF=SUMEMF+(EMF*(DPP  )   )
      SUMLHR=SUMLHR+(DPF(IT)*(DPP/GRAVIT) )

      if (debug_ijt) then
        write (unit_dc(n), '(a, i4, e20.12)')  &
                                 'in mulsub: IT,SUMLHR= ',IT,SUMLHR
      endif


      if (tcc(it) .le. tfre) then
         summel=summel+dpp*dpf(it)/gravit
      end if
      SUMEFC=SUMEFC+(EHF*(DPP  )   )

      if (debug_ijt) then
        write (unit_dc(n), '(a, 2e20.12)')  &
                    'in mulsub: SUMTHET,SUMEFC= ',SUMTHET,SUMEFC
      endif


      EMFHR(IT)=EMF

 !     Compute eddy fluxes of tracers.
      do kcont=1,ncont
        etfh = ((rcl(ith)/rr)**2)*wv(ith)*dpdzh*(xclo(ith,kcont)- &
               xtrae(ith,kcont))
        etfl = ((rcl(itl)/rr)**2)*wv(itl)*dpdzl*(xclo(itl,kcont)- &
                xtrae(itl,kcont))
        IF (IT .EQ. ITH) THEN
           etfhr(it+1,kcont) = etfh/dp
           sumetf(kcont) = sumetf(kcont)+etfhr(it+1,kcont)  &
                           *dp/2.
            etfh=0.
        end if
        etf = (etfl-etfh)/(2.*dpp)
        etfhr(it,kcont) = etf
        sumetf(kcont) = sumetf(kcont)+etf*dpp
      end do

      CMFHR(IT)=QLHR  +EMF
      IF (IT .EQ. ITH) THEN
         RLHR(IT+1)=0.
         CMFHR(IT+1)=EMFHR(IT+1)
         CTFHR(IT+1)=EFCHR(IT+1)
      END IF

      if (debug_ijt) then
!      IF (kou .eq. kpar) THEN
       IF (kou .eq. 4) THEN
          write (unit_dc(n), '(a, i4, 2f19.10)') &
                                  'in mulsub: IT,PL,PH= ',IT,PL,PH
          write (unit_dc(n), '(a, 4e20.12)')  &
                       'in mulsub: EHFH,EHFL,EXF,EHF= ',EHFH,EHFL,   &
                                                        EXF,EHF
         write (unit_dc(n), '(a, 3e20.12)')   &
                          'in mulsub: EMFH,EMFL,EMF= ',EMFH,EMFL,EMF
         write (unit_dc(n), '(a, 3e20.12)')   &
                               'in mulsub: ETFH,ETFL,ETF=         ', &
                                                 etfh,etfl,etf
         write (unit_dc(n), '(a, 3e20.12)')   &
                               'etfh diag: rcl,wv,dpdzh=         ', &
                                  rcl(ith),wv(ith),dpdzh
         write (unit_dc(n), '(a, 3e20.12)')   &
                                'etfl diag: rcl,wv,dpdzl=         ', &
                                 rcl(itl),wv(itl),dpdzl
         do kcont=1,ncont
           write (unit_dc(n), '(a, 3e20.12)')  &
                                 'etfh diag: xclo,xtrae= ',         &
                                   xclo(ith,kcont),xtrae(ith,kcont)
           write (unit_dc(n), '(a, 3e20.12)')  &
                                 'etfl diag: xclo,xtrae= ',         &
                                   xclo(itl,kcont),xtrae(itl,kcont)
         end do
         write (unit_dc(n), '(a, 3e20.12)')   &
                            'in mulsub: WV,RH,QE= ',WV(ITH),RH,QE(ITH)
         write (unit_dc(n), '(a, 2e20.12)')  &
                             'in mulsub: RLHR,DPF= ',RLHR(IT),DPF(IT)
         write (unit_dc(n), '(a, 3e20.12)')   &
                             'in mulsub: WV,RL,QE= ',WV(ITL),RL,QE(ITL)
        END IF
      endif


      IF (P .EQ. PT) THEN
         CTFHR(IT)=0.
         CMFHR(IT)=0.
         APT=(rcl(it)/rcl(1))**2
      END IF
 511  CONTINUE
!     if (mpp_pe() == 3) then
!print *, 'after 511 loop  '
!     endif
 502  CONTINUE
!     if (mpp_pe() == 3) then
!print *, 'after 502     '
!     endif

      if (debug_ijt) then
         write (unit_dc(n), '(a, 2e20.12)') &
                            'in mulsub: SUMLHR,SUMEMF= ',SUMLHR,SUMEMF
         write (unit_dc(n), '(a, e20.12)')'in mulsub: SUMTHET=',SUMTHET
         do kcont=1,ncont
           write (unit_dc(n), '(a, e20.12)')'in mulsub: SUMETF=', &
                                          sumetf(kcont)
         end do
      endif


      SBL=0.
!

      if (debug_ijt) then
         write (unit_dc(n), '(a)')'in mulsub: VERAV DFR'
      endif

      lcons=.false.
      CALL VERAV(DFR,PB,PT,lcons, THETL,FREA,SBL,PS,PR,PHR,CAPPA,   &
                 debug_ijt,n)

      if (debug_ijt) then
        do jk=1,nlev
!          print *,   'DONNER_DEEP/mulsub: jk,fres,frea= ',jk,fres(jk),frea(jk)
           write (unit_dc(n), '(a, i4, 2e20.12)')  &
                       'in mulsub: jk,fres,frea= ',jk,fres(jk),frea(jk)
        end do
!      print *,   'summel= ',summel
       write (unit_dc(n), '(a, e20.12)')  'summel= ',summel
      endif

!     if (mpp_pe() == 3) then
!print *, 'before summel '
!     endif

! if (mpp_pe() == 3) then
! print *, 'summel, rc, cu', summel, rc, cu, dint
! endif
      summel=summel*rc/cu
      dints=dint+summel
!
!     No MESOSCALE, Next 6 Lines
!
! if (mpp_pe() == 3) then
! print *, 'lmeso', lmeso
! endif
      if (.not. lmeso) then 
! if (mpp_pe() == 3) then
! print *, 'dint ', dint 
! endif
      if (dint .ne. 0.) then
! if (mpp_pe() == 3) then
! print *, 'pmel,p1, pb', pmel, p1, pb
! endif
      if (pb > pmel) then
        p1=pmel
        dmela=(summel+dint)*gravit/(p1-pb)
        dmela=-dmela*8.64e07

        call ver(dmela,pb,p1,pr,phr,elta, debug_ijt,n)
      endif

      end if
      end if
      summel=-summel*86400.

      if (debug_ijt) then
        write (unit_dc(n), '(a, 3e20.12,a)')   &
                     'in mulsub: summel,rc,cu= ',summel,rc,cu,' mm/day'
        write (unit_dc(n), '(a)')  'in mulsub: VERAV QLW'
      endif

!     if (mpp_pe() == 3) then
!print *, 'before verav - qlw'
!     endif

      lcons = .false.
      CALL VERAV(QLW,PB,PT,lcons, THETL,EVAP,SBL,PS,PR,PHR,CAPPA,  &
                 debug_ijt,n)
      
      if (debug_ijt) then
        write (unit_dc(n), '(a)')  'in mulsub: VERAV RLHR'
      endif

!     if (mpp_pe() == 3) then
!print *, 'before verav - rlhr'
!     endif

      lcons = .false.
      CALL VERAV(RLHR,PB,PT, lcons, THETL,H1,SBL,PS,PR ,PHR ,CAPPA,   &
                  debug_ijt,n)

!
!      Calculate subcloud tracer-flux convergence (kg(tracer)/kg/sec)
!
       do kcont=1,ncont
         stbl(kcont)=0.
         sbl=(stbl(kcont)*gravit)/(ps-pb)
         do it=1,ncap
            etfhrv(it)=etfhr(it,kcont)
         end do
         if (debug_ijt) then
           write (unit_dc(n), '(a)')'in mulsub: VERAV qtrv'
         endif
         lcons = .true.
         call verav(etfhrv,pb,pt,lcons, thetl,qtrv,sbl,ps,pr,phr, &
                    cappa, debug_ijt,n)
         do jk=1,nlev
           qtr(jk,kcont)=qtrv(jk)
         end do
       end do
!     if (mpp_pe() == 2) then
!if (i == 33 ) then
!print *, 'kou, rlhr', kou, (rlhr(k), k=1,100)
!print *, 'kou, h1', kou, (h1(k), k=1,40)
!endif
!     endif
!
!
      SSBL=0.
!
!     CALCULATE SUBCLOUD MOISTURE-FLUX CONVERGENCE (KG(H2O)/KG/SEC)
!
      SBL=(SSBL*GRAVIT       )/(PS-PB)
!
      if (debug_ijt) then
!       print *,   'DONNER_DEEP/mulsub: VERAV EMFHR'
        write (unit_dc(n), '(a)') 'in mulsub: VERAV EMFHR'
      endif

!     if (mpp_pe() == 3) then
!print *, 'before verav - emfhr'
!     endif

      lcons = .true.
      CALL VERAV(EMFHR,PB,PT,lcons,   THETL,Q1,SBL,PS,PR,PHR,CAPPA,  &
                 debug_ijt, n)

!     if (mpp_pe() == 3) then
!print *, 'after  verav - emfhr'
!     endif

      DO 517 JK=1,nlev
         fre_v(i,j,JK)=0.
         EM(JK)=Q1(JK)*8.64E07
         disd_v(i,j,jk)=disd_v(i,j,jk)+em(jk)*arat_v(i,j,kou)
         CMF(JK)=(-H1(JK)+Q1(JK))*8.64E07
         if (t(jk) .ge. tfre) then
           convrat=latvap/CP_AIR
         end if
         if (t(jk) .lt. tfre) then
           convrat=latsub/CP_AIR
         end if
         RLH(JK)=H1(JK)*86400.*convrat
 517  CONTINUE
      THETL=.TRUE.
!
!
      SSBL=0.
!
!     CALCULATE SUBCLOUD ENTROPY-FLUX CONVERGENCE (K/S) DUE ONLY
!     TO SURFACE HEAT FLUX
!
      SBL=GRAVIT*SSBL/((PS-PB)*CP_AIR   )
!

      if (debug_ijt) then
       write (unit_dc(n), '(a)')  'in mulsub: VERAV EFCHR'
      endif

!     if (mpp_pe() == 3) then
!print *, 'before mesub  '
!     endif

      lcons = .true.
      CALL VERAV(EFCHR,PB,PT,lcons, THETL,H1,SBL,PS,PR,PHR,CAPPA, &
                 debug_ijt,n)
      SBL=0.
      THETL=.FALSE.
      PBMA(KOU)=PB
      PTMA(KOU)=PT
      CUTO(KOU)=CU
      PRETO(KOU)=RC

!     if (mpp_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v pre   mesub   ', (elt_v(33,1,k), k=1,40)
!print *, 'elta  pre   mesub   ', (elta (     k), k=1,40)
!endif
!endif

      if (lmeso) then
        elt_vk(:) = elt_v(i,j,:)
        fre_vk(:) = fre_v(i,j,:)
        CALL MESub(CU,RC,DINTS,plzb, PR,PHR,PS,PB,PT,GRAVIT,    &
                   CAPPA, T,CA, tmel,ECD,ECE,fre_vk,elt_vk, debug_ijt,n)
        elt_v(i,j,:) = elt_vk(:)
        fre_v(i,j,:) = fre_vk(:)
        APTSUM=APTSUM+APT*arat_v(i,j,kou)
      else
        ca=0.
        do jk=1,nlev
          ecd(jk)=0.
          ece(jk)=0.
        end do
      end if
      CATOT=CATOT+CA*arat_v(i,j,kou)

!     if (mpp_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v pre   518 loop', (elt_v(33,1,k), k=1,40)
!print *, 'elta  pre   518 loop', (elta (     k), k=1,40)
!endif
!endif


!ljdtest
   qtrsum=0.
!ljdtest

      DO 518 JK=1,nlev
!
!     MESOSCALE, NEXT 2 LINES
!
      if (lmeso) then
      FRES(JK)=FRES(JK)+fre_v(i,j,JK)*arat_v(i,j,kou)

        if (debug_ijt) then
          write (unit_dc(n), '(a, i4, 2e20.12)')  &
                   'in mulsub: jk,fres,fre= ',jk,fres(jk),fre_v(i,j,jk)
          do kcont=1,ncont
            write (unit_dc(n), '(a, 2i4, 2e20.12)')  &
                   'in mulsub: jk,kou,qtr,qtren=         ', &
                   jk,kou,qtr(jk,kcont),qtren_vk(jk,kcont)

!ljdtest
  qtrsum=qtrsum+qtr(jk,kcont)*(phr(jk)-phr(jk+1))
  write (unit_dc(n), '(a, 2i4, e20.12)') 'in mulsub: jk,kou,qtrsum= ', &
   jk,kou,qtrsum
!ljdtest
          end do
        endif

        ELTS(JK)=ELTS(JK)+elt_v(i,j,JK)*arat_v(i,j,kou)
      end if
      FRES(JK)=FRES(JK)+FREA(JK)*arat_v(i,j,kou)

      if (debug_ijt) then
        write (unit_dc(n), '(a, i4, 2e20.12)')  &
                     'in mulsub: jk,fres,frea= ',jk,fres(jk),frea(jk)
      endif


       elts(jk)=elts(jk)+arat_v(i,j,kou)*elta(jk)
       EFC(JK)=H1(JK)*86400.
       CTF(JK)=EFC(JK)+RLH(JK)
       disb_v(i,j,jk)=disb_v(i,j,jk)+efc(jk)*arat_v(i,j,kou)
       disc_v(i,j,jk)=disc_v(i,j,jk)+rlh(jk)*arat_v(i,j,kou)
       DISN(JK)=DISN(JK)+RLH(JK)
       ecds_v(i,j,JK)=ecds_v(i,j,JK)+ECD(JK)*arat_v(i,j,kou)
       eces_v(i,j,JK)=eces_v(i,j,JK)+ECE(JK)*arat_v(i,j,kou)
       ENCTF(JK)=ENCTF(JK)+arat_v(i,j,kou)*CTF(JK)
       ENCMF(JK)=ENCMF(JK)+arat_v(i,j,kou)*CMF(JK)
       ENEV(JK)=ENEV(JK)+arat_v(i,j,kou)*EVAP(JK)
       qtren_vk(jk,:)=qtren_vk(jk,:)+arat_v(i,j,kou)*qtr(jk,:)
       IF (DINT .EQ. 0.) DISG(JK)=DISG(JK)-arat_v(i,j,kou)*((ECD(JK) &
                                  +ECE(JK))*latvap/(CP_AIR   *1000.))
       IF (DINT .NE. 0) THEN
         DISG(JK)=DISG(JK)-arat_v(i,j,kou)*(ECE(JK)*LATSUB/  &
                  (CP_AIR   *1000.))
         DISG(JK)=DISG(JK)-arat_v(i,j,kou)*(ECD(JK)*latvap/  &
                  (CP_AIR   *1000.))
       END IF
 518  CONTINUE
      DO 519 k=1,ncap
         RLSM(k)=(RLHR(k)*arat_v(i,j,kou))+RLSM(k)
         EMSM(k)=(EMFHR(k)*arat_v(i,j,kou))+EMSM(k)
         etsm(k,:)=(etfhr(k,:)*arat_v(i,j,kou))+etsm(k,:)
         DIS(k)=PB+(k-1)*DP
 519  CONTINUE
      NDIA=1

      if (debug_ijt) then
        write (unit_dc(n), '(a)') 'in mulsub: P & QLW'
      endif

!     if (mpp_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v after 518 loop', (elt_v(33,1,k), k=1,40)
!endif
!endif

      DO 716 k=1,ncap-1
        if (debug_ijt) then
          write (unit_dc(n), '(a, i4, 2e20.12)')  &
                        'in mulsub: k, P & QLW', k, DIS(k), QLW(k)
        endif


      IF (WV(k+1) .LE. 0.) GO TO 707
      NDIA=NDIA+1
 716  CONTINUE
 707  CONTINUE


      if (debug_ijt) then
        do k=1,nlev
          write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, p & ctf', k, pr(k), ctf(k)
        end do
        do k=1,nlev
          write (unit_dc(n), '(a, i4, f19.10, e20.12)') &
                       'in mulsub: k, p & cmf', k, pr(k), cmf(k)
        end do
      endif


       if (debug_ijt) then
         write (unit_dc(n), '(a, i4, 2f19.10)')  &
                            'in mulsub: kou,pb,pt= ',kou,pb,pt
       endif

 31   CONTINUE

     
!     if (mpp_pe() == 3) then
!       print *, 'after   31 loop'
!     endif


!
!     Select pt so it refers to the most penetrative sub-ensemble.
!     This is frequently, but not always, the sub-ensemble with
!     the lowest entrainment.
!
      pt=min(ptma(1),ptma(2),ptma(3),ptma(4),ptma(5),ptma(6),ptma(7))
!
!     Convert normalized convective cloud fraction and mass fluxes 
!     to GCM resolution.
!
      sbl=0.

      cual_vk(:) = cual_v(i,j,:)
      lcons = .false.
      call verav(cuah,pb,pt,lcons, thetl,cual_vk,sbl,ps,pr,phr,cappa,  &
                 debug_ijt,n)
      cual_v(i,j,:) = cual_vk(:)
      call verav(cuql,pb,pt,lcons,thetl,cuq_vk,sbl,ps,pr,phr,cappa, &
                 debug_ijt,n)
      cuq_v(i,j,:) = cuq_vk(:)
      call verav(cuqli,pb,pt,lcons,thetl,cuqll_vk,sbl,ps,pr,phr,cappa, &
                 debug_ijt,n)
      cuql_v(i,j,:) = cuqll_vk(:)

      uceml_vk(:) = uceml_v(i,j,:)
      call verav(ucemh,pb,pt,lcons,thetl,uceml_vk,sbl,ps,pr,phr,cappa, &
                 debug_ijt,n)
      uceml_v(i,j,:) = uceml_vk(:)

      apt=aptsum
      ampt_v(i,j)=5.*apt
!
!     Calculate quantities required for realizability check on
!     cloud fraction. See a bounds notes (7/6/97).
!
      al=maxval(alp)
      aal=0.
      aalm=0.
      if (al .gt. 0.) then
        aal=1./al
        aalm=1./(al+ampt_v(i,j))
      end if
      if (lmeso) then
         amax_v(i,j)=aalm
      else
        amax_v(i,j)=aal
      end if

      if (debug_ijt) then
        write (unit_dc(n), '(a, e20.12, a, e20.12)')  &
                         'in mulsub: CUTOT=',CUTOT,' PRETOT=',PRETOT
        write (unit_dc(n), '(a, e20.12)') 'in mulsub: CATOT=',CATOT
      endif


      CATOT=-CATOT
      PB=PBMA(1)

      if (debug_ijt) then     
       write (unit_dc(n), '(a, 3f19.10, l4)')  &
                        'in mulsub: ps,pb,pt,lmeso= ',ps,pb,pt,lmeso
      endif


      if (lmeso) then

        if (debug_ijt) then
         write (unit_dc(n), '(a, e20.12)')  &
                                    'in mulsub: ampt= ',ampt_v(i,j)
        endif


        dmeml_vk(:) = dmeml_v(i,j,:)
        umeml_vk(:) = umeml_v(i,j,:)
        elt_vk(:) = elt_v(i,j,:)
!       qmes_vk(:) = qmes_v(i,j,:)
!       tmes_vk(:) = tmes_v(i,j,:)
        do kcont=1,ncont
!!BUGFIX:
!         xgcm(:,kcont)=xgcm_v(i,j,:,ncont)
          xgcm(:,kcont)=xgcm_v(i,j,:,kcont)
        end do
        CALL MEens(CUTOT,etfhr,PRETOT,PR,PHR,PS,PB,plzb,PT,GRAVIT,    &
                   CAPPA,EPSILO, latvap, T,Q,APT,RLSM, EMSM,   &
                   catot,tmel,xgcm, cuml,CMU,cmuxxx,dmeml_vk,EMD, &
                   emdi,EME,emexxx, pmd, pztm, WMM,WMP, wtp, &
                   elt_vk,contotxx,tmes_vk,qmes_vk,qtmes_vk, &
                   umeml_vk, debug_ijt, n)

        Don_conv%pmd_v(i,j) = pmd
        Don_conv%pztm_v(i,j) = pztm
        cmui_v(i,j) = cmuxxx
        emei_v(i,j) = emexxx
        Don_conv%contot(i,j) = contotxx
        Don_conv%emdi_v(i,j) = emdi
        dmeml_v(i,j,:) = dmeml_vk(:)
        umeml_v(i,j,:) = umeml_vk(:)
        elt_v(i,j,:) = elt_vk(:)
        qmes_v(i,j,:) = qmes_vk(:)
        tmes_v(i,j,:) = tmes_vk(:)
        qtmes_v(i,j,:,:) = qtmes_vk(:,:)
        qtren_v(i,j,:,:) = qtren_vk(:,:)
        wtp_v(i,j,:,:) =wtp(:,:)

       end if
       if (pretot .ne. 0.) then
         if (lmeso) tpre_v(i,j)=pretot/Don_conv%contot(i,j)
         if (.not. lmeso) tpre_v(i,j)=pretot
       end if
       sumwmp=0.

!     if (mpp_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v after meens', (elt_v(33,1,k), k=1,40)
!endif
!endif
!ljdtest
     qtmesum=0.
!ljdtest

         DO 132 JK=1,nlev

           if (debug_ijt) then
             write (unit_dc(n), '(a, i4, f19.10, 2e20.12)') &
                       'in mulsub: jk, pr,cual,cuml= ',jk,pr(jk),  &
                         cual_v(i,j,jk), cuml(jk)
             do kcont=1,ncont
               write (unit_dc(n), '(a, 2i4, f19.10, e20.12)')  &
                    'in mulsub: jk, pr,wtp= ',jk,kcont,pr(jk), &
                                        wtp_v(i,j,jk,kcont)
               write (unit_dc(n), '(a, 2i4, f19.10, e20.12)')  &
                     'in mulsub: jk, pr,qtmes= ',jk,kcont,pr(jk),  &
                                               qtmes_v(i,j,jk,kcont)
               write (unit_dc(n), '(a, 2i4, f19.10, e20.12)')  &
                      'in mulsub: jk, pr,qtren= ',jk,kcont,pr(jk), &
                                                qtren_v(i,j,jk,kcont)

!ljdtest
   qtmesum=qtmesum+qtmes_v(i,j,jk,kcont)*(phr(jk)-phr(jk+1))
    write (unit_dc(n), '(a, i4, e20.12)') 'in mulsub: jk,qtmesum= ', &
      jk,qtmesum
!ljdtest
             end do
           endif


           cual_v(i,j,jk)=cual_v(i,j,jk)+cuml(jk)

           if (debug_ijt) then
             write (unit_dc(n), '(a, i4, 2e20.12)')  &
                          'in mulsub: jk,cuml,cual= ',jk,cuml(jk),   &
                            cual_v(i,j,jk)
             write (unit_dc(n), '(a, i4, 3e20.12)')  &
                     'in mulsub: jk,cmu,emd,eme= ',jk,cmu(jk),emd(jk), &
                                    eme(jk)
             write (unit_dc(n), '(a, i4, 3e20.12)') &
                           'in mulsub: jk,wmm,wmp,elt= ',jk,wmm(jk),  &
                               wmp(jk),  elt_v(i,j,jk)
             write (unit_dc(n), '(a, i4, f20.14, e20.12)')  &
                      'in mulsub: jk,tmes,qmes= ',jk,tmes_v(i,j,jk),  &
                        qmes_v(i,j,jk)
           endif


           cmus_v(i,j,JK)=CMU(JK)-WMM(JK)
           emds_v(i,j,JK)=EMD(JK)
           emes_v(i,j,JK)=EME(JK)
           wmms_v(i,j,JK)=wmms_v(i,j,JK)+WMM(JK)
           wmps_v(i,j,JK)=wmps_v(i,j,JK)+WMP(JK)
           sumwmp=sumwmp+wmps_v(i,j,jk)*(phr(jk)-phr(jk+1))
           if ((phr(jk+1) .le. pb) .and.       &
                (phr(jk) .ge. pb)) psmx=phr(jk+1)
!
!      MESOSCALE, Next Line
!
           if (lmeso) then
             ELTS(JK)=ELTS(JK)+elt_v(i,j,JK)
           end if
 132      CONTINUE
          sumwmp=sumwmp/(gravit*1000.)

          if (debug_ijt) then
            write (unit_dc(n), '(a, e20.12,a)')  &
                      'in mulsub:sumwmp= ',sumwmp,' mm/day'
          endif


 131  CONTINUE
 134  CONTINUE

        if (debug_ijt) then
          write (unit_dc(n), '(a,e20.12, a, e20.12)')  &
              'in mulsub: CATOT= ',CATOT,' contot=',Don_conv%contot(i,j)
         endif


!
!     surface heat flux (W/(m**2))
!
!     if (mpp_pe() == 3) then
!print *, 'ps, psmx', ps, psmx
!     endif
      ssbl=sfcsf_v(i,j)
      sbl=gravit*ssbl/((ps-psmx)*CP_AIR   )
      if (ps > psmx) then
        call ver(sbl,ps,psmx,pr,phr,sfch, debug_ijt, n)
      endif

!
!     surface moisture flux (kg(H2O)/((m**2) sec)
!
      ssbl=sfcqf_v(i,j)
      sbl=(ssbl*gravit)/(ps-psmx)
      if (ps > psmx) then
        call ver(sbl,ps,psmx,pr,phr,sfcq, debug_ijt,n)
      endif
      ESUMB=0.
      ESUMC=0.
      SUMF=0.
      SUMM=0.
      sumqme=0.
      DO 513 JK=1,nlev
      EVAP(JK)=ENEV(JK)*8.64E07
      DISH(JK)=-cmus_v(i,j,JK)
      DISF(JK)=DISH(JK)+ecds_v(i,j,JK)+eces_v(i,j,JK)+wmps_v(i,j,JK)
      disf(jk)=disf(jk)+emes_v(i,j,jk)+emds_v(i,j,jk)
      disf(jk)=disf(jk)+qmes_v(i,j,jk)
      CMF(JK)=ENCMF(JK)
!
      DISM(JK)=CMF(JK)
!
!     MESOSCALE, NEXT LINE
!
      if (lmeso)    &
            CMF(JK)=CMF(JK)+DISF(JK)
!
!     NO MESOSCALE, NEXT 2 LINES
!
      if (.not. lmeso) then
        CMF(JK)=CMF(JK)+EVAP(JK)
        DISF(JK)=EVAP(JK)
      end if
!
      SUMF=SUMF+DISF(JK)*(PHR(JK)-PHR(JK+1))
      sumqme=sumqme+qmes_v(i,j,jk)*(phr(jk)-phr(jk+1))
      SUMM=SUMM+DISM(JK)*(PHR(JK)-PHR(JK+1))
      ESUMB=ESUMB+CMF(JK)*(PHR(JK)-PHR(JK+1))
 513  CONTINUE
!
      ESUMB=ESUMB/(GRAVIT*1000.)
      SUMF=SUMF/(GRAVIT*1000.)
      sumqme=sumqme/(gravit*1000.)
      SUMM=SUMM/(GRAVIT*1000.)

      if (debug_ijt) then
        write (unit_dc(n), '(a, e20.12, a)') &
               'in mulsub: SUMF= ',SUMF, ' MM/DAY'
        write (unit_dc(n), '(a, e20.12, a)') &
              'in mulsub: SUMM= ',SUMM, ' MM/DAY'
        write (unit_dc(n), '(a, e20.12, a)') &
              'in mulsub: ESUMB=',ESUMB, ' MM/DAY'
        write (unit_dc(n), '(a, e20.12, a)') &
               'in mulsub: sumqme= ',sumqme, ' mm/day'
      endif


      SUMG=0.
      SUMN=0.
      sumelt=0.
      sumfre=0.
      summes=0.
      DO 514 JK=1,nlev

        if (debug_ijt) then
!         print *, 'DONNER_DEEP/mulsub: JK,H1= ',JK,H1(JK)
          write (unit_dc(n), '(a, i4, e20.12)')  &
                                 'in mulsub: JK,H1= ',JK,H1(JK)
        endif



      CTF(JK)=ENCTF(JK)
!
      ESUMC=ESUMC+CTF(JK)*(PHR(JK)-PHR(JK+1))
!     DISL(JK)=cmus_v(i,j,JK)*LATSUB/(cpair_mul*1000.)
      DISL(JK)=cmus_v(i,j,JK)*LATSUB/(CP_AIR   *1000.)
      sumelt=sumelt+elts(jk)*(phr(jk)-phr(jk+1))
      sumfre=sumfre+fres(jk)*(phr(jk)-phr(jk+1))
!     fre_v(i,j,JK)=FRES(JK)*LATICE/(cpair_mul*1000.)
      fre_v(i,j,JK)=FRES(JK)*LATICE/(CP_AIR   *1000.)

      if (debug_ijt) then
        write (unit_dc(n), '(a, i4, e20.12)')  &
                            'in mulsub: jk,fres= ', jk, fres(jk)
      endif


!     elt_v(i,j,JK)=-ELTS(JK)*LATICE/(cpair_mul*1000.)
      elt_v(i,j,JK)=-ELTS(JK)*LATICE/(CP_AIR   *1000.)
      DISN(JK)=CTF(JK)+DISG(JK)+fre_v(i,j,JK)+elt_v(i,j,JK)
!     disga=((emes_v(i,j,jk)+emds_v(i,j,jk))*latsub/(cpair_mul    &
      disga=((emes_v(i,j,jk)+emds_v(i,j,jk))*latsub/(CP_AIR       &
             *1000.))

      if (debug_ijt) then
        write (unit_dc(n), '(a, i4, f19.10, 2e20.12)')  &
                    'in mulsub: jk,pr,emds,disga= ',jk,pr(jk),   &
                     emds_v(i,j,jk), disga
      endif


      disg_sv(jk) = disg(jk)
      DISG(JK)=DISL(JK)+DISG(JK)+fre_v(i,j,JK)+elt_v(i,j,jk)-disga+   &
               tmes_v(i,j,jk)
      if (t(jk) .ge. tfre) then
!     DISO(JK)=EVAP(JK)*latvap/(cpair_mul*1000.)
      DISO(JK)=EVAP(JK)*latvap/(CP_AIR   *1000.)
      end if
      if (t(jk) .le. tfre) then
!     DISO(JK)=EVAP(JK)*LATSUB/(cpair_mul*1000.)
      DISO(JK)=EVAP(JK)*LATSUB/(CP_AIR   *1000.)
      end if
!
!     MESOSCALE, NEXT LINE
!
      if (lmeso)     &
            CTF(JK)=CTF(JK)+DISG(JK)
!
!     NO MESOSCALE, NEXT 2 LINES
!
       if (.not. lmeso) then
      CTF(JK)=CTF(JK)-DISO(JK)
      DISG(JK)=-DISO(JK)
      end if
!
      disa_v(i,j,JK)=CTF(JK)
      dise_v(i,j,JK)=CMF(JK)
!
      SUMG=SUMG+DISG(JK)*(PHR(JK)-PHR(JK+1))
      SUMN=SUMN+DISN(JK)*(PHR(JK)-PHR(JK+1))
      summes=summes+tmes_v(i,j,jk)*(phr(jk)-phr(jk+1))
 514  CONTINUE

!     if (mpp_pe() == 2) then
! if (i == 33) then
!print *, 'disa_v in mulsub', (disa_v(33,1,k), k=1,40)
!print *, 'lmeso', lmeso
!print *, 'ctf    in mulsub', (ctf   (     k), k=1,40)
!print *, 'enctf  in mulsub', (enctf (     k), k=1,40)
!print *, 'diso   in mulsub', (diso  (     k), k=1,40)
!print *, 'disg   in mulsub', (disg  (     k), k=1,40)
!print *, 'disg_svin mulsub', (disg_sv(     k), k=1,40)
!print *, 'disl   in mulsub', (disl  (     k), k=1,40)
!print *, 'disga  in mulsub',  disga                  
!print *, 'fre_v  in mulsub', (fre_v (33,1,k), k=1,40)
!print *, 'elt_v  in mulsub', (elt_v (33,1,k), k=1,40)
!print *, 'tmes_v in mulsub', (tmes_v(33,1,k), k=1,40)
!  endif
!      endif
       sumelt=sumelt/(gravit*1000.)
       sumfre=sumfre/(gravit*1000.)
  


!     ESUMC=(ESUMC*cpair_mul)/(GRAVIT*latvap)
      ESUMC=(ESUMC*CP_AIR   )/(GRAVIT*latvap)


!     summes=summes*cpair_mul/(gravit*latvap)
      summes=summes*CP_AIR   /(gravit*latvap)


!     SUMG=SUMG*cpair_mul/(GRAVIT*latvap)
!     SUMN=SUMN*cpair_mul/(GRAVIT*latvap)
      SUMG=SUMG*CP_AIR   /(GRAVIT*latvap)
      SUMN=SUMN*CP_AIR   /(GRAVIT*latvap)


      if (debug_ijt) then
        write (unit_dc(n), '(a, 2e20.12,a)')  &
                     'in mulsub: sumelt,sumfre= ',sumelt,sumfre,  &
                                                ' mm/day'
        write (unit_dc(n), '(a, e20.12)') 'in mulsub: ESUMC=',ESUMC
        write (unit_dc(n), '(a,e20.12,a)')   &
                                'in mulsub: summes=',summes,' mm/day'
        write (unit_dc(n), '(a, e20.12,a)')  &
                                'in mulsub: SUMG=',SUMG,' MM/DAY'
        write (unit_dc(n), '(a,e20.12,a)')  &
                                'in mulsub: SUMN= ',SUMN,' MM/DAY'
      endif

      ESUM=0.
      SUMEV=0.
      ESUMA=0.
      DO 516 k=1,nlev
      ESUMA=ESUMA+CTF(k)*(PHR(k)-PHR(k+1))
      ESUM=ESUM+CMF(k)*(PHR(k)-PHR(k+1))
      SUMEV=SUMEV+EVAP(k)*(PHR(k)-PHR(k+1))
 516  CONTINUE
      ESUM=ESUM/(GRAVIT*1000.)
      SUMEV=SUMEV/(GRAVIT*1000.)
!     ESUMA=(ESUMA*cpair_mul)/(GRAVIT*latvap)
      ESUMA=(ESUMA*CP_AIR   )/(GRAVIT*latvap)

      if (debug_ijt) then
        write (unit_dc(n), '(3(a,e20.12))')  &
                  'in mulsub: ESUM= ',ESUM,' ESUMA=',ESUMA,    &
                  '  SUMEV=',SUMEV
        do k=1,nlev
          if (disc_v(i,j,k) /= 0.00 ) then
            write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                    'in mulsub: k, P & LHR =',  k, PR(k),disc_v(i,j,k)
          endif
       end do
       do k=1,nlev
         if (disb_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)') &
                     'in mulsub: k, P & EFC =',  k, PR(k),disb_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (disd_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & EMF =',  k, PR(k),disd_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (disn(k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & cell thermal forcing =',    &
                        k, PR(k),disn(k)
         endif
       end do
       do k=1,nlev
         if (dism(k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & cell moisture forcing =',    &
                        k, PR(k),dism(k)
         endif
       end do
       do k=1,nlev
         if ( fre_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)') &
                      'in mulsub: k, P & meso up freeze        =',    &
                        k, PR(k),fre_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if ( elt_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                     'in mulsub: k, P & meso down melt        =',    &
                       k, PR(k),elt_v(i,j,k)
         endif
       end do
       do k=1,ncc 
         write (unit_dc(n), '(a, i4, 2e20.12)')  &
                    'in mulsub: k, P & cond/efc              =',    &
                     k, dis(k),ctfhr(k)
         IF (WV(k+1) .LE. 0.) exit      
       end do
       do k=1,nlev
         write (unit_dc(n), '(a, i4, 2e20.12)')  &
                    'in mulsub: k, P & cond/mfc              =',    &
                     k, dis(k),cmfhr(k)
         IF (WV(k+1) .LE. 0.) exit      
       end do
       do k=1,nlev
         if (cmus_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)') &
                      'in mulsub: k, P & meso up con           =',    &
                         k, PR(k),cmus_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (emds_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                     'in mulsub: k, P & meso down evap        =',    &
                       k, PR(k),emds_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (emes_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                        'in mulsub: k, P & meso up evap        =',    &
                          k, PR(k),emes_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (wmms_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & meso cell con       =',    &
                       k, PR(k),wmms_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (wmps_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & meso vap redist     =',    &
                          k, PR(k),wmps_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (tmes_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & meso efc            =',    &
                        k, PR(k),tmes_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (tmes_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)') &
                        'in mulsub: k, P & meso mfc            =',    &
                         k, PR(k),qmes_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (eces_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & up con evap         =',    &
                          k, PR(k),eces_v(i,j,k)
         endif
       end do
       do k=1,nlev
         if (ecds_v(i,j,k) /= 0.00 ) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & down con evap         =',    &
                       k, PR(k),ecds_v(i,j,k)
         endif
       end do
     endif




     if (debug_ijt) then
       do k=1,nlev
         if (ctf(k) /= 0.0) then
           write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                           'in mulsub: k, p & ens thermal forc', &
                              k, pr(k), ctf(k)
          endif
        end do
        do k=1,nlev
          if (cmf(k) /= 0.0) then
            write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                          'in mulsub: k, p & ens moisture forc', &
                           k, pr(k), cmf(k)
          endif
        end do
        do k=1,nlev
          if (disg(k) /= 0.0) then
            write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                         'in mulsub: k, p & thermal modifications', &
                          k, pr(k), disg(k)
          endif
        end do
        do k=1,nlev
          if (disf(k) /= 0.0) then
            write (unit_dc(n), '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, p & moisture modifications', &
                       k, pr(k), disf(k)
          endif
        end do
      endif






!     Convert to MKS units.
!
       do k=1,nlev
         disa_v(i,j,k)=disa_v(i,j,k)/86400.
         disb_v(i,j,k)=disb_v(i,j,k)/86400.
         disc_v(i,j,k)=disc_v(i,j,k)/86400.
         disd_v(i,j,k)=disd_v(i,j,k)/8.64e07
         dise_v(i,j,k)=dise_v(i,j,k)/8.64e07
         fre_v(i,j,k)=fre_v(i,j,k)/86400.
         elt_v(i,j,k)=elt_v(i,j,k)/86400.
         cmus_v(i,j,k)=cmus_v(i,j,k)/8.64e07
         ecds_v(i,j,k)=ecds_v(i,j,k)/8.64e07
         eces_v(i,j,k)=eces_v(i,j,k)/8.64e07
         emds_v(i,j,k)=emds_v(i,j,k)/8.64e07
         emes_v(i,j,k)=emes_v(i,j,k)/8.64e07
         wmms_v(i,j,k)=wmms_v(i,j,k)/8.64e07
         wmps_v(i,j,k)=wmps_v(i,j,k)/8.64e07
         tmes_v(i,j,k)=tmes_v(i,j,k)/86400.
         qmes_v(i,j,k)=qmes_v(i,j,k)/8.64e07
       end do

     endif



165   continue





    endif




  end do
  end do

  Don_conv%pb_v = pb_v

! call toc ('mulsub', '2')

end subroutine mulsub_vect

!####################################################################


subroutine cloudm_vect(i, j, TB_v,PB_v,TCC_v,PT_v,WV_v,RR_v,rcl_v,    &
                       TE_v,QE_v,PRESS_v,T_v,Q_v,SIG_v,       cappa,  &
                       GRAVIT,LATICE,KOU,dfre,TFRE, xba_v,xgcm_v, &
                       PRECIP_v,CONINT_v, DPF_v, DFR_v,DINT_v,flux_v,&
                       QLWA_v, xclo_v, xtrae_v, debug_ijt, n2)

logical, intent(in)                :: debug_ijt
integer, intent(in)                :: n2
integer, intent(in)                :: i, j, kou
 real,    intent(in)           ::                       gravit, &
                                       latice, dfre, tfre
real, dimension(:,:), intent(in)   ::  tb_v, pb_v, rr_v, press_v
real, dimension(:,:), intent(out)  ::  pt_v, precip_v, conint_v,dint_v
real, dimension(:,:,:), intent(in)  :: t_v, q_v, sig_v
real, dimension(:,:,:), intent(in)  :: xba_v
real, dimension(:,:,:,:), intent(in)  :: xgcm_v
real, dimension(:,:,:), intent(out) :: tcc_v, wv_v, rcl_v, te_v, qe_v, &
                                       dpf_v, dfr_v, flux_v, qlwa_v
real, dimension(:,:,:,:), intent(out)  ::  xclo_v,xtrae_v

!
!     ONE-DIMENSIONAL CLOUD MODEL 
!     L. DONNER     NCAR     3 OCT 1984
!     Modified to eliminate dependence of cloud properties on
!     cloud fractional area.
!
!
!     CLOUD MODEL
!
!     ON INPUT:
!
!     TB     cloud temperature at cloud base (K)
!     PB     pressure at cloud base (Pa)
!     rr     cloud radius at base (m)
!     press  surface pressure (Pa)
!     t      large-scale temperature at GCM resolution (K)
!            Index 1 at physical bottom.
!     q      large-scale mixing ratio at GCM resolution (K)
!            Index 1 at physical bottom.
!     sig    (pressure/surface pressure) GCM vertical coordinate
!            Index 1 at physical bottom.
!     rdgas  gas constant for dry air (J/(kg K))
!     epsilo ratio of molecular weights, water vapor to dry air
!     cpair  specific heat for dry air (J/(kg K))
!     gravit gravity constant (m/(s**2))
!     LATICE LATENT HEAT OF FUSION (J/KG)
!     KOU    INDICATOR FOR ENTRAINMENT COEFFICEINT IN ENSEMBLE
!     dfre   freezing range (K)
!     TFRE   freezing temperature (K)
!     xba_v  tracer concentration at cloud base
!            (lon,lat,tracer index)
!     xgcm_v tracer concentration at GCM resolution
!            Vertical index 1 at physical bottom.
!            (lon,lat,vert,tracer index)
!
!     ON OUTPUT:
!     
!     DPF    CONDENSATION RATE*(CLOUD RADIUS**2)  ((M**2)*KG(H2O)/KG/S)
!     precip precipiation integral*(cloud radius**2) (kg/s)
!     CONINT CONDENSATION INTEGRAL*(CLOUD RADIUS**2)  (KG/S)
!     QLWA   CLOUD LIQUID WATER CONTENT (KG(H2O)/KG)
!     DFR    MOISTURE TENDENCY DUE TO FREEZING IN
!            CONVECTIVE UPDRAFT*(CLOUD RADIUS**2) ((M**2)*G(H2O)/(KG DAY))
!     DINT   INTEGRATED WATER MASS FROZEN IN CONVECTIVE
!            UPDRAFT*(CLOUD RADIUS**2) (KG(WATER) /SEC)
!     flux   updraft density*(radius**2)*vert vel (kg/s)
!            Index 1 at physical base of cloud.
!     tcc    cloud temperature (K)
!            Index 1 at physical base of cloud.
!     pt     pressure at cloud top (Pa)
!     wv     cloud vertical velocity (m/s)
!            Index 1 at physical base of cloud.
!     rcl    cloud radius (m)
!            Index 1 at physical base of cloud.
!     te     large-scale temperature at cloud-model resolution (K)
!            Index 1 at physical base of cloud.
!     qe     large-scale mixing ratio at cloud-model resolution (kg/kg)
!            Index 1 at physical base of cloud.
!     xclo_v cloud tracer
!            Vertical index 1 at physical bottom.
!            (lon,lat,vert,tracer index)
!
!     PARAMETER(NCM=100,NCM1=NCM-1)
!     PARAMETER(NLEV=40      )

      real, dimension (ncap) :: sub1, sub2, sub3, sub4, sub5, sub6, &
                                sub7, pf, dpf, qllw, qe, te, dfr, &
                                tcc, wv, rsc, qlwa, disc, disd, disp, &
                                dis, disb, rcl, flux
      real, dimension (ncap, size(xgcm_v,4))  :: xclo,xtrae,clsou
      real, dimension (size(xgcm_v,4)) :: xba
      real, dimension ( size(t_v,3) ) :: sig, t, q
      real, dimension ( size(t_v,3),size(xgcm_v,4) ) :: xgcm
      real  :: xdu
      real  :: wdet  ! vertical velocity at which cloud detrainment 
                     ! begins (m/s)
      real  :: rbound ! value of cumulus radius (m) at which cloud
                      ! model calculation stops
      real  :: wbound ! value of cumulus vertical velocity (m/s) at
                      ! which cloud model calculation stops

      logical :: testlc, test2
      integer :: kcont
      integer :: k, ilvm, n, ih, jk, kk
      real    :: cappa, accond, acpre, actot, sub1sum, pcsave
      real    :: tb, pb, rr, press, pstop,pt, alp, dint, dtfr, fp, dp, &
                 alpp, qcw, qlw, qrw, es, precip, conint, dtupa, dfrac,&
                 sumfrea, p, pp, dtdp, dpd, dwdp, rmu, tcest, west, &
                 rest, estest, qcest, test, qest, qcwt, qlwt, dcw1, &
                 dqrw3, dtdp2, dpd2, dwdp2, dfraca, dtupb, rbar, rmub, &
                 tvb, tvc, dzdp, dz, tav, tavp1, rm, rma, cpm, cpma, &
                 qrwt, &
                 wvb, pdab, pdam, pdap, pi, pip, pim, rmup, pdanp, tv


      save sub1sum, pcsave




  
       
        tb = tb_v(i,j)
        pb  = pb_v(i,j)
        rr  = rr_v(i,j)
        press = press_v(i,j)
         t(:) = t_v(i,j,:)
         q(:) = q_v(i,j,:)
       sig(:) = sig_v(i,j,:)
      xba(:) = xba_v(i,j,:)
      xgcm(:,:) = xgcm_v(i,j,:,:)
  
      TESTLC=.TRUE.
      test2=.false.

      if (debug_ijt) then
          write (unit_dc(n2), '(a, i4)')  'in cloudm: kou= ',kou
      endif

!
!     Set pstop (in Pa) as lowest presure to which cloud can extend.
!
      pstop=4.e03
      wdet = .1
      rbound = .01
      wbound = .01
      IF (PB .le. pstop) PT=pstop
      IF (PB .le. pstop) go to 307

      ALP=.5
      DINT=0. 
      DTFR=0.
      FP=0.
      PF(1)=0.
      DP=-1000.
      if (kou .eq. 1) sub1sum=0.
      if (kou .eq. 1) pcsave=press
!
!     GATE:
!     These values are used for the GATE case in Donner (1993, JAS).
!     The most penetrative sub-ensemble must have KOU=KPAR.
!
      ALPP=.0915
      IF (KOU .EQ. 2) ALPP=ALPP/(1.3)
      IF (KOU .EQ. 3) ALPP=ALPP/(1.8)
      IF (KOU .EQ. 4) ALPP=ALPP/(2.5)
      IF (KOU .EQ. 5) ALPP=ALPP/(3.3)
      IF (KOU .EQ. 6) ALPP=ALPP/(4.5)
!     IF (KOU .EQ. 7) ALPP=ALPP/(40.)
      IF (KOU .EQ. 7) ALPP=ALPP/(10.)
!
!     KEP:
!
!      ALPP=.0915
!      IF (KOU .EQ. 2) ALPP=ALPP/(1.22)
!      IF (KOU .EQ. 3) ALPP=ALPP/(1.56)
!      IF (KOU .EQ. 4) ALPP=ALPP/(2.05)
!      IF (KOU .EQ. 5) ALPP=ALPP/(2.6)
!      IF (KOU .EQ. 6) ALPP=ALPP/(3.21)
!      IF (KOU .EQ. 7) ALPP=ALPP/(7.84)
!      if (kou .eq. 1) alpp=alpp/(5.5)
!
      DO 3 k=1,ncap
      flux(k)=0.
      DPF(k)=0.
      DFR(k)=0.
      DISC(k)=0.
      DISD(k)=0.
      DIS(k)=0.
      QLWA(k)=0.
      SUB1(k)=0.
      SUB2(k)=0.
      SUB3(k)=0.
      SUB4(k)=0.
      SUB5(k)=0.
      SUB6(k)=0.
      SUB7(k)=0.
      PF(k)=0.
      rcl(k)=0.
      DPF(k)=0.
      DISB(k)=PB+(k-1)*DP
      TCC(k)=0.
      QLLW(k)=0.
      QE(k)=0.
      TE(k)=0.
      xclo(k,:)=0.
      xtrae(k,:)=0.
      clsou(k,:)=0.
 3    WV(k)=0.
      if (do_donner_tracer) then
         xclo(1,:)=xba(:)
      endif
      WV(1)=0.5
      rcl(1)=rr
      QCW=0.
      QLW=0.
      QRW=0.
      TCC(1)=TB
      CALL lookup_es(TB, ES)
      RSC(1)=ES*EPSILO/(PB-ES)
      PRECIP=0.
      CONINT=0.
      ILVM=NLEV-1
      N=0

      if (debug_ijt) then
       write (unit_dc(n2), '(a, e20.12, f19.10, f20.14)')  &
                                'in cloudm: RR,PB,TB= ',RR,PB,TB
      endif


      ACCOND=0.
      ACPRE=0.
      dtupa=0.
      dfrac=0.
      sumfrea=0.

      if (debug_ijt) then
      do k=1,nlev-kstart_diag+1
       write (unit_dc(n2), '(a, i4, e20.12)')  &
                      'in cloudm: k,sig   = ',k,sig(k)    
      end do
      endif

      DO 1 k=1,ncap-1
      IF (k .NE. 1) GO TO 10
      CALL lookup_es(TCC(k), ES)
      RSC(k)=EPSILO*ES/(PB-ES)
      IH=0

      if (debug_ijt) then
       write (unit_dc(n2), '(a, i4)')  'in cloudm: ILVM=',ILVM
      endif


      DO 2 JK=1,ILVM


      IF ( (PRESS*SIG(JK) .GE. PB) .AND. (PRESS*SIG(JK+1) .LE. PB) )  &
        IH=JK
 2    CONTINUE

      IF (IH .EQ. 0) THEN
         QE(1)=Q(1)
         TE(1)=T(1)
         if (do_donner_tracer) then
               xtrae(1,:)=xgcm(1,:)
         endif
         GO TO 10
      END IF
      QE(1)=Q(IH)+ ( (Q(IH+1)-Q(IH) )*ALOG(PB/(SIG(IH)*PRESS))/  &
        ALOG(SIG(IH+1)/SIG(IH)) )
      TE(1)=T(IH)+ ( (T(IH+1)-T(IH))*ALOG(PB/(SIG(IH)*PRESS))/   &
        ALOG(SIG(IH+1)/SIG(IH)))
         if (do_donner_tracer) then
      xtraE(1,:)=xgcm(IH,:)+ ( (xgcm(IH+1,:)  &
                 -xgcm(IH,:))*ALOG(PB/(SIG(IH)*PRESS))/  &
                 ALOG(SIG(IH+1)/SIG(IH)))
         endif

       if (debug_ijt) then
         write (unit_dc(n2), '(a, 2e20.12)')  &
                             'in cloudm: QE,TE= ',QE(1),TE(1)
       endif


 10   IH=0

       if (debug_ijt) then
         write (unit_dc(n2), '(a, 2f20.14)')  &
                           'in cloudm: TE,TCC= ',TE(1),TCC(1)
       endif


      P=PB+k*DP
!
      if ( p .lt. pstop) then

       if (debug_ijt) then
         write (unit_dc(n2), '(a, f19.10)')  &
                            'in cloudm: pstop in Clouda= ',pstop
       endif


!      if (p .lt. pstop) test2=.true.
        if (p .lt. pstop) go to 4
      end if
!
      IF ( (P .LE. 50.E03) .AND. (TESTLC) ) THEN
         N=0
         GO TO 4
      END IF
      DO 11 JK=1,ILVM
      IF ( (PRESS*SIG(JK) .GE. P) .AND. (PRESS*SIG(JK+1) .LE. P) )  &
       IH=JK
  11  CONTINUE
      IF (IH .EQ. 0) THEN
         QE(k+1)=Q(1)
         TE(k+1)=T(1)
         if (do_donner_tracer) then
         xtrae(k+1,:) = xgcm(1,:)
         END IF
         GO TO 12
      END IF
      QE(k+1)=Q(IH)+ ( (Q(IH+1)-Q(IH))*ALOG(P/(PRESS*SIG(IH)))/  &
       ALOG(SIG(IH+1)/SIG(IH)))
      TE(k+1)=T(IH)+( (T(IH+1)-T(IH))*ALOG(P/(PRESS*SIG(IH)))/   &
       ALOG(SIG(IH+1)/SIG(IH)))
         if (do_donner_tracer) then
      XtraE(k+1,:)=XGCM(IH,:)+( (XGCM(IH+1,:)  &
                  -XGCM(IH,:))*ALOG(P/(PRESS*SIG(IH)))/ALOG(SIG(IH+1) &
                  /SIG(IH)))
         END IF
 12   CONTINUE
      CALL lookup_es(TCC(k), ES)
      PP=P-DP
      DISP(k)=PP

      if (debug_ijt) then
 
        write (unit_dc(n2), '(a, 3e20.12)')  &
                            'in cloudm: WV,PP,TE= ',WV(k),PP,TE(k)
        write (unit_dc(n2), '(a,f20.14, 2e20.12)')  &
                     'in cloudm: TCC(k),qe(k),QLW= ',TCC(k),WV(k),QLW
      endif





      if   ((pp .le. pcsave) .and. (.not. testlc) .and. &
            (wv(k) .le. wdet)) then
          rmu=0.
      else
          rmu=2.*alpp/rcl(k)
      end if ! ((pp .le. pcsave) .and. (.not. testlc) .and.
             ! (wv(k) .le. wdet))

      CALL SIMULT(TCC(k),rcl(k),WV(k),PP,TE(k),QE(k),QLW,    &
                  TESTLC,pcsave,DTDP,DPD ,DWDP, LATICE,TFRE  &
                  , cappa, GRAVIT, RMU, wdet, debug_ijt, n2)

      if (debug_ijt) then
       write (unit_dc(n2), '(a, 3e20.12)')   &
                               'in cloudm: QE,QLW,RR= ',QE(k),QLW,RR
        write (unit_dc(n2), '(a, 2e20.12)')  &
                                'in cloudm: DPD,DWDP= ',DPD,DWDP
      endif


      TCEST=TCC(k)+DTDP*DP
      WEST=WV(k)+DWDP*DP
      rEST=rcl(k)+DPD*DP
      rbar=(rest+rcl(k))/2.

       if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')  &
                         'in cloudm: rest,west= ',rest,west
       endif


      IF (rbar .LE. rbound) GO TO 4
      IF (WEST .LE. wbound) GO TO 4
      if (rest .le. rbound) go to 4
      CALL lookup_es(TCEST, ESTEST)
      QCEST=EPSILO*ESTEST/(P-ESTEST)
      TEST=TE(k+1)
      QEST=QE(k+1)
      QRWT=QRW
      QCWT=QCW
      QLWT=QLW
      CALL MICRO(TCC(k),TCEST,PP,P,TE(k),TEST,QE(k),QEST,WV(k),   &
                 WEST,rbar,RMU,QRWT,QCWT,QLWT,DCW1,DQRW3,   &
                 debug_ijt, n2)

       if (debug_ijt) then
         write (unit_dc(n2), '(a, f19.10, f20.14, e20.12)') &
                              'in cloudm: P,TEST,QEST= ',P,TEST,QEST
         write (unit_dc(n2), '(a, f20.14, 2e20.12)')   &
                       'in cloudm: TCEST,rest,QLWT= ',TCEST,REST,QLWT
       endif

       if ((p .le. pcsave) .and. (.not. testlc) .and. &
           (west .le. wdet)) then
          rmu=0.
       else
         rmu=2.*alpp/rest
       end if ! ((p .le. pcsave) .and. (.not. testlc) .and.
              !  (west .le. wdet))

      CALL SIMULT(TCEST,rEST,WEST,P,TEST,QEST,QLWT,TESTLC,    &
                  pcsave,DTDP2,DPD2,DWDP2, LATICE,TFRE  &
                  , cappa, GRAVIT, RMU, wdet, debug_ijt, n2)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, l4, 2e20.12)')  &
                   'in cloudm: TESTLC,DTDP2,DPD2= ',TESTLC,DTDP2,DPD2
        write (unit_dc(n2), '(a, e20.12)')  'in cloudm: DWDP2= ',DWDP2
      endif


      rcl(k+1)=rcl(k)+DPD2*DP
      WV(k+1)=WV(k)+DWDP2*DP
      IF ((WV(k+1) .LE. wbound) .OR. (rcl(k+1) .LE. rbound))  THEN
         rcl(k+1) = 0.0
         WV(k+1)=0.
         GO TO 4
      END IF
      DTDP=(DTDP+DTDP2)/2.
      DWDP=(DWDP2+DWDP)/2.
      DPD=(DPD+DPD2)/2.
      TCC(k+1)=TCC(k)+DTDP*DP
!
!     ADD EFFECT OF FREEZING
!
      IF ((TCC(k) .GE. TFRE) .AND. (TCC(k+1) .LE. TFRE) .AND.    &
        (DTFR .EQ. 0.)) THEN
         DTFR=QLWT*LATICE/CP_AIR
!
!     Multiply by factor to take account that not all this
!     water will freeze before falling out. Use Leary and Houze
!     (JAS,1980) rc/cu ratio to estimate this ratio.
      dtfr=.52*dtfr
!
      end if
      dfraca=(tfre-tcc(k+1))/dfre
      if (dfraca .gt. 1.) dfraca=1.
      dfrac=amax1(dfrac,dfraca)
      dtupb=dtfr*dfrac 
      dfr(k+1)=dtupb-dtupa
      sumfrea=sumfrea+dfr(k+1)
      dtupa=dtupb
      TCC(k+1)=TCC(k+1)+dfr(k+1)

      if (debug_ijt) then
        write (unit_dc(n2), '(4(a, e20.12),a, i4)') &
             'in cloudm: DTFR=',DTFR,' dfr=',dfr(k+1),   &
               'LATICE=',LATICE, 'CPAIR=',CP_AIR,'k= ',k
      endif


      WV(k+1)=WV(k)+DWDP*DP
      rcl(k+1)=rcl(k)+DPD*DP
      IF ((WV(k+1) .LE. wbound) .OR. (rcl(k+1) .LE. rbound))  THEN
         rcl(k+1) = 0.0
         WV(k+1)=0.
         GO TO 4
      end if
!
!     Ensure that all ensemble members have pt .le. pressure at which
!     ensemble member 1 becomes buoyant. Simult subroutine will use
!     value of pcsave to do so.
!
      if (kou .eq. 1) then
        if ((testlc) .and. (WV(k+1) .GT. WV(k))) pcsave=p
        IF (WV(k+1) .GT. WV(k)) TESTLC=.FALSE.
      else
        IF (WV(k+1) .GT. WV(k)) TESTLC=.FALSE.
      end if
      CALL lookup_es(TCC(k+1), ES)
      RSC(k+1)=EPSILO*ES/(P-ES)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, f20.14, 2e20.12)')  &
                           'in cloudm: TE,QE,RSC= ',TEST,QEST,RSC(k+1)
      endif

 
      RBAR=(rcl(k)+rcl(k+1))/2.
      if (rbar .le. rbound) go to 4
      RMUB=2.*ALPP/RBAR
      TVB=TE(k)*(1.+.61*QE(k))
      TVB=TVB+TEST*(1.+.61*QEST)
      TVB=TVB/2.
      TVC=TCC(k)*(1.+.61*RSC(k))
      TVC=TVC+TCC(k+1)*(1.+.61*RSC(k+1))
      TVC=TVC/2.
      DZDP=-(1.+ALP)*TVC*TVB*RDGAS   
      DZDP=2.*DZDP/((P+PP)*GRAVIT*(ALP*TVB+TVC))
      sub1(k)=wv(k)/dzdp
      sub1(k+1)=wv(k+1)/dzdp
      DZ=DP*DZDP
!
!     CALL MICROPHYSICS ROUTINE.
!
      CALL MICRO(TCC(k),TCC(k+1),PP,P,TE(k),TE(k+1),QE(k),QE(k+1),  &
                 WV(k),WV(k+1),RBAR,RMU,QRW,QCW,QLW,DCW1,DQRW3, &
                 debug_ijt, n2)
      DISC(k+1)=QCW
      DISD(k+1)=QRW
      TAV=TCC(k)*(1.+.61*RSC(k))
      TAVP1=TCC(k+1)*(1.+.61*RSC(k+1))
      QLWA(k+1)=QLW

      if (debug_ijt) then
        write(unit_dc(n2), '(a, f20.14, 2e20.12)') &
                     'in cloudm: TCC(k+1),WV(k+1),QLW= ',TCC(k+1),  &
                      WV(k+1),QLW
        write (unit_dc(n2), '(a, 2e20.12)')  &
                        'in cloudm: DCW1,DQRW3= ',DCW1,DQRW3
      endif


      RM=RDGAS   *(1.+.609*RSC(k))
      RMA=RDGAS   *(1.+.609*RSC(k+1))
      CPM=CP_AIR*(1.+.87*RSC(k))
      CPMA=CP_AIR*(1.+.87*RSC(k+1))
      WVB=(WV(k)+WV(k+1))/2.
!  
!    Calculate in-cloud tracer distribution.
!
       if (do_donner_tracer) then
      do kcont=1,ncont
        xdu=xclo(k,kcont)
        call Clotr(alp,rmub,gravit,pp,p,rdgas,clsou(k,kcont),  &
                   clsou(k+1,kcont),tcc(k),tcc(k+1),te(k),te(k+1), &
                   xtrae(k,kcont), xtrae(k+1,kcont),xdu, debug_ijt, n2)
        xclo(k+1,kcont) = xdu
      end do
       endif
!
!
      pdab=rbar**2
      pdam=(rcl(k))**2
      pdap=(rcl(k+1))**2
!
      DIS(k)=DQRW3*PDAB*.001*WVB*GRAVIT/DP
      PF(k)=.001*DCW1*WVB*PDAB*GRAVIT/DP
!
!     CALCULATE MOISTURE SUBJECT TO FREEZING
!
      IF (  (DFR(k+1) .NE. 0.) ) THEN
        DFR(k+1)=DFR(k+1)*P*WV(k+1)*pdap*GRAVIT/(RDGAS   *TAVP1*DP)
        DFR(k+1)=-DFR(k+1)*CP_AIR/LATICE
        DINT=DINT-DFR(k+1)*DP/GRAVIT
        DFR(k+1)=DFR(k+1)*8.64E07
      END IF
!
      ACCOND=ACCOND+DCW1
      ACPRE=ACPRE+DQRW3
      ACTOT=ACCOND+ACPRE

      if (debug_ijt) then
        write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                          'in cloudm: k,ACCOND,ACPRE= ',k,ACCOND,ACPRE
        write (unit_dc(n2), '(a, i4, e20.12)')  &
                           'in cloudm:  k,ACTOT= ',k,ACTOT
      endif


      IF (k .EQ. 1) FLUX(k)=(rcl(k)**2)*WV(k)*PP/(RM*TCC(k))
      FLUX(k+1)=(rcl(k+1)**2)*WV(k+1)*P /(RMA*TCC(k+1))

      if (debug_ijt) then
        write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                     'in cloudm: k,rcl,RMU= ',k,rcl(k),RMU
      endif


      IF (QLW .lt. 0.) GO TO 4
      PI=(1.0E05/PP)**(RDGAS   /CP_AIR)
      PIP=(1.0E05/P)**(RDGAS   /CP_AIR)
      PIM=PP-DP
      PIM=(1.0E05/PIM)**(RDGAS   /CP_AIR)
      IF (k .EQ. 1) THEN
         SUB1(k)=pdam*sub1(k)*(RSC(k+1)-RSC(k))/DP
         SUB3(k)=-(PIP*TCC(k+1)-PI*TCC(k))/DP
      END IF
      IF (k .NE. 1) THEN
         SUB1(k)=pdam*sub1(k)*(RSC(k+1)-RSC(k-1))/(2.*DP)
         SUB3(k)=-(PIP*TCC(k+1)-PIM*TCC(k-1))/(2.*DP)
      END IF
      sub1(k)=sub1(k)-pf(k)
      RMUP=RMUB*DZDP
      SUB2(k)=RMUP*PI*latvap*(QE(k)-RSC(k))/CP_AIR
      SUB4(k)=RMUP*PI*(TE(k)-TCC(k))
      N=k
 1    CONTINUE
 4    CONTINUE

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')  &
                  'in cloudm: DINT IN CLOUDM,sumfrea= ',DINT,sumfrea
      endif


      KK=N
      IF (KK .EQ. 0.) GO TO 23

      if (debug_ijt) then
        write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                            'in cloudm: n,rsc(n+2),rsc(n)= ',n,   &
                             rsc(n+1),rsc(n)
      endif

 
      pdanp=rcl(n+1)**2
      SUB1(n+1)=pdanp*sub1(n+1)*(RSC(n+1)-RSC(n))/(2.*DP)
 24   FORMAT(2X,'PRECIP=',E15.5,'KG/(M**2)*SEC')
      DO 22 K=1,KK
        P=PB+(K-1)*DP
        sub1sum=sub1sum+sub1(k)*dp
        IF (K .EQ. 1) THEN
          DPF(1)=PF(1)
        ELSE
          DPF(K)=(PF(K)+PF(K-1))/2.
        END IF
        IF (K .EQ. KK) DPF(KK+1)=PF(KK)/2.

        if (debug_ijt) then
          write (unit_dc(n2), '(a, i4, 2e20.12)')   &
                           'in cloudm: K,PF,DIS= ',K,PF(K),DIS(K)
        endif


       PI=(1.0E05/P)**(RDGAS   /CP_AIR)
       TV=TCC(K)*(1.+.61*RSC(K))
       SUB5(K)=-PI*latvap*DPF(K)*RDGAS   *TV/(CP_AIR*   &
       WV(K)*P*GRAVIT*pdam)
       SUB6(K)=SUB1(K)+SUB2(K)
       SUB7(K)=SUB3(K)+SUB4(K)
       SUB6(K)=-SUB6(K)*WV(K)*P*GRAVIT*PDAM/(RDGAS   *TV)
       SUB7(K)=-SUB7(K)*WV(K)*P*GRAVIT*PDAM/(RDGAS   *TV)

       if (debug_ijt) then
         write (unit_dc(n2), '(a, i4, f19.10, e20.12)')  &
                                'in cloudm: K,P,SUB1= ',K,P,SUB1(K)
       endif


       if (k .eq. n) then

         if (debug_ijt) then
           write (unit_dc(n2), '(a, i4, f19.10, e20.12)') &
                            'in cloudm: N,P,SUB1+= ',n,P,SUB1(n+1)
         endif


       end if
       DIS(K)=-DIS(K)
       CONINT=CONINT+PF(K)*DP/GRAVIT
 22   CONTINUE
      if (kou .eq. 7) sub1sum=86400.*sub1sum/gravit

      if (debug_ijt) then
        write (unit_dc(n2), '(a, i4, a, e20.12,a)')   &
                   'in cloudm: kou= ',kou,' Vapr Adv= ',sub1sum,  &
                   ' mm/day'
        write (unit_dc(n2), '(a, e20.12, a)')  &
                   'in cloudm: CONINT= ',CONINT,' KG/(M**2)/SEC'
      endif


      PRECIP=CONINT*ACPRE/ACCOND
      SUB6(1)=0.
      SUB7(1)=0.
 23   PT=PB+KK*DP

       if (debug_ijt) then
         write (unit_dc(n2), '(a, i4, f19.10)') &
                        'in cloudm: kou,pcsave= ',kou,pcsave
       endif





        dpf_v(i,j,:) = dpf(:)
        precip_v(i,j) = precip
        conint_v(i,j) = conint
        dfr_v(i,j,:) = dfr(:)
        dint_v(i,j) = dint
        tcc_v(i,j,:) = tcc(:)
        wv_v(i,j,:) = wv(:)
        rcl_v(i,j,:) = rcl(:)
        te_v(i,j,:) = te(:)
        qe_v(i,j,:) = qe(:)
        flux_v(i,j,:) = flux(:)
        qlwa_v(i,j,:) = qlwa(:)
        xtrae_v(i,j,:,:) = xtrae(:,:)
        xclo_v(i,j,:,:) = xclo(:,:)

    307  continue

        pt_v(i,j) = pt


end subroutine cloudm_vect

!#####################################################################

subroutine lcl(tb,constab,pb,qb,press,sig,t,qin)

real, intent(inout)   ::  pb
real, intent(out)   :: tb, qb
logical, intent(inout) :: constab
real, dimension(:), intent(in) :: t, sig, qin
real, intent(in)     :: press

     integer :: istart, i
     real :: dp, q, tc, p, pt, gam, es, rs


      istart=2
      dp=-100.
      q=qin(istart)
      tc=t(istart)
      do i=1,600
      p=press*sig(istart)+dp*(i-1.)
      pt=press*sig(istart)+dp*i
      GAM=.286*(1.-.26*Q)*TC/(P/      press       )
      TC=TC+GAM*DP/    press
!     CALL ESTABL(ES,TC)
!     CALL escomp(TC, ES)
      CALL lookup_es(TC, ES)
      RS=.622*ES/(PT-ES)
      IF (RS .LE. Q) PB=PT
      IF (RS .LE. Q) TB=TC
      IF (RS .LE. Q) QB=qin(istart)
      IF (RS .LE. Q) GO TO 6
      end do
      IF (PB .EQ. 0.) CONSTAB=.TRUE.
 6    CONTINUE



end subroutine lcl


!#####################################################################


subroutine satad(t,p,lat,dtp)
!!  IT APPEARS THAT THIS ROUTINE IS NOT CURRENTLY ACCESSED BY THE MODEL

real, intent(in) :: t, p,lat
real, intent(out) :: dtp


!
!      computes saturated adiabatic lapse rate
!
!      on input:
!        t    temperature (K)
!        p    pressure (Pa)
!        lat  latent heat (J/kg)
!      on output:
!        dtp  saturated adiabatic lapse rate (K/m)
!
!
!     constants
!
      real :: rair_satad, gravit, cpair_satad, epsilo_satad, &
              rh2o_satad, es, qs, tv, desdt, dtpd
      rair_satad=287.04
      gravit=9.80616
       cpair_satad=1004.64
      epsilo_satad=.622
      rh2o_satad=461.
!     write(6,*) 'satad t,p= ',t,p
!
!     calculate saturation specific humidity
!
!     call establ(es,t    )
!     call escomp(t, es    )
      call lookup_es(t, es    )
      qs=epsilo_satad*es/(p-es)
      tv=t*(1.+.61*qs)
!
!     calculate saturated adiabatic lapse rate
!
      dtp=tv+(lat*qs/rair_satad)
      dtp=-dtp*gravit/(cpair_satad*t)
      desdt=lat*es/(rh2o_satad*(t**2))
      dtpd=1.+(epsilo_satad*lat*desdt/(cpair_satad*p))
      dtp=dtp/dtpd
!     return

end subroutine satad





!####################################################################


subroutine simult(tc,rda,wv,p,te,qe,qlw,testlc,pcsave,   &
                  dtdp,drdp,dwdp, latice,tfre, cappa, gravit, &
                  rmu, wdet, debug_ijt, n2)

!--------------------------------------------------------------------
logical, intent(in) :: testlc, debug_ijt
integer, intent(in) :: n2
real, intent(in) ::  tc, rda, wv, p, te, qe, qlw, pcsave, &
                       cappa,   latice, tfre, gravit, &
                      wdet,rmu
real, intent(out) :: dtdp, drdp, dwdp
!--------------------------------------------------------------------
!
!     Generates cloud profiles.
!     See LJD "Cloud Model 89" notes on Generalized mu (10/1/89)
!     and dwdz (10/2/89). The value of epsilon is taken as 1.
!     Version where cloud propeties independent of cloud area.
!
!     On input:
!
!        tc   cu temperature (K) at pressure p (Pa)
!        rda  cu radius at pressure p
!        rmu  entrainment coefficient (/m)
!        te   environmental temperature (K) at pressure p
!        qe   environmental mixing ratio at pressure p
!        wv   cu vertical velocity (m/s) at pressure p
!        qlw  liquid water
!        testlc  indicator(logical)
!        pcsave  pressure at which cloud ensemble 1 becomes
!                buoyant (Pa)
!        wdet    vertical speed at which cumulus detrainment begins 
!                (m/s)
!
!     On Output:
!
!        dtdp    temperature derivative (K/Pa)
!        drdp    cu radius derivative (m/Pa)
!        dwdp    cu vertical velocity derivative (m/s/Pa)
!
!



      real lat
      real  :: alp, c2, c3, dp, epm, es, rsc, htve, htv, da, der, &
               fm, rmup, fmp, rst, teae, tcae, dz, c4, test, c5, dwdz, &
               wtest, west, pdap, dadp



!
!      assign constants
!
      alp=.5
      dp=-1000.
      epm=0.
      lat=latvap
      if (tc .lt. tfre) lat=latvap+latice
!
      call lookup_es(tc, es)
      rsc=epsilo*es/(p+(epsilo-1.)*es)
       htve=te   *(1.+.609*qe   )
      HTV=TC    *(1.+.609*RSC   )

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')  &
                         'in simult: htve,htv= ',htve,htv
      endif

      DA=RDGAS   *(1.+ALP)*htv*htve    /(GRAVIT*p * (ALP*htve +htv   ))
 102  FORMAT(2X,3E15.5,I5,2E15.5)
      DER=latvap*ES/(RH2O*(TC **2   ))
      fm=p*(rda**2)*wv/(RDGAS   *htv)
      rmup=-da*rmu
      fmp=fm*exp(rmup*dp)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, f19.10, 2e20.12)')  &
                           'in simult: p,fm,fmp= ',p,fm,fmp
      endif


      rst=RDGAS   *(1.+.608*rsc)
      c2=(htv+(lat*rsc/rst))*gravit/(CP_AIR*tc)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, e20.12)')  'in simult: c2= ',c2
!       print *,  'DONNER_DEEP/simult: te,p,qe= ',te,p,qe
        write (unit_dc(n2), '(a, f20.14, f19.10, e20.12)')  &
                                     'in simult: te,p,qe= ',te,p,qe
      endif


      call tae(te,p,qe,lat,       dp,             cappa,teae)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, f20.14)') 'in simult: teae= ',teae
      endif


      tcae=tc*exp(lat*rsc/(CP_AIR*tc))
      dz=-dp*da
      c3=(tcae-teae)*(1.-exp(rmu*dz))/(exp(rmu*dz)*dz   &
         *exp(lat*rsc/(CP_AIR*tc)))

      c4=1.+(epsilo*lat*der/(CP_AIR*p))

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')  'in simult: c3,c4= ',c3,c4
      endif


      dtdp=(c2-c3)*da/c4
      test=tc+dtdp*dp
      c5=gravit*(htv-htve)/(htve*(1.+alp))
!
!     No  Liq  Remove Next Line
!
      c5=c5-gravit*qlw
      dwdz=2.*c5
      wtest=(wv**2)+dwdz*dz
      if (wtest .lt. 0.) wtest=0.
      wtest=sqrt(wtest)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, l4, e20.12)')  &
                         'in simult: testlc,wtest= ',testlc,wtest
      endif


        dwdp=(wtest-wv)/dp
!
!     Allow for entrainment to change vertical velocity if
!     parcel has (1) achieved initial acceleration and pressure
!     is less than (or equal to) that where most entraining parcel init        ially
!     accelerates or (2) parcel is accelerating. The first 
!     condition ensures that a shallow, lower-entrainment cloud
!     will not develop below the pressure where the most entraining
!     parcel initially accelerates.

      if ( ((.not. testlc) .and. (p .le. pcsave)) .or.   &
          (dwdp .le. 0.) ) then


      dwdz=wv*(1.-exp(rmu*dz))
      dwdz=dwdz/(dz*exp(rmu*dz))
      dwdp=(-dwdz*da)+dwdp
      end if
      
      if (debug_ijt) then
        write (unit_dc(n2), '(a,l4,  2e20.12)')  &
                     'in simult: testlc,dwdp,test= ',testlc,dwdp,test
      endif


!
!     If parcel has not initially accelerated, do not allow
!     vertical velocity to fall. BL turbulence or other
!     sub-grid mechanisms assumed to maintain upward motion.
!
      if ( (testlc) .and. (dwdp .gt. 0.) ) dwdp=0.
!
!     If parcel between bottom surface and pressure at which
!     most entraining parcel initially accelerates, do not
!     allow vertical velocity to fall. This ensures that 
!     less entraining clouds will not develop between the ground
!     and the pressure at which the most entraining parcel initally
!     accelerates. BL turbulence or other sub-grid motions assumed
!     to maintain upward motion.
!
      if ( (.not. testlc) .and. (p .gt. pcsave) .and.    &
           (dwdp .gt. 0.) ) dwdp=0.
!
      west=wv+dwdp*dp
      if (west .lt. wdet) then
        drdp=0.
        return
      end if
      call lookup_es(test, es)
      rsc=epsilo*es/(p+(epsilo-1.)*es)
      htv=test*(1.+.609*rsc)
      pdap=fmp*RDGAS   *htv/((p+dp)*west)
      dadp=(pdap-(rda**2))/dp
      drdp=dadp/(2.*rda)

end subroutine simult




subroutine meens (cu,etfhr,rc,pr,phr,ps,pb, plzb, pt,gravit,cappa, &
                  epsilo, lat, t,q,apt,rlhr,emfhr,ca,tmel,xgcm,cuml, &
                  cmu,cmui,dmeml,emd,emdi,eme,emei,pmd,pztm,wmm, &
                  wmp,wtp, elt,contot, tmes,qmes,qtmes, &
                  umeml, debug_ijt, n2)


!------------------------------------------------------------------
logical, intent(in)          ::  debug_ijt        
integer, intent(in)          :: n2
real, dimension(:), intent(in) :: pr, phr, t, q, rlhr, emfhr
real, dimension(:), intent(out) :: cuml, cmu, dmeml, emd, eme, wmm, &
                                   wmp, umeml, tmes, qmes
real, dimension(:,:), intent(inout)  :: wtp,qtmes
real,         intent(in)  :: cu, rc, ps, pb, plzb, pt, gravit, cappa, &
                              lat,       apt, tmel
real, dimension(:,:), intent(in) :: etfhr
real, dimension(:,:), intent(in) :: xgcm

!! ADD THE FOLLOWING:
real, intent(in) :: epsilo
 


real, intent(inout) ::  ca
real, dimension(:), intent(inout) ::  elt
real,   intent(out) :: cmui, emei, contot, emdi, pmd, pztm
!------------------------------------------------------------------




!     Calculate mesoscale heat and moisture sources, using
!     variation on Leary and Houze (JAS, 1980).
!     Performs calculations after all sub-ensemble calculations complete.
!     For notation, see "Cu Closure A notes," 2/97
!
!     On Input:
!       cu      condenstation integral
!               sigma(i=1,N) (a(i,p_b)/a(1,p_b))*
!               (sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*c_u(i,j))
!               (mm/day)
!       etfhr   tracer flux convergence (kg/kg/s)-high resolution
!               sigma(i=1,N)(a(i,p_b)/a(1,p_b))*d(((r(i,j)**2)/(r(i,p_b        )**2))*
!               omega*'(i,j) * tracer*'(i,j))/dp
!       rc      precipitation integral
!               sigma(i=1,N) (a(i,p_b)/a(1,p_b))*
!               (sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*r_c(i,j))
!               (mm/day)
!       pr,phr  low-resolution pressure levels (Pa)
!       ps      surface pressure (Pa)
!       pb      cloud-base pressure (Pa)
!       plzb    level of zero buoyancy (Pa)
!       pt      cloud-top pressure (Pa)
!       gravit  gravity constant (m/s**2)
!       cappa   ratio of gas constant to specific heat for dry air
!       epsilo  ratio of gas constants, dry air to water vapor
!       rair    gas constant for dry air (J/kg/K)
!       cpair   specific heat for dry air (J/kg/K)
!       lat     latent heat for phase change (J/kg)
!       rh2o    gas constant for water vapor (J/kg/K)
!       t       low-resolution temperature (K)
!       q       low-resolution specific humidity (kg(H2O)/kg)
!       apt     mesoscale fraction, normalized by [a(1,p_b)*5]
!       rlhr    latent heat release (kg/kg/s)-high resolution
!               sigma(i=1,N)(a(i,p_b)/a(1,p_b))*((r(i,j)**2)/(r(i,p_b)**2))*
!               (sigma(k=1,4) gamma^*(k,i,j))
!       emfhr   moisture flux convergence (kg/kg/s)-high resolution
!               sigma(i=1,N)(a(i,p_b)/a(1,p_b))*d(((r(i,j)**2)/(r(i,p_b)**2))*
!               omega*'(i,j) * q*'(i,j))/dp
!       ca      condensed water X-fer from cells to anvil (mm/day)
!               weighted as cu,rc
!       tmel    melting temperature (K)
!       xgcm    low-resolution tracer mixing ratio (kg/kg)
!
!     On Output:
!       cuml    fractional mesoscale area, normalized by
!               a(1,p_b) at resolution of GCM
!       cmu     water mass condensed in mesoscale updraft
!               (g/kg/day) (normalized by a(1,p_b))
!       cmui    vertical integral of mesoscale-updraft deposition
!               (kg(H2O)/((m**2)*sec) 
!       dmeml   mass flux in mesoscale downdraft (kg/((m**2) s))
!               (normalized by a(1,p_b)) (index 1 at atmosphere bottom)
!               (resolution of GCM)
!       emd     water mass evaporated in mesoscale
!               downdraft (g/kg/day) (normalized by a(1,p_b))
!       emdi    vertical integral of mesoscale-downdraft sublimation
!               (mm/d)
!       eme     water mass evaporated from mesoscale
!               updraft (g/kg/day) (normalized by a(1,p_b))
!       emei    vertical integral of mesoscale-updraft sublimation
!               (kg(h2O)/((m**2)*sec)
!       pmd     pressure at top of mesoscale downdraft (Pa)
!       pztm    pressure at top of mesoscale updraft (Pa)
!       wmm     water vapor removal by condensation of
!               cell vapor source (g/kg/day) (normalized by a(1,p_b))
!       wmp     water vapor redistributed from cell vapor source
!               (g/kg/day) (normalized by a(1,p_b))
!       wtp     tracer redistributed by mesoscale processes
!               (kg/kg/s) (normalized by a(1,p_b))
!       elt     melting of ice in mesoscale updraft-
!               equivalent (g/kg/day)-which falls as meso sfc precip
!               (normalized by a(1,p_b))
!       contot  ratio of convective to total precipitation
!       tmes    temperature tendency due to mesoscale entropy-flux-
!               convergence (K/day) (normalized by a(1,p_b))
!       qmes    moisture tendency due to mesoscale moisture-flux
!               convergence (g/kg/day) (normalized by a(1,p_b))
!       qtmes   tracer tendency due to mesoscale tracer-flux
!               convergence (kg/kg/s) (normalized by a(1,p_b))
!       umeml   mass flux in mesoscale updraft (kg/((m**2) s))
!               (normalized by a(1,p_b)) (index 1 at atmosphere bottom)
!               (resolution of GCM)
!

!


      real, dimension (size(t(:))+1) :: emt, emq
      real, dimension (size(t(:))) :: owm, omv, emdx, tempq, tempqa, &
                                   tempt, otmv
      real, dimension (ncap) ::  wmhr, wphr, p, cumh, dmemh, wthrv
      real, dimension(size(t(:)),ncont) :: temptr,otm
      real, dimension(ncap,ncont) :: wthr
      real, dimension(ncont)  :: q1t

      logical :: thetl, lcons
      integer :: kcont
      integer :: jk, i, jj, jsave, jkm, jkp, jksave, j, k
      real :: dp, ptt, po, gnu, alrat, berat, etrat, alpha, beta, eta, &
              gnum, pzm, cmfhr, ome, ampt, pc1, pc2, omer, pctm, q1, &
              q4, es, qsat, q3, q5, anv, qref, qs, alp, pp, pm, pre, &
              qprip, qprim, eqfp, eqfm, hfmin, qmu, hflux, pfmin, omd, &
              c2, c3, c1, fjk, fjkm, qb, fjkb, qbm, qmd, qsmd, fjkmd, &
              qmmd, pi, qten, psa, owms, a, b, p1, emea, p3, pst, emda,&
              p2, tpri, sum, sbl, tme, wmpt, ta, wmmt, te, tep, sumq,  &
              tmu, targ, tprimd, tb, qsb, wa, wb, tmd, rin, tten, wmc, &
              wpc, rm, rmm, rma
      real :: wtpt, qtprip, qtprim, eqtfp, eqtfm






      if (debug_ijt) then
        write (unit_dc(n2), '(a)') 'in meens: entering meens'
      endif


!
!
!      define constants
!
      dp=-1000.
      ptt=pt+dp
      pztm=plzb
      tpri=1.
      po=100.e03
!
!     Restrict pztm to .ge. 10 kPa, cf Ackerman et al (JAS,1988)
!     (unless pt .le. 10kPa)
!     Stratospheric water vapor too high in AM2p9 with this pztm. 
!     use pztm >= plzb+dp
!
      if (pztm .lt. plzb) pztm=plzb
      if (ptt .lt. 10.e03  ) pztm=plzb+dp
      thetl=.false.
!
!      water-budget constants calculated from flux-
!      condensation parameterization
! 
       if (debug_ijt) then
         write (unit_dc(n2), '(a, 2e20.12)') 'in meens: rc,cu= ',rc,cu
         write (unit_dc(n2), '(a,  e20.12)') 'in meens: pztm = ',pztm 
       endif


      gnu=rc/cu
!
!     maintain Leary and Houze ratios of alpha, beta, and
!     eta to their sum
!
      sum=1.-gnu
      alrat=.25
      berat=.13
      etrat=.62
      alpha=alrat*sum
      beta=berat*sum
      eta=etrat*sum
!
!     use Leary and Houze ratios for gnu-sub-m
!
      gnum=.5

       if (debug_ijt) then
         write (unit_dc(n2), '(a, e20.12)') 'in meens: gnu= ',gnu
       endif

!
!     calculate mass of cu incorporated into mesoscale anvil,
!     integrated water evaporated in convective downdraft,
!     and integrated water evaporated from convective updraft
 
       if (debug_ijt) then
         write (unit_dc(n2), '(a, 3e20.12)')  &
                           'in meens: alpha,beta,eta= ',alpha,beta,eta
       endif


      if (ca .ge. 0.) ca=eta*cu
      if (ca .lt. 0.) ca=-ca

       if (debug_ijt) then
         write (unit_dc(n2), '(a, e20.12)') 'in meens: ca= ',ca
       endif

!
!     condensation of tower vapor source
!
      emt(nlev+1)=0.
      emq(nlev+1)=0.
      do 18 jk=1,nlev
         emt(jk)=0.
         emq(jk)=0.
         emd(jk)=0.
         eme(jk)=0.
         qmes(jk)=0.
         tmes(jk)=0.
         omv(jk)=0.
         tempqa(jk)=q(jk)
         wmp(jk)=0.
         wmm(jk)=0.
         cmu(jk)=0.
         tempq(jk)=q(jk)
         tempt(jk)=t(jk)
         cuml(jk)=0.
         dmeml(jk)=0.
         umeml(jk)=0.
        if (do_donner_tracer) then
         do kcont=1,ncont
          wtp(jk,kcont)=0.
          qtmes(jk,kcont)=0.
         temptr(jk,kcont)=xgcm(jk,kcont)
       end do
       endif
 18   continue
      do 12 i=1,ncap
      wmhr(i)=0.
      wphr(i)=0.
      cumh(i)=0.
      dmemh(i)=0.
      wthr(i,:)=0.
      p(i)=pb+(i-1)*dp
 12   continue
      pzm=0.
      do 11 i=1,ncap



      cmfhr=-rlhr(i)+emfhr(i)
      if (cmfhr .gt. 0.) then
         wmhr(i)=-cmfhr
         if (do_donner_tracer) then
         do kcont=1,ncont
           if (etfhr(i,kcont) .gt. 0.) then
             wthr(i,kcont)=-etfhr(i,kcont)
           end if
         end do
         endif
         if (pzm .eq. 0.) then
            pzm=p(i)
         end if
      end if


 11   continue
      do i=1,ncap

       if (debug_ijt) then
         write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                     'in meens: i,rlhr,emfhr= ',i,rlhr(i),emfhr(i)
       endif

       if (debug_ijt) then
         write (unit_dc(n2), '(a, i4, e20.12)')  &
                      'in meens: i,wmhr= ',i,wmhr(i)
       endif

      end do

      if (pzm .eq. 0.) pzm=pt
      if (pzm .le. pztm) pzm=pztm-dp
      sbl=0.
      lcons = .false.
      call verav(wmhr,pb,pt,lcons, thetl,owm,sbl,ps,pr,phr,cappa,  &
                 debug_ijt,n2)
      if (do_donner_tracer) then
      do kcont=1,ncont
         wthrv(:)=wthr(:,kcont)
         call verav(wthrv,pb,pt,lcons, thetl,otmv,sbl,ps,pr,phr,cappa, &
                    debug_ijt,n2) 
         otm(:,kcont)=otmv(:)
      end do
       endif

       if (debug_ijt) then
         write (unit_dc(n2), '(a, f19.10)') 'in meens: pzm= ',pzm
         do jk=1,nlev
!          print *, 'DONNER_DEEP/meens: jk,owm= ',jk,owm(jk)
           write (unit_dc(n2), '(a, i4, e20.12)')  &
                                  'in meens: jk,owm= ',jk,owm(jk)
         end do
       endif


      ome=-.463
      ampt=5.*apt
      tme=6.48e04
!
!     calculate redistribution of Cu H2O source
!
      do 13 jk=1,nlev
         if (owm(jk) .ge. 0.) go to 16
         pc1=phr(jk)
         pc2=phr(jk+1)
         omer=ome
         pctm=pc2+ome*tme
         if (pctm .le. pztm) omer=(pztm-pc2)/tme
         pctm=pc2+omer*tme
         q1=owm(jk)*(pc2-pc1)*tme/(pc1-pc2-omer*tme)
         q4=q1/2.
         if (do_donner_tracer) then
         do kcont=1,ncont
           q1t(kcont)=otm(jk,kcont)*(pc2-pc1)*tme/(pc1-pc2-  &
                      omer*tme)
         end do
         endif

         if (debug_ijt) then
           write (unit_dc(n2), '(a, 3e20.12)')  &
                         'in meens: pctm,pztm,q4= ',pctm,pztm,q4
         endif

         do 17 jj=jk,nlev



            if (phr(jj) .lt. pctm) go to 16
            tempq(jj)=tempq(jj)+(q1/ampt)
            tempqa(jj)=tempqa(jj)+(q4/ampt)
            wmpt=(q1/tme)
            if (do_donner_tracer) then
            do kcont=1,ncont
              temptr(jj,kcont)=temptr(jj,kcont)+(q1t(kcont)/  &
                               (2.* ampt))
              wtpt=q1t(kcont)/tme
              if (phr(jj+1) .le. pctm) wtpt=wtpt*  &
                                            (phr(jj)-pctm)/  &
                                             (phr(jj)-phr(jj+1))
              wtp(jj,kcont)=wtpt+wtp(jj,kcont)
            end do
            endif
            if (phr(jj+1) .le. pctm) wmpt=wmpt*(phr(jj)-pctm)/   &
                                          (phr(jj)-phr(jj+1))
            wmp(jj)=wmpt+wmp(jj)


 17    continue

       do jj=jk,nlev
         if (debug_ijt) then
           write (unit_dc(n2), '(a, i4, f19.10)') &
                            'in meens: jj,pr= ',jj,pr(jj)
         endif
         if (phr(jj) .lt. pctm) go to 216
         if (debug_ijt) then
           write (unit_dc(n2), '(a, i4, 3e20.12)')  &
                  'in meens: jj,q1,tempq,wmm= ',jj,q1,tempq(jj),wmm(jj)
           write (unit_dc(n2), '(a, e20.12)')  &
                   'in meens: wmp= ',wmp(jj)
           write (unit_dc(n2), '(a, i4, e20.12)')  &
                   'in meens: jj,tempqa= ',jj,tempqa(jj)
         endif

       end do

       if (debug_ijt) then
         write (unit_dc(n2), '(a, i4, 3e20.12)')  &
                  'in meens: jk,q1,tempq,wmm= ',jk,q1,tempq(jk),wmm(jk)
         write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                   'in meens: jk,wmp,owm= ',jk,wmp(jk),owm(jk)
       endif

 216   continue
 16   continue


      if (debug_ijt) then
        write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                        'in meens: jk,wmp,owm= ',jk,wmp(jk),owm(jk)
      endif


      wmp(jk)=wmp(jk)+owm(jk)
       if (do_donner_tracer) then
      wtp(jk,:)=wtp(jk,:)+otm(jk,:)
        endif
 13   continue
!
!     calculate condensed portion of redistributed H2O
!
      do 23 jk=1,nlev
        if ( (phr(jk+1) .le. pzm) .and. (phr(jk) .ge. pztm) ) then
          pc2=phr(jk+1)
          omer=ome
          pctm=pc2+ome*tme
          if (pctm .le. pztm) omer=(pztm-pc2)/tme
          pctm=pc2+omer*tme
          ta=t(jk)+tpri
          call lookup_es(ta, es)
          qsat=epsilo*es/(pr(jk)+(epsilo-1.)*es)
          q3=qsat-tempq(jk)
          q5=qsat-tempqa(jk)
          if (q3 .le. 0.) then
            wmmt=q3*ampt/tme
            if (phr(jk+1) .le. pctm) wmmt=wmmt*(phr(jk)-pctm)/  &
                                          (phr(jk)-phr(jk+1))
            wmm(jk)=wmmt
          end if
          if (q5 .le. 0.) tempqa(jk)=qsat
        end if
 23   continue
!
!     calculate meso up con portion due to meso lift
!
         anv=0.
         do 14 jk=1,nlev
            if (pr(jk) .gt. pzm) go to 19
            if (anv .eq. 0.) qref=tempqa(jk)
            anv=1.
            if (pr(jk) .lt. pztm) go to 21
            te=t(jk)+tpri
            call lookup_es(te, es)
            qs=epsilo*es/(pr(jk)+(epsilo-1.)*es)
            jsave=jk

            if (debug_ijt) then
              write (unit_dc(n2), '(a, 2e20.12)')  &
                          'in meens: qs,tempqa= ',qs,tempqa(jk)
            endif


            if (qref .ge. qs) go to 21
 19         continue
 14      continue
 21     continue
        alp=6.*ome/((pzm-pztm)**2)
        do 22 jk=1,nlev
          if (jk .eq. 1) then
            jkm=jk
          else
            jkm=jk-1
          end if
          if (jk .eq. nlev) then
            jkp=jk
          else
            jkp=jk+1
          end if
          if (pr(jk) .lt. pztm) go to 20
          if (pr(jk) .gt. pzm) go to 24
          pp=phr(jk+1)
          pm=phr(jk)
          if (phr(jk+1) .lt. pztm) pp=pztm
          if (phr(jk) .gt. pzm) pm=pzm
!
!        Calculate mesoscale entropy-flux convergence
!        Analytic integration used, possible only because
!        mesoscale temperature perturbation is not function of pressure
!        See "Vertical Velocity in Mesoscale Cloud" notes, 11/12/91
!
         tmes(jk)=(pzm+pztm)*(RDGAS   -CP_AIR)*(pp-pm)/CP_AIR
         tmes(jk)=tmes(jk)+((2.*CP_AIR-RDGAS   )*((pp**2)-(pm**2))/   &
                  (2.*CP_AIR))
         tmes(jk)=tmes(jk)-(RDGAS   *pztm*pzm/CP_AIR)*alog(pp/pm)
         tmes(jk)=tmes(jk)/(phr(jk+1)-phr(jk))
         tmes(jk)=tmes(jk)*ampt*tpri*alp
!
!         Calculate mesoscale vertical velocities
!
            omv(jk)=(pzm+pztm)*((pp**2)-(pm**2))/2.
            omv(jk)=omv(jk)-(((pp**3)-(pm**3))/3.)
            omv(jk)=omv(jk)-pztm*pzm*(pp-pm)
            omv(jk)=omv(jk)/(phr(jk+1)-phr(jk))
            omv(jk)=omv(jk)*alp
            if (jk .lt. jsave) go to 24
            pre=pr(jk)
            te=t(jk)+tpri
            call lookup_es(te, es)
            tempqa(jk)=epsilo*es/(pr(jk)+(epsilo-1.)*es)
            if (qref .ge. tempqa(jk)) then
              tep=t(jkp)+tpri
              if (pr(jkp) .le. pztm) then
                cmu(jk)=-omv(jk)*(tempqa(jk)-tempqa(jkm))/(pr(jk)   &
                         -pr(jkm))
              else if (jk .eq. jsave) then
                call lookup_es(tep, es)
                tempqa(jkp)=epsilo*es/(pr(jkp)+(epsilo-1.)*   &
                            es)
                cmu(jk)=-omv(jk)*(tempqa(jkp)-tempqa(jk))/(pr(jkp)   &
                        -pr(jk))
                qref=tempqa(jkp)
              else
                call lookup_es(tep, es)
                tempqa(jkp)=epsilo*es/(pr(jkp)+(epsilo-1.)*es)
                cmu(jk)=-omv(jk)*(tempqa(jkp)-tempqa(jkm))/(pr(jkp)   &
                       -pr(jkm))
                qref=tempqa(jkp)
              end if
              if (cmu(jk) .lt. 0.) cmu(jk)=0.
            else
              cmu(jk)=0.
            end if
            cmu(jk)=cmu(jk)*ampt*8.64e07

            if (debug_ijt) then
              write (unit_dc(n2), '(a, i4, f20.14, e20.12)') &
                        'in meens: jk,t,omv= ',jk,t(jk),omv(jk)
            endif


 24       continue
 22     continue
 20    continue
!
!        Calculate mesoscale moisture-flux convergence
!
        sumq=0.
        do 25 jk=1,nlev 
          if (jk .eq. 1) then
            jkm=jk
          else
            jkm=jk-1
          end if
          if (jk .eq. nlev) then
            jkp=jk
          else
            jkp=jk+1
          end if
          if (phr(jk) .lt. pztm) go to 26
          if (phr(jk+1) .gt. pzm) go to 27
          qprip=(tempqa(jkp)+tempqa(jk)-q(jkp)-q(jk))/2.
          qprim=(tempqa(jk)+tempqa(jkm)-q(jk)-q(jkm))/2. 
          if (do_donner_tracer) then
          do kcont=1,ncont
            qtprip=(temptr(jkp,kcont)+temptr(jk,kcont)  &
                    -xgcm(jkp,kcont)-xgcm(jk,kcont))/2.
            qtprim=(temptr(jk,kcont)+temptr(jkm,kcont) &
                  -xgcm(jk,kcont)-xgcm(jkm,kcont))/2.
            eqtfp=ampt*qtprip*alp*(phr(jk+1)-pztm)*  &
                  (pzm-phr(jk+1))
            eqtfm=ampt*qtprim*alp*(phr(jk)-pztm)*  &
                  (pzm-phr(jk))
            if ((phr(jk) .le. pzm) .and. (phr(jk+1) .ge. pztm)) then
              qtmes(jk,kcont)=(eqtfm-eqtfp)/(phr(jk+1)-phr(jk))
            end if
            if ((pzm .le. phr(jk)) .and. (pzm .ge. phr(jk+1))) then
              qtmes(jk,kcont)=eqtfp/(phr(jk)-phr(jk+1))
            end if
            if ((pztm .ge. phr(jk+1)) .and. (pztm .le. phr(jk))) then
              qtmes(jk,kcont)=eqtfm/(phr(jk+1)-phr(jk))
              if ((pzm .le. phr(jk)) .and. (pzm .ge. phr(jk+1))) then
                qtmes(jk,kcont)=0.
              end if
            end if ! ((pztm .ge. phr(jk+1)) .and. (pztm .le. phr(jk)))
          end do
          endif
          eqfp=ampt*qprip*alp*(phr(jk+1)-pztm)*(pzm-phr(jk+1))
          eqfm=ampt*qprim*alp*(phr(jk)-pztm)*(pzm-phr(jk))
          if ((phr(jk) .le. pzm) .and. (phr(jk+1) .ge. pztm)) then
            qmes(jk)=(eqfm-eqfp)/(phr(jk+1)-phr(jk))
          end if
          if ((pzm .le. phr(jk)) .and. (pzm .ge. phr(jk+1))) then
            qmes(jk)=eqfp/(phr(jk)-phr(jk+1))
          end if
          if ((pztm .ge. phr(jk+1)) .and. (pztm .le. phr(jk))) then
            qmes(jk)=eqfm/(phr(jk+1)-phr(jk))
            if ((pzm .le. phr(jk)) .and. (pzm .ge. phr(jk+1))) then
              qmes(jk)=0.
            end if
          end if ! ((pztm .ge. phr(jk+1)) .and. (pztm .le. phr(jk)))
 27     continue
 25   continue
 26   continue
!
!     Calculate eddy flux of moist static energy in mesoscale
!     updraft and identify its minimum.
!
      hfmin=0.
      do 41 jk=1,nlev
         if (pr(jk) .lt. pztm) go to 42
         if (pr(jk) .gt. pzm) go to 43
         tmu=t(jk)+tpri
         call lookup_es(tmu, es)
         qmu=epsilo*es/(pr(jk)+(epsilo-1.)*es)
         hflux=omv(jk)*(CP_AIR*tpri+lat*(qmu-q(jk)))
         if (hflux .lt. hfmin) then
            hfmin=hflux
            pfmin=pr(jk)
         end if
 43      continue
 41   continue
 42   continue

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')  &
                        'in meens: hfmin,pfmin= ',hfmin,pfmin
      endif


!
!     Calculate mesoscale downdraft speed at pmd. Follow Leary
!     and Houze (1980,JAS) and set magnitude to half that in
!     mesoscale updraft; vertical pressure velocity constant with ht. 
!
      omd=-alp*((pzm-pztm)**2)/8.
      omd=omd/2.
      pmd=pzm+20.e03
      if (pmd .gt. ps) pmd=ps
      if (pmd .gt. ps) pmd=ps
!
!     Calculate temperature and specific humidity in mesoscale
!     downdraft. 
!
      do 44 jk=1,nlev
         if (pr(jk) .lt. pmd) go to 45
         if (pr(jk) .gt. pb) go to 46
!
!     Calculate c2, the relative humidity in the mesoscale downdraft,
!     after Table 3 of Leary and Houze (1980, JAS).
!
         c2=1.-(.3*(pr(jk)-pmd)/(pb-pmd))
!
!     Calculate c3, the factor which yields the eddy flux of moist
!     static energy when multiplied by the minimum of moist static
!     energy in the mesoscale updraft. Multiply by 1.3 to take account
!     of convective downdrafts. See Fig. 7 of Leary and Houze
!     (1980,JAS).
!
         c3=(pr(jk)-pmd)/(pb-pmd)
         c3=1.3*c3
!
!     See "Moist Static Energy A, 1/26/91" notes.
!
         targ=t(jk)
         call lookup_es(targ, es)
         qs=epsilo*es/(pr(jk)+(epsilo-1.)*es)
         c1=epsilo*lat*es/(pr(jk)*rh2o*(t(jk)**2))
         tprimd=c3*hfmin/omd
         tprimd=tprimd-lat*(c2*qs-q(jk))
         tprimd=tprimd/(CP_AIR+lat*c1*c2)
         tempt(jk)=t(jk)+tprimd
         targ=tempt(jk)
         call lookup_es(targ, es)
         tempqa(jk)=c2*es*epsilo/(pr(jk)+(epsilo-1.)*es)

         if (debug_ijt) then
           write (unit_dc(n2), '(a, 4e20.12)')  &
                     'in meens: tprimd,tempqa,q,qs= ',tprimd,   &
                     tempqa(jk),q(jk),qs
!          print *, 'DONNER_DEEP/meens: pr,rh,factr= ',pr(jk),c2,c3
           write (unit_dc(n2), '(a, f19.10, 2e20.12)')  &
                     'in meens: pr,rh,factr= ',pr(jk),c2,c3
         endif


 46    continue
 44   continue
 45   continue
!
!     Calculate eddy fluxes of potential temperature and specific
!     humidity in mesoscale downdraft.
!
      do 28 jk=2,nlev-1
         if (phr(jk) .lt. pmd) go to 47
         if (phr(jk) .gt. pb) go to 29
         if ((pr(jk-1) .le. pb) .and. (pr(jk) .ge. pmd)) then
            fjk=ampt*omd*((po/pr(jk))**(RDGAS   /CP_AIR))  &
            *(tempt(jk)-t(jk))    
            fjkm=ampt*omd*((po/pr(jk-1))**(RDGAS   /CP_AIR))*   &
            (tempt(jk-1)-t(jk-1))
            emt(jk)=(fjk+fjkm)/2.
            fjk=ampt*omd*(tempqa(jk)-q(jk))
            fjkm=ampt*omd*(tempqa(jk-1)-q(jk-1))
            emq(jk)=(fjk+fjkm)/2.
         end if
         if (pr(jk-1) .ge. pb) then
            fjk=ampt*omd*((po/pr(jk))**(RDGAS   /CP_AIR))*   &
            (tempt(jk)-t(jk))
            call polat(q,pr,qb,pb, debug_ijt, n2)
            call polat(t,pr,tb,pb, debug_ijt, n2)
            call lookup_es(tb, es)
            qsb=epsilo*es/(pb+(epsilo-1.)*es)
            tprimd=hfmin/omd
            tprimd=tprimd-lat*(.7*qsb-qb)
            c1=epsilo*lat*es/(pb*rh2o*(tb**2))
            tprimd=tprimd/(CP_AIR+.7*lat*c1)
            fjkb=ampt*omd*((po/pb)**(RDGAS   /CP_AIR))*tprimd
            wa=(phr(jk)-pr(jk))/(pb-pr(jk))
            wb=(pb-phr(jk))/(pb-pr(jk))
            emt(jk)=wa*fjkb+wb*fjk
            fjk=ampt*omd*(tempqa(jk)-q(jk))
            targ=tb+tprimd
            call lookup_es(targ, es)
            qbm=.7*epsilo*es/(pb+(epsilo-1.)*es)
            fjkb=ampt*omd*(qbm-qb)
            emq(jk)=wa*fjkb+wb*fjk
         end if
         if (pr(jk) .le. pmd) then
            fjkm=ampt*omd*((po/pr(jk-1))**(RDGAS   /CP_AIR))*   &
            (tempt(jk-1)-t(jk-1))
            call polat(q,pr,qmd,pmd, debug_ijt,n2)
            call polat(t,pr,tmd,pmd, debug_ijt, n2)
            call lookup_es(tmd, es)
            qsmd=epsilo*es/(pmd+(epsilo-1.)*es)
            c1=epsilo*lat*es/(pmd*rh2o*(tmd**2))
            tprimd=-lat*(qsmd-qmd)/(CP_AIR+lat*c1)
            fjkmd=ampt*omd*((po/pmd)**(RDGAS   /CP_AIR))*tprimd
            wa=(pr(jk-1)-phr(jk))/(pr(jk-1)-pmd)
            wb=(phr(jk)-pmd)/(pr(jk-1)-pmd)
            emt(jk)=fjkmd*wa+fjkm*wb
            targ=tmd+tprimd
            call lookup_es(targ, es)
            qmmd=epsilo*es/(pmd+(epsilo-1.)*es)
            fjkm=ampt*omd*(tempqa(jk-1)-q(jk-1))
            fjkmd=ampt*omd*(qmmd-qmd)
            emq(jk)=fjkmd*wa+fjkm*wb
          end if

          if (debug_ijt) then
            write (unit_dc(n2), '(a, i4, 3e20.12)')  &
                        'in meens: jk,phr,emt,emq= ',jk,phr(jk),   &
                         emt(jk),emq(jk)
          endif


          emt(jk)=((po/pr(jk))**(RDGAS   /CP_AIR))*emt(jk)
 29       continue
 28   continue
 47   continue
!
!     Calculate temperature and specific humidity tendencies due
!     to eddy-flux convergences in mesoscale downdraft.
!
       rin=0.
       do 31 jk=nlev,1, -1
         if ((phr(jk+1) .le. pzm) .and. (phr(jk) .ge. pzm))  &
         jksave=jk+1
         pi=(po/pr(jk))**(RDGAS   /CP_AIR)
         if ((emt(jk+1) .ne. 0.) .and. (emt(jk) .eq. 0.)   &
              .and. (rin .eq. 0.)) then
            tten=-emt(jk+1)/(phr(jk+1)-ps)
            qten=-emq(jk+1)/(phr(jk+1)-ps)
            rin=1.
         end if
         if (rin .eq. 1.) then
            tmes(jk)=tmes(jk)+(tten/pi)
            qmes(jk)=qmes(jk)+qten
         end if
         if ((rin .eq. 0.) .and. (emt(jk+1) .ne. 0.) .and.   &
              (emt(jk) .ne. 0.)) then
            tten=(emt(jk+1)-emt(jk))/(phr(jk+1)-phr(jk))
            tten=-tten/pi
            qten=(emq(jk+1)-emq(jk))/(phr(jk+1)-phr(jk))
            qten=-qten
            tmes(jk)=tmes(jk)+tten
            qmes(jk)=qmes(jk)+qten
         end if

         if (debug_ijt) then
           write (unit_dc(n2), '(a, i4, f19.10, f20.14, e20.12)')   &
                    'in meens: jk,pr,tmes,qmes= ',jk,pr(jk),  &
                     tmes(jk),qmes(jk)
         endif


 31   continue
!
!     Apply flux at top of mesoscale downdraft to level between
!     pzm and pmd, as flux at bottom of mesoscale downdraft is
!     applied to PBL.
!
      psa=0.
      do 32 jk=1,nlev
         if ((emt(jk) .ne. 0.) .and. (emt(jk+1) .eq. 0.)) then
            tten=emt(jk)/(phr(jksave)-phr(jk))
            qten=emq(jk)/(phr(jksave)-phr(jk))
            psa=phr(jk)
         end if
 32      continue

         if (debug_ijt) then
           write (unit_dc(n2), '(a, 2f19.10)')  &
                                 'in meens: pmd,pb= ',pmd,pb
         endif

         do jk=1,nlev
           if ((pr(jk) .le. psa) .and. (pr(jk) .ge. phr(jksave))) then

             if (debug_ijt) then
               write (unit_dc(n2), '(a, 3e20.12)')  &
                   'in meens: po,psa,phr(jksave)= ',po,psa,phr(jksave)
!              print *, 'DONNER_DEEP/meens: jk,qmes,qten= ',jk, &
!                        qmes(jk),qten
               write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                         'in meens: jk,qmes,qten= ',jk,qmes(jk),qten
!              print *, 'DONNER_DEEP/meens: jk,tmes,tten= ',jk, &
!                         tmes(jk),tten
               write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                           'in meens: jk,tmes,tten= ',jk,tmes(jk),tten
             endif

             qmes(jk)=qmes(jk)+qten
             pi=(po/pr(jk))**(RDGAS   /CP_AIR)
             tmes(jk)=tmes(jk)+(tten/pi)
           end if
         end do
      cmui=0.
      wmc=0.
      wpc=0.
      owms=0.
      do 15 jk=1,nlev
        tmes(jk)=tmes(jk)*86400.
        qmes(jk)=qmes(jk)*8.64e07
        wmm(jk)=wmm(jk)*8.64e07
        owm(jk)=owm(jk)*8.64e07
        wmp(jk)=wmp(jk)*8.64e07
        wmc=wmc+wmm(jk)*(phr(jk)-phr(jk+1))
        owms=owms+owm(jk)*(phr(jk)-phr(jk+1))
        wpc=wpc+wmp(jk)*(phr(jk)-phr(jk+1))
        cmui=cmui+cmu(jk)*(phr(jk)-phr(jk+1))
 15   continue
      wmc=wmc/(gravit*1000.)
      wpc=wpc/(gravit*1000.)
      owms=owms/(gravit*1000.)
      cmui=cmui/(gravit*1000.)

       if (debug_ijt) then
         write (unit_dc(n2), '(a, e20.12,a,a,e20.12,a)')  &
                    'in meens: wmc=',wmc,' mm/day',' wpc=',wpc, 'mm/day'
         write (unit_dc(n2), '(a, e20.12, a, a, e20.12, a)')  &
                      'in meens: owms= ',owms,' mm/day',' cmui= ',   &
                        cmui,'mm/day'
       endif


! 
!
!     calculate mesoscale precipitation
!
      cmui=cmui-wmc
      rm=gnum*(cmui+ca)
      contot=rc/(rm+rc)

       if (debug_ijt) then
         write (unit_dc(n2), '(a, 2e20.12)')  &
                                 'in meens: cmui,ca=', cmui,ca
!        print *, 'DONNER_DEEP/meens: rm= ',rm,'rc= ',rc
         write (unit_dc(n2), '(a, e20.12, a, e20.12)')  &
                                    'in meens: rm= ',rm,'rc= ',rc
       endif
!
!      calculate integrated water evaporated in mesoscale  
!      downdraft using Leary and Houze coefficients
!
      a=.4
      emdi=a*(cmui+ca)
      b=.1
      emei=b*(cmui+ca)
!
!     calculate vertical tendency distributions (g/kg/day)
!
!      evaporation from mesoscale updraft
!
      p1=pzm
      emea=emei*gravit*1000./(p1-pztm)

       if (debug_ijt) then
         write (unit_dc(n2), '(a, 2e20.12)')  &
                          'in meens: emea,emei= ',emea,emei
       endif


      if (p1 > pztm) then
        call ver(emea,p1,pztm,pr,phr,eme, debug_ijt, n2)
      endif
!
!
      p3=pzm
      pst=pzm+30.e03
      pst=ps
      p1=pst
      emda=emdi*gravit*1000./(p1-p3)
!
!      evaporation in mesoscale downdraft
!
       do 35 jk=1,nlev
          if (phr(jk) .lt. p3) go to 36
          if (phr(jk+1) .gt. pst) go to 37
          pm=phr(jk)
          pp=phr(jk+1)
          if ((phr(jk) .ge. pst) .and. (phr(jk+1) .le.   &
          pst)) pm=pst
          if ((phr(jk) .ge. p3) .and. (phr(jk+1) .le. p3))  &
          pp=p3
          emd(jk)=emda*(pm-pp)*(pm+pp-2.*pst)

          if (debug_ijt) then
            write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                               'in meens: jk,emda,emd= ',jk,emda,emd(jk)
          endif

          emd(jk)=emd(jk)/((phr(jk)-phr(jk+1))*(p3-pst))
 37   continue
 35   continue
 36   continue

         if (debug_ijt) then
           write (unit_dc(n2), '(a, e20.12, f19.10)')  &
                                      'in meens: emdi,ps= ',emdi,ps
         endif


!
!      evaporation in mesoscale downdraft
!

       if (debug_ijt) then
         write (unit_dc(n2), '(a, e20.12, f19.10)')  &
                                 'in meens: emdi,ps= ',emdi,ps
       endif



       if (debug_ijt) then
         write (unit_dc(n2), '(a, e20.12)') 'in meens: pzm= ',pzm
       endif



!      if (debug_ijt) then
!        print *, 'DONNER_DEEP/meens: factsum= ',factsum
!      endif


!
!     equivalent for melting of mesoscale sfc precip
!     generated by deposition in mesoscale updraft
!
      p2=-10.
      do 51 j=1,nlev-1
         if (phr(j+1) < pztm ) go to 52
         if ((t(j) .ge. tmel) .and. (t(j+1) .le. tmel))  then
            p2=phr(j+1)
            go to 52
         end if
 51   continue
 52   continue
      if (p2 .ne. -10.) then
        rm=(rc/contot)-rc
      else
        rm=0.
      end if
      rmm=rm
      rma=rmm*gravit*1000./(pb-p2)
      if (pb > p2) then
        call ver(rma,pb,p2,pr,phr,elt, debug_ijt,n2)
      endif
!
!     Convert cmui and emei from mm/day to kg(H2O)/((m**2) sec)
      cmui=cmui/86400.
      emei=emei/86400.
      do j=1,ncap
        if ( (p(j) .le. pzm) .and. (p(j) .ge. pztm)) cumh(j)=ampt
        if ( (p(j) .le. pb) .and. (p(j) .ge. pmd) ) dmemh(j)=    &
                                                    -omd*ampt/gravit
      end do
      lcons = .false.
      call verav(cumh,pb,pztm,lcons,thetl,cuml,sbl,ps,pr,phr,cappa,  &
                 debug_ijt,n2)
      call verav(dmemh,pb,pmd,lcons,thetl,dmeml,sbl,ps,pr,phr,cappa, &
                 debug_ijt,n2)
      do j=1,nlev
        umeml(j)=-omv(j)*ampt/gravit
      end do

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 4e20.12)') &
                          'in meens: pzm,pztm,pb,p2= ',pzm,pztm,pb,p2
        do k=1,nlev-kstart_diag+1
          if (omv(k) /= 0.00) then
!            print *, 'DONNER_DEEP/meens: j,omv= ',k,omv(k)
             write (unit_dc(n2), '(a, i4, e20.12)')  &
                            'in meens: j,omv= ',k,omv(k)
          endif
          if (elt(k) /= 0.0) then 
!           print *, 'DONNER_DEEP/meens: j,elt= ',k,elt(k)
            write (unit_dc(n2), '(a, i4, e20.12)')  &
                                   'in meens: j,elt= ',k,elt(k)
          endif
        end do
      endif



end subroutine meens





subroutine mesub (cu,rc,dint,plzb,pr,phr,ps,pb,pt,gravit,cappa, t,ca, &
                  tmel,ecd,ece,fre,elt, debug_ijt, n2)

!--------------------------------------------------------------------
integer, intent(in) :: n2
real,               intent(inout) :: dint
real,               intent(in) :: cu, rc,plzb,  ps, pb, pt, gravit, &
                                  cappa,                   tmel
real, dimension(:), intent(in) :: pr, phr, t
logical,            intent(in) :: debug_ijt
real, dimension(:), intent(inout) :: ecd, ece, fre, elt
real,               intent(inout) :: ca
!---------------------------------------------------------------------

!
!     Calculate mesoscale heat and moisture sources, using
!     variation on Leary and Houze (JAS, 1980).
!     Performs calculations related to individual subensembles only.
!     For notation, see "Cu Closure A notes," 2/97
!
!     On Input:
!
!       cu      condenstation integral
!               sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*c_u(i,j)
!               for ensemble i
!               (mm/day)
!       rc      precipitation integral
!               sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*r_c(i,j)
!               for ensemble i
!               (mm/day)
!       pr,phr  low-resolution pressure levels (Pa)
!       ps      surface pressure (Pa)
!       pb      cloud-base pressure (Pa)
!       pt      cloud-top pressure (Pa)
!       gravit  gravity constant (m/s**2)
!       cappa   ratio of gas constant to specific heat for dry air
!       epsilo  ratio of gas constants, dry air to water vapor
!       rair    gas constant for dry air (J/kg/K)
!       t       low-resolution temperature (K)
!       tmel    melting temperature (K)
!       dint    water mass frozen in convective updraft
!               plus ice deposited convective updraft
!               (kg(H2O)/((m**2) sec)
!               weighted as cu,rc
!
!     On Output:
!
!       ca      condensed water X-fer from cells to anvil (mm/day)
!               weighted as cu,rc
!       ecd     water mass evaporated in convective
!               downdraft (g/kg/day) (vertical integral ecdi
!               weighted as rc,cu)
!       ece     water mass evaporated from convective updraft
!               (g/kg/day) (vertical integral ecei weighted as cu,rc)
!       fre     equivalent for freezing in mesoscale updraft
!               (g/kg/day) (vertical integral caa weighted as cu,rc)
!       elt     equivalent for melting in mesoscale downdraft
!               (g/kg/day) (vertical integral elt weighted as cu,rc)
!

      integer :: i, nlevm, jk
      real  :: dp, ptt, pztm, pzm, gnu, sum, alrat, berat, etrat, &
               alpha, beta, eta, gnum, ecei, ecdi, caa, pz0, p1, ecda, &
               p2, elta, ecea

!
!      define constants
!
      dp=-1000.
      ptt=pt+dp
      pztm=ptt-30.e03
!     Restrict pztm to .ge. 10 kPa, cf Ackerman et al (JAS,1988)
!     (unless pt .le. 10 kPa)
!     Stratospheric water vapor too high in AM2p9 with this pztm. 
!     use pztm >= plzb+dp
!
!     if (pztm .lt. 10.e03) pztm=10.e03
      if (pztm .lt. plzb) pztm=plzb
!     if (ptt .lt. 10.e03) pztm=ptt+dp
      if (ptt .lt. plzb) pztm=plzb+dp
      pzm=ptt
      if (pzm .le. pztm) pzm=pztm-dp
!
!      water-budget constants calculated from flux-
!      condensation parameterization
! 

       if (debug_ijt) then
         write (unit_dc(n2), '(a, 2e20.12)') 'in mesub: rc,cu= ',rc,cu
         write (unit_dc(n2), '(a,  e20.12)') 'in mesub: plzb = ',plzb
       endif

      gnu=rc/cu
!
!     maintain Leary and Houze ratios of alpha, beta, and
!     eta to their sum
!
      sum=1.-gnu
      alrat=.25
      berat=.13
      etrat=.62
      alpha=alrat*sum
      beta=berat*sum
      eta=etrat*sum
!
!     use Leary and Houze ratios for gnu-sub-m
!
      gnum=.5

       if (debug_ijt) then
         write (unit_dc(n2), '(a, e20.12)')  'in mesub: gnu= ',gnu
       endif


!
!     calculate mass of cu incorporated into mesoscale anvil,
!     integrated water evaporated in convective downdraft,
!     and integrated water evaporated from convective updraft
!


      ca=eta*cu


       if (debug_ijt) then
         write (unit_dc(n2), '(a, 3e20.12)')  &
                          'in mesub: alpha,beta,eta= ',alpha,beta,eta
         write (unit_dc(n2), '(a, e20.12)') 'in mesub: ca= ',ca
       endif

      ecei=beta*cu
      ecdi=alpha*cu
!
!
!     calculate vertical tendency distributions (g/kg/day)
!
!      calculate equivalent for freezing in mesoscale updraft
!
      dint=dint*8.64e07*gravit/(pzm-pztm)
      if (dint .ne. 0.)    &
       caa=(ca+ecei)*gravit*1000./(pzm-pztm)
      if (dint .eq. 0.)   &
       caa=ca*gravit*1000./(pzm-pztm)

       if (debug_ijt) then
         write (unit_dc(n2), '(a, 2e20.12)')  &
                                 'in mesub: caa,dint=',caa,dint
       endif


      caa=caa-dint

      if (debug_ijt) then
!       if (caa .gt. 0.) then
        write (unit_dc(n2), '(a, 3e20.12)')  &
                            'in mesub: caa,pzm,pztm= ',caa,pzm,pztm
        if (caa .gt. 0.) then
          do i=1,nlev-kstart_diag+1
            if (fre(i) /= 0.0) then
              write (unit_dc(n2), '(a, i4, e20.12)')  &
                                     'in mesub: i,fre= ',i,fre(i)
            endif
          end do
        endif
      endif

      if ( pzm > pztm) then
        if (caa .gt. 0.)  &
           call ver(caa,pzm,pztm,pr,phr,fre, debug_ijt, n2)
!       if (caa .gt. 0.)  call ver(caa,pzm,pztm,   phr,fre, debug_ijt)
      endif
!
!     evaporation in convective downdraft
!
      pz0=ptt
      p1=ps
      ecda=ecdi*gravit*1000./(p1-pz0)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 3e20.12)')  &
                         'in mesub: ecda,p1,pz0= ',ecda,p1,pz0
        do i=1,nlev-kstart_diag+1
          if (ecd(i) /= 0.0) then
            write (unit_dc(n2), '(a, i4, e20.12)')  &
                                    'in mesub: i,ecd= ',i,ecd(i)
          endif
        end do
      endif


      if (p1 > pz0)then
        call ver(ecda,p1,pz0,pr,phr,ecd, debug_ijt,n2)
      endif
!
!     calculate melting due to excess freezing, over
!     ecei and ca
!
      nlevm=nlev-1
      p2=0.
      do 51 jk=1,nlev
         if (phr(jk+1) < pztm) exit
         if ((t(jk) .ge. tmel) .and. (t(jk+1) .le. tmel)) then
           p2=phr(jk+1)
           go to 52
         end if
 51   continue
 52   continue
      elta=0.
      if (caa .le. 0.) then 
         caa=-caa*(pzm-pztm)/(pb-p2)
         elta=caa
      end if
      if (p2 .eq. 0.) elta=0.

      if (debug_ijt) then
        write (unit_dc(n2), '(a, e20.12, 2f19.10)') &
                           'in mesub: elta,pb,p2= ',elta,pb,p2
        do i=1,nlev-kstart_diag+1
          if (elt(i) /= 0.0) then
            write (unit_dc(n2), '(a, i4, e20.12)')  &
                                    'in mesub: i,elt= ',i,elt(i)
          endif
        end do
      endif


      if (pb > p2) then
        call ver(elta,pb,p2,pr,phr,elt, debug_ijt, n2)
!       call ver(elta,pb,p2,   phr,elt, debug_ijt)
      endif
!
!     evaporation from convective updraft
!
      p1=pt+5.0e03
      p2=ptt
      ecea=ecei*gravit*1000./(p1-p2)

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')  &
                         'in mesub: ecea,ecei= ',ecea,ecei
        write (unit_dc(n2), '(a, e20.12, 2f19.10)')  &
                         'in mesub: ecea,p1,p2= ',ecea,p1,p2
        do i=1,nlev-kstart_diag+1
          if (ece(i) /= 0.0) then
            write (unit_dc(n2), '(a, i4, e20.12)')  &
                                     'in mesub: i,ece= ',i,ece(i)
          endif
        end do
      endif


      if (p1 > p2) then
        call ver(ecea,p1,p2,pr,phr,ece, debug_ijt, n2)
!       call ver(ecea,p1,p2,   phr,ece, debug_ijt)
      endif


end subroutine mesub



!####################################################################

subroutine micro(tc1,tc2,p1,p2,te1,te2,qe1,qe2,w1,w2, rr,rmu,  &
                 qrw,qcw,qlw,dcw1, dqrw3, debug_ijt, n2)

 !--------------------------------------------------------------------
integer, intent(in) :: n2
real,      intent(inout) :: qcw, qrw, rmu
real,        intent(in) :: tc1, tc2, p1, p2, te1, te2, qe1, qe2, &
                             w1, w2, rr
real,        intent(inout) :: qlw, dcw1, dqrw3
logical,      intent(in) :: debug_ijt
!--------------------------------------------------------------------




!
!     Kessler microphysics
!
!     On Input:
!        tc1    cloud temperature (K) at pressure p1 (Pa)
!        tc2    cloud temperature (K) at pressure p2 (Pa)
!        te1    environmental temperature (K) at pressure p1 (Pa)
!        te2    environmental temperature (K) at pressure p2 (Pa)
!        qe1    environmental mixing ratio (kg(H2O)/kg) at pressure p1
!        qe2    environmental mixing ratio (kg(H2O)/kg) at pressure p2
!        w1     cloud vertical velocity (m/s) at pressure p1
!        w2     cloud vertical velocity (m/s) at pressure p2
!        rr     cloud radius (m)
!        qrw    rain water (g(H2O)/m**3)
!        qcw    cloud water (g(H2O)/m**3)
!        rmu    entrainment coefficient (/m)
!
!     On Output:
!        qlw    total liquid water (kg(H2O)/kg)
!        dcw1   condensation increment (g(H2O)/m**3)
!        dqrw3  precipitation increment (g(H2O)/m**3)
!
      real :: ep, g, rd, alp, dp, epm, hr, es1, es2, rs1, rs2, tcb, &
              qcb, d1, d2, dz, qeb, rmusa, cond, w, dcw2, dqcw3, ent, &
              qcwa, red, qrwa, tv, teb, tve
!     CONSTANTS
!

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')   &
                                 'in micro: qrw,qcw= ',qrw,qcw
        write (unit_dc(n2), '(a, e20.12)')  'in micro: rr= ',rr
      endif

      ep=.622
      g=9.81
      rd=287.05
      alp=.5
      dp=-1000.
!
!     epm is defined in "Generalized mu 10/1/89"
!
      epm=0.
      hr=rr
!
!     call establ(es1,tc1)
!     call establ(es2,tc2)
!     call escomp(tc1, es1)
!     call escomp(tc2, es2)
      call lookup_es(tc1, es1)
      call lookup_es(tc2, es2)
      rs1=ep*es1/(p1-es1)
      rs2=ep*es2/(p2-es2)
      tcb=(tc1+tc2)/2.

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12, 2f19.10)') &
                         'in micro: rs1,rs2,p1,p2= ',rs1,rs2,p1,p2
        write (unit_dc(n2), '(a, 2e20.12, 2f20.14)')  &
                         'in micro: es1,es2,tc1,tc2= ',es1,es2,tc1,tc2
        write (unit_dc(n2), '(a, 2f20.14)')  &
                          'in micro: te1,te2= ',te1,te2
      endif


      qcb=(rs1+rs2)/2.
      d1=rd*(1.+alp)*tc1*te1/(g*p1*(alp*te1+tc1))
      d2=rd*(1.+alp)*tc2*te2/(g*p2*(alp*te2+tc2))
      dz=-(d1+d2)*dp/2
!
!     calculate condensation
!
      qeb=(qe1+qe2)/2.

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)')  &
                   'in micro: qeb,rmu= ',qeb,rmu
      endif


      rmusa=rmu
!     rmu=0.
      cond=rs2*exp(rmu*dz)-rs1+(1.-exp(rmu*dz))*qeb
      cond=-cond/(1.-epm+epm*exp(rmu*dz))
      if (cond .lt. 0.) cond=0.

      if (debug_ijt) then
        write (unit_dc(n2), '(a, e20.12)') 'in micro: cond= ',cond
      endif


      dcw1 =(cond)*(p1+p2)*500./(tcb*rd*(1.+.61*qcb))
      w=(w1+w2)/2.
      IF (QCW .GE. .5) DCW2=DZ*(QCW-.5)/W
      IF (QCW .GE. .5) DCW2=DCW2*1.E-03
      IF (QCW .LT. .5) DCW2=0.
      DQCW3=0.
      IF ((QCW .EQ. 0.) .OR. (QRW .EQ. 0.)) GO TO 7
      DQCW3=5.26E-03*QCW*(QRW**.875)*DZ/W
!
!     calculate effect of entrainment on cloud water
!
 7    continue
      ent=rmu*dz
      qcw=qcw/exp(ent)
      qcwa=qcw+dcw1-dcw2-dqcw3
      if (qcwa .lt. 0.) then
        red=(qcw+dcw1)/(dcw2+dqcw3)
        dcw2=dcw2*red
        dqcw3=dqcw3*red
      end if
      QCW=QCW+DCW1-DCW2-DQCW3

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 3e20.12)')  &
                             'in micro: ent,qcw,dcw1= ',ent,qcw,dcw1
      endif


      IF (QCW .LT. 0.) QCW=0.
      DQRW3=0.
      IF (QRW .EQ. 0.) GO TO 8
      DQRW3=(QRW**1.125)*5.1*DZ/(W    *HR)
!
!     calculate effect of entrainment on rain water
!
 8    continue
      qrw=qrw/exp(ent)
      qrwa=qrw+dcw2+dqcw3-dqrw3
      if (qrwa .lt. 0.) then
        red=(qrw+dcw2+dqcw3)/dqrw3
        dqrw3=red*dqrw3
      end if
      QRW=QRW+DCW2+DQCW3-DQRW3

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 2e20.12)') &
                               'in micro: ent,qrw= ',ent,qrw
      endif


      IF (QRW .LT. 0.) QRW=0.
      QLW=QRW+QCW

      if (debug_ijt) then
        write (unit_dc(n2), '(a, 3e20.12)')   &
                      'in micro: exit micro dcw1,dqrw3,qlw= ',   &
                       dcw1,dqrw3,qlw
      endif


!     QLW IN UNITS OF G/M**3
      TV=TCB*(1.+.61*((RS1   +RS2     )/2.))
      QEB=(QE1  +QE2    )/2.
      TEB=(TE1  +TE2    )/2.
      TVE=TEB*(1.+.61*QEB)
      QLW=2.0E-03*QLW*TV*RD/(P1+p2    )
      rmu=rmusa



!-------------------------------------------------------------------

end  subroutine micro




subroutine clotr (alp, ent, gravit, p1, p2, rd, sou1, sou2, tc1, tc2, &
                  te1, te2, xe1, xe2, xc, debug_ijt, n2)

!-----------------------------------------------------------------------
!     Calculates in-cloud tracer profile.
!
!     Leo Donner
!     GFDL
!     14 Jan 00
!
!-----------------------------------------------------------------------
!
!     IMPLICIT NONE
!
!
!
!---------------------Arguments-----------------------------------------
!
!     Input arguments
!
      real alp            !  virtual mass coefficient (dimensionless)
      real ent            !  entrainment coefficient (m**-1)
      real gravit         !  gravity constant (m/s2)
      real p1             !  pressure at level nearer earth surface (Pa)
      real p2             !  pressure at level farther from surface
                          !  (Pa)
      real rd             !  gas constant (J/kg K)
      real sou1           !  in-cloud source of x at pressure p1 (Pa)
      real sou2           !  in-cloud source of x at pressure p2 (Pa) 
      real tc1            !  cloud temperature (K) at pressure p1 (Pa)
      real tc2            !  cloud temperature (K) at pressure p2 (Pa)
      real te1            !  environmental temperature (K) at pressure  
                          !  p1 (Pa)
      real te2            !  environmental temperature (K) at pressure  
                          !  p2 (Pa)
      real xe1            !  large-scale tracer concentration at
                          !  pressure p1 (Pa)
      real xe2            !  large-scale tracer concentration at
                          !  pressure p2 (Pa)

      integer n2
      logical debug_ijt
!
!     Input/output arguments
!
      real xc             !  cloud-tracer concentration at pressure
                          !  p1 (Pa) on input, p2 (Pa) on output
!
!---------------------Local Workspace-----------------------------------
!
      real d1                  !  -dz/dp from p1 side  (m/Pa) 
      real d2                  !  -dz/dp from p2 side (m/Pa)
      real dz                  !  height increment (m)
      real epm                 !  defined in "Generalized mu, 10/1/89"
      real seb                 !  layer average of sou
      real xeb                 !  layer average of xe
!
!     NOTES:  xc and xe must have same units.
!             sou must have the (units of xc and xe)/sec.
!
!-----------------------------------------------------------------------
      integer ktest
      integer ktdiag
      integer idiag
      integer jdiag
      real tediag
!#include "cudiag.H"

      if (debug_ijt) then
        !print *, 'DONNER_DEEP/clotr: entering clotr'
        write (unit_dc(n2), '(a)') 'in clotr: entering clotr'
      endif


!
!     Initialization.
!
      epm=0.
!
!     Calculate in-cloud profile of tracer.
!
      d1=rd*(1.+alp)*tc1*te1/(gravit*p1*(alp*te1    &
         +tc1))
      d2=rd*(1.+alp)*tc2*te2/(gravit*p2*(alp*te2   &
         +tc2))
      dz=(d1+d2)*(p1-p2)/2.
      xeb=(xe1+xe2)/2.
      seb=(sou1+sou2)/2.
      xc=xc/exp(ent*dz)
      xc=xc+(((exp(ent*dz)-1.)*xeb)/    &
         exp(ent*dz))
      xc=xc+(seb*(1.+epm*(exp(ent*dz)-1.))   &
         /exp(ent*dz))
      if (xc .lt. 0.) xc=0.

      if (debug_ijt) then
        write (unit_dc(n2), '(a, e20.12)')   &
           'in clotr: xc= ',xc
        write (unit_dc(n2), '(a, 3e20.12)')   &
              'in clotr: d1,d2,dz= ',d1,d2,dz
        write (unit_dc(n2), '(a, e20.12)')   &
               'in clotr: xeb= ',xeb
        write (unit_dc(n2), '(a, e20.12)')   &
               'in clotr: seb= ',seb
        write (unit_dc(n2), '(a, e20.12)')   &
              'in clotr: ent= ',ent
      endif




end subroutine clotr




subroutine tae (t,p,q, lat, dp, cappa, ta)

real,     intent(in)   :: t, p, q, lat, dp, cappa
real,     intent(out)  :: ta

!
!     calculates adiabatic equivalent temperature
!
!     On Input:
!        t     temperature (K)
!        p     pressure (Pa)
!        q     specific humidity (kg(H2O)/kg)
!              various constants
!
!     On Output:
!        ta    adiabatic equivalent temperature
!
       integer :: k
       real    :: pstop, pr, te, es, qe

!
!     define constants
!
!     cappa=rair_mul/cpair_mul
!
      pstop = 4.0E03
      ta=t
!     do 1 i=1,nc
!     do 1 i=1,ncap
      do k=1,ncap
      pr=p+(k-1)*dp
      if (pr .le. 0.) exit      
!! CHANGE in connecting to fez
      if ( pr < pstop) exit       
      te=t*((pr/p)**cappa)
!     call establ(es,te)
!     if (te > 373.0 .or. te < 114.0) then
!       print *, 'mpp_pe, i, t, te, p, pr', mpp_pe(), i,   &
!                                       t, te, p, pr
!       go to 2
!     call escomp(te, es)
      call lookup_es(te, es)
      qe=epsilo*es/(p-es)
      if (q .ge. qe) then
!        ta=t*exp(latvap*q/(te*cpair_mul))
!        ta=t*exp(lat*q/(te*cpair_mul))
         ta=t*exp(lat*q/(te*CP_AIR))
!        go to 2
         exit
      end if
      end do
!1    continue
!2    continue

end subroutine tae





subroutine verav (qrnh,pb,pt,lcons,thetl,qrnl,qbl,ps,pr,phr,cappa,   &
                  debug_ijt,n2)


!--------------------------------------------------------------------
integer, intent(in) :: n2
logical,      intent(in) :: debug_ijt, thetl, lcons
real, dimension(:), intent(in) :: qrnh, pr, phr
real,         intent(in)  ::  qbl
real,         intent(in)  :: pb, pt,      ps, cappa
real, dimension(:), intent(inout) :: qrnl
!--------------------------------------------------------------------

!
!      layer averaging for large-scale       source due to cumulus
!      convection
!      Verav notes 1/7/04 (available from Leo Donner) explain
!      this routine in more detail, especially the procedures
!      use to enforce tracer conservation, which are invoked when
!      lcons=.true.
!
!      on input:
!        lcons     if .true., vertical integral of qrnl obeys
!                  conservation constraint. integral is
!                  zero if thetl is .false. integral of
!                  qrnl x (ratio of potential temperature to
!                  temperature) is zero if thetl is .true.
!                  no constraint on vertical integral of qrnl
!                  if lcons is .false.
!        qrnh(ncm) large-scale       source at cloud-model resolution
!        qbl  large-scale PBL    source due to Cu
!        pb        cloud-base pressure (Pa)
!        pt        cloud-top pressure (Pa)
!        ps        surface pressure (Pa)
!        thetl     indicator for temperature calculation
!        pr        large-scale pressure levels (Pa)
!        phr       large-scale pressure half levels (Pa)
!        cappa     constant
!
!     on output:
!
!        qrnl(nlev)large-scale       source at CCM resolution
!

real, dimension (ncap) :: p
logical                :: test2
real                   :: qblimp  !  BL tendency required
                                  !  for convervation in absence of
                                  !  of BL source


      integer :: i, j1, j
      real :: dp, ptt,rintsum, ph, pl, rint, rkou, phrh, phrl, qlsum,&
              qlsu, pi, qbl0, sb, wta, wtb, thetsum






      test2=.false.

      if (debug_ijt) then
        write (unit_dc(n2), '(a)') 'in verav: entering verav'
      endif


      dp=-1000.
      ptt=pt+dp
      do 1 i=1,ncap
        p(i)=pb+(i-1)*dp
 1    continue
      do 8 i=1,nlev
        qrnl(i)=0.
  8   continue
      j1=1
      rintsum=0.
      do 2 i=1,nlev
        ph=phr(i+1)
        pl=phr(i)

        if (debug_ijt) then
          write (unit_dc(n2), '(a, 4e20.12)')  &
                           'in verav: pl,ph,pb,pt',pl,ph,pb,pt
        endif


        if (ph .ge. pb) go to 30

        if (debug_ijt) then
          write (unit_dc(n2), '(a, l4)') 'in verav: thetl= ',thetl
        endif


        if (pl .le. ptt) go to 7
        rint=0.0
        rkou=0.
         do 4 j=j1,ncap-1
         phrh=(p(j)+p(j+1))/2.
         phrl= p(j)-(dp/2.)
         if (j .eq. 1 .and. (pl .ge. pb)) rkou=pl-pb
         if (phrl .gt. pl) phrl=pl
         if (phrl .gt. pb) phrl=pb
         if (phrh .lt. ph) phrh=ph
         if (phrh .le. ptt) then
            rkou=rkou+phrl-ph
            rint=rint+qrnh(j)*(phrl-ptt)
            qrnl(i)=rint/rkou
            rintsum=rintsum+qrnh(j)*(phrl-ptt)
            go to 7
         end if

         if (debug_ijt) then
           write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                     'in verav: j,phrl,phrh= ',j,phrl,phrh
         endif


         rkou=rkou+phrl-phrh
         rint=rint+qrnh(j)*(phrl-phrh)
         rintsum=rintsum+qrnh(j)*(phrl-phrh)

         if (debug_ijt) then
           write (unit_dc(n2), '(a, e20.12)')  &
                             'in verav: rintsum= ',rintsum
         endif

! This routine assumes that the thicknesses of all layers in
! the model calling the cumulus parameterization are greater
! than that of the cloud model in the cumulus paramterization.
! If not, uncomment the following lines and seriously consider
! whether this approach is physically appropriate under such
! conditions.
!    if ((phr(i)-phr(i+1)) .lt. (-dp)) j1=j+1

         j1=j

         if (debug_ijt) then
!          print *, 'DONNER_DEEP/verav: j,qrnh,p,rkou',j, &
!                           qrnh(j),p(j),rkou
           write (unit_dc(n2), '(a, i4, 3e20.12)')  &
                 'in verav: j,qrnh,p,rkou',j,qrnh(j),p(j),rkou
!          print *, 'DONNER_DEEP/verav: rint= ',rint
           write (unit_dc(n2), '(a, e20.12)') 'in verav: rint= ',rint
         endif


         if (phrh   .le. ph       ) then
            rint=rint/rkou
            qrnl(i)=rint
            go to 30
         end if
 4       continue
 30   continue

      if (debug_ijt) then
        write (unit_dc(n2), '(a, i4, 3e20.12)') &
                      'in  verav: i,pl,ph,qrnl(i)',i,pl,ph,qrnl(i)
      endif


 2    continue
 7    continue
      if (lcons) then
        qlsum=0.
        qlsu=0.
        do 9 i=1,nlev

          if (debug_ijt) then
            write (unit_dc(n2), '(a, i4, f19.10, l4)')  &
                            'in verav: i,pr,thetl= ',i,pr(i),thetl
          endif


          if (thetl) then
          pi=(1.0e05/pr(i) )**cappa
          qlsum=qlsum+qrnl(i)*pi*(phr(i)-phr(i+1))
          if (phr(i) .le. pb)    &
                    qlsu=qlsu+qrnl(i)*pi*(phr(i)-phr(i+1))
          end if  ! (thetl)
          if (.not. thetl) qlsum=qlsum+qrnl(i)*(phr(i)-phr(i+1))
 9      continue
        qbl0=qbl
        qblimp=-qlsum/(ps-pb)

        if (debug_ijt) then
          write (unit_dc(n2), '(a, e20.12)') 'in verav: qblimp= ',qblimp
          write (unit_dc(n2), '(a, e20.12)') 'in verav: qlsu=',qlsu
          write (unit_dc(n2), '(a, e20.12)')  &
                                    'in verav: thetl qlsum= ',qlsum
        endif ! (debug_ijt)



      if (debug_ijt) then
        write (unit_dc(n2), '(a, e20.12)') 'in verav: ps= ',ps
      endif


      do 3 i=1,nlev

        if (debug_ijt) then
          write (unit_dc(n2), '(a, i4, 3e20.12)')  &
                       'in verav: i,qblimp,phr,pr= ',i,qblimp,  &
                                                     phr(i+1),pr(i)
        endif


        if (phr(i+1) .ge. pb) then
          qrnl(i)=qblimp
          if (thetl) then
            pi=(1.0e05/pr(i))**cappa
            qrnl(i)=qblimp/pi
          end if ! (thetl)
          if (.not. thetl) qrnl(i)=qblimp
         end if ! (phr(i+1) .ge. pb)
         if (phr(i+1) .lt. pb) then

           if (debug_ijt) then
             write (unit_dc(n2), '(a, 3e20.12)')  &
                   'in verav: phr,phr+,qrnl= ',phr(i),phr(i+1),qrnl(i)
           endif


           if (phr(i) .le. pb) go to 5
             wta=qrnl(i)*(phr(i)-phr(i+1))
             wtb=qblimp*(phr(i)-pb)
             qrnl(i)=(wta+wtb)/(phr(i)-phr(i+1))
             if (thetl) then
               pi=(1.0e05/pr(i))**cappa
               thetsum=qlsum*(ps-phr(i))/(pb-ps)

               if (debug_ijt) then
                 write (unit_dc(n2), '(a, 4e20.12)')  &
                   'in verav: qbl0,qlsum,qlsu,thetsum= ',qbl0,qlsum,  &
                                                      qlsu,thetsum
               endif


               qrnl(i)=(qlsu)+thetsum
               qrnl(i)=-qrnl(i)
               qrnl(i)=qrnl(i)/(pi*(phr(i)-phr(i+1)))
               qrnl(i)=qrnl(i) + ( (qbl0/  (phr(i)-phr(i+1)))   &
                       *(phr(i)-pb ))
             end if
           end if
 3    continue
 5    continue
      qlsum=0.
      do 11 i=1,nlev

        if (debug_ijt) then
          write (unit_dc(n2), '(a, i4, 2e20.12)') &
                          'in verav: i,phr,phr+= ',i,phr(i),phr(i+1) 
          write (unit_dc(n2), '(a, f19.10, a, e20.12)')  &
                           'in verav: pr= ',pr(i),'qrnl= ',qrnl(i)
        endif


        pi=(1.0e05/pr(i))**cappa
        if (.not. thetl) qlsum=qlsum+qrnl(i)*(phr(i)-phr(i+1))
        if (thetl) qlsum=qlsum+qrnl(i)*pi*(phr(i)-phr(i+1))
 11   continue

      if (debug_ijt) then
        write (unit_dc(n2), '(a, e20.12)') 'in verav: qlsum= ',qlsum
      endif
     end if ! (lcons)


end subroutine verav



subroutine ver (xav, p1, p2, pr, phr, x, debug_ijt, n2)

!------------------------------------------------------------------
 integer, intent(in) :: n2
 logical,   intent(in) :: debug_ijt
 real,     intent(in) :: xav,p1, p2
  real, dimension(:), intent(in) :: phr
!  real, dimension(:), intent(in) :: pr
!!! COMPILER BUG ?????
!!! pr must be dimensioned below to work properly ???  
 real, dimension(:), intent(out) :: x
!------------------------------------------------------------------

!
!     vertical averaging subroutine for mesoscale tendencies
!
!     On Input:
!       xav     average tendency
!       p1      high pressure boundary(Pa)
!       p2      low pressure boundary(Pa)
!       pr,phr  low-resolution pressure levels
!
!     On Output:
!       x       vertically averaged tendency (g/kg/day)
!

!      dimension         pr(nlev)               




!******************************************************************


!---- NOTICE:
!
! THE FOLLOWING LINE IS NEEDED, BUT NEED NOT BE CHANGED IF VERTICAL|
! RESOLUTION IS CHANGED
! IT LIKELY REFLECTS THE PRESENCE OF A BUG SOMEWHERE WHICH HAS YET TO
! BE LOCATED
!
!
       real pr
       dimension         pr(127 )               
!
       integer :: i
!
!
!---- END OF NOTICE



!******************************************************************


!      dimension         pr( size(phr)-1 )               

    

      if (debug_ijt) then
        write (unit_dc(n2), '(a, e20.12, 2f19.10)')  &
                               'in ver: xav,p1,p2= ',xav,p1,p2
      endif

      if ( p1 < p2 ) then
       call error_mesg ('DONNER_DEEP/ver',  &
           ' input pressure p1 is less than input pressure p2', FATAL)
      endif



!
      do 1 i=1,nlev
        if (p1 .le. phr(i+1)) x(i)=0.
        if (p2 .ge. phr(i)) x(i)=0.
        if ((p1 .ge. phr(i)) .and. (p2 .le. phr(i+1)) )   &
                                                      x(i)=xav
        if ( (p1 .le. phr(i)) .and. (p1 .ge. phr(i+1))     &
             .and. (p2 .le. phr(i+1)) ) x(i)=xav*(p1-phr(i+1))/  &
                                             (phr(i)-phr(i+1))
        if ( (p2 .ge. phr(i+1)) .and. (p2 .le. phr(i))   &
              .and. (p1 .ge. phr(i)) ) x(i)=xav*(phr(i)-p2)/  &
                                            (phr(i)-phr(i+1))
        if ((p1 .le. phr(i)) .and. (p2 .ge. phr(i+1))) x(i)=   &
                                          xav*(p1-p2)/(phr(i)-phr(i+1))
 1    continue

      if (debug_ijt) then
        do i=1,nlev
          if (x(i) /= 0.0) then
            write (unit_dc(n2), '(a, i4, e20.12)') &
                                           'in ver: i,x= ',i,x(i)
          endif
        end do
      endif



end subroutine ver




subroutine polat (xv, pv,x,p, debug_ijt, n2)

!---------------------------------------------------------------------
integer, intent(in) :: n2
logical,            intent(in) :: debug_ijt
real, dimension(:), intent(in) :: xv, pv
real,              intent(in)  :: p
real,               intent(inout)  :: x
!---------------------------------------------------------------------



!
!     INTERPOLATES XV TO PRESSURE P.
!
!     ON INPUT:
!
!         XV(N)   DATA AT RESOLUTION N
!         PV(N)   PRESSURE AT N LEVELS
!
!     ON OUTPUT:
!
!         X       DATA AT PRESSURE P
!


       integer :: i

       if (debug_ijt) then
          write (unit_dc(n2), '(a, f19.10)') 'in polat: p=', p
       endif


       IF (P .GE. PV(1)) THEN

         if (debug_ijt) then
           write (unit_dc(n2), '(a, f19.10)') 'in polat: p= ', p
         endif


          X=(XV(2)-XV(1))/(PV(2)-PV(1))
          X=X*(P-PV(1))+XV(1)
          GO TO 1
       ENDIF
       IF (P .LE. PV(nlev)) THEN

         if (debug_ijt) then
           write (unit_dc(n2), '(a, f19.10)') 'in polat: p= ', p
         endif


          X=(XV(nlev)-XV(nlev-1))/(PV(nlev)-PV(nlev-1))
          X=X*(P-PV(nlev))+XV(nlev)
          GO TO 1
       ENDIF
       DO 10 I=1,nlev-1
       IF ((P .GE. PV(i+1)) .AND. (P .LE. PV(I))) THEN

         if (debug_ijt) then
            write (unit_dc(n2), '(a, f19.10)')  'in polat: p= ', p
!           print *, 'DONNER_DEEP/polat: i,xv(i),xv(i+1)= ',i, &
!                                                  xv(i),xv(i+1)
            write (unit_dc(n2), '(a, i4, 2e20.12)')  &
                         'in polat: i,xv(i),xv(i+1)= ',i,xv(i),xv(i+1)
         endif


          X=(XV(I+1)-XV(I))/(PV(I+1)-PV(I))
          X=X*(P-PV(I+1))+XV(I+1)
          GO TO 1
       ENDIF
 10    CONTINUE
 1     CONTINUE

end subroutine polat


!####################################################################

!######################################################################

subroutine donner_deep_sum( is, js, Don_conv) 
!------------------------------------------------------------------
! this subroutine puts into global storage the fields needed by the
! radiation scheme, after cleaning them up a little bit by checks
! put in by D. Schwarzkopf and L. Donner.
!------------------------------------------------------------------

integer, intent(in)                   :: is,js
type(donner_conv_type), intent(inout) :: Don_conv
                                                            

!------------------------------------------------------------------
! declare local variables


integer :: k, idim, jdim, kdim, unit
integer :: ie, je

real, dimension(size(Don_conv%xice,1), &
                size(Don_conv%xice,2),  &
                size(Don_conv%xice,3)) :: &
                                               meso_area
   
!------------------------------------------------------------------
! code
  

!-------------------------------------------------------------------
! compute useful integers, the dimensions of the array (idim,jdim,kdim)
! and the ending integers

  idim = SIZE(Don_conv%xice,1)
  jdim = SIZE(Don_conv%xice,2)
  kdim = SIZE(Don_conv%xice,3)

  ie = is + idim - 1
  je = js + jdim - 1


!-------------------------------------
! increment counter

  if (do_average) then
       nsum(is:ie,js:je)   =  nsum(is:ie,js:je)   +  1
  else
       nsum(is:ie,js:je)   =  1
  end if

!      unit = open_file ('fort.149', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' donner_deep_1'
! idim = size(nsum,1)
! jdim = size(nsum,2)
! write (unit,*) idim,jdim,is,ie,js,je
! write (unit,*) nsum
!      call close_file (unit)

!================================
!
!   DEFINE MESOSCALE PROPERTIES
!
!================================
  

!-------------------------------------------
! create 3d meso anvil area 
!
 
  do k=1,nlev  
     where (Don_conv%xice(:,:,k) == 0.0)
            meso_area(:,:,k) =  0.0
     elsewhere
            meso_area(:,:,k) =  max(0.0,Don_conv%ampta1(:,:))
     endwhere
  enddo  

!------------------------------------
! define mesoscale anvil area
!
    
  if (do_average) then
     meso_cloud_frac(is:ie,js:je,:) = &
     meso_cloud_frac(is:ie,js:je,:) + meso_area(:,:,:)
  else
     meso_cloud_frac(is:ie,js:je,:) = meso_area(:,:,:)
  end if

!------------------------------------
! define mesoscale liquid amt
!

  if (do_average) then 
       meso_liquid_amt(is:ie,js:je,:) = &
       meso_liquid_amt(is:ie,js:je,:) + 0.0
  else
       meso_liquid_amt(is:ie,js:je,:) = 0.0
  end if

!------------------------------------
! define mesoscale liquid size
!
 
  if (do_average) then 
        meso_liquid_size(is:ie,js:je,:) = &
                              meso_liquid_size(is:ie,js:je,:) + 0.0
  else
        meso_liquid_size(is:ie,js:je,:) = 0.0
  end if
  
!------------------------------------
! define mesoscale ice amt
!
 
  if (do_average) then 
        meso_ice_amt(is:ie,js:je,:) = &
                  meso_ice_amt(is:ie,js:je,:) + Don_conv%xice(:,:,:)
  else
        meso_ice_amt(is:ie,js:je,:) = Don_conv%xice(:,:,:)
  end if

!------------------------------------
! define mesoscale ice size
!
 
  if (do_average) then 
        meso_ice_size(is:ie,js:je,:) = &
                  meso_ice_size(is:ie,js:je,:) + Don_conv%dgeice(:,:,:)
  else
        meso_ice_size(is:ie,js:je,:) = Don_conv%dgeice(:,:,:)
  end if
 

!================================
!
!   DEFINE CELL PROPERTIES
!
!================================
  
!------------------------------------
! define cell area
!
 
  if (do_average) then
        cell_cloud_frac(is:ie,js:je,:) = &
                           cell_cloud_frac(is:ie,js:je,:) + &
!                          max( 0.0, cual_3d(:,:,:) - meso_area(:,:,:) )
                     max( 0.0, Don_conv%cual(:,:,:) - meso_area(:,:,:) )
  else
        cell_cloud_frac(is:ie,js:je,:) = &
!                        max( 0.0, cual_3d(:,:,:) - meso_area(:,:,:) )
                     max( 0.0, Don_conv%cual(:,:,:) - meso_area(:,:,:) )
  end if
  
!------------------------------------
! define cell liquid amt
!
  if (do_average) then 
        cell_liquid_amt(is:ie,js:je,:) = &
                   cell_liquid_amt(is:ie,js:je,:) + Don_conv%cuql(:,:,:)
  else
        cell_liquid_amt(is:ie,js:je,:) = Don_conv%cuql(:,:,:)
  end if
  
!------------------------------------
! define cell liquid size
!
 
  if (do_input_cell_liquid_size) then
  
  
    if (do_average) then  
       cell_liquid_size(is:ie,js:je,:) = &
       cell_liquid_size(is:ie,js:je,:) + cell_liquid_eff_diam_input
    else  
       cell_liquid_size(is:ie,js:je,:) = cell_liquid_eff_diam_input
    end if
    
  else if (do_bower_cell_liquid_size) then
  
    if (do_average) then 
       cell_liquid_size(is:ie,js:je,:) = &
       cell_liquid_size(is:ie,js:je,:) + Don_conv%cell_liquid_eff_diam
    else 
       cell_liquid_size(is:ie,js:je,:) = Don_conv%cell_liquid_eff_diam
    end if

  endif

!------------------------------------
! define cell ice amt
!
 
  if (do_average) then
     cell_ice_amt(is:ie,js:je,:) = &
     cell_ice_amt(is:ie,js:je,:) + Don_conv%cuqi(:,:,:)
  else
     cell_ice_amt(is:ie,js:je,:) = Don_conv%cuqi(:,:,:)
  end if
  
!------------------------------------
! define cell ice size
!
  
  if (do_default_cell_ice_size) then

     if (do_average) then
       cell_ice_size(is:ie,js:je,:) = &
       cell_ice_size(is:ie,js:je,:) + cell_ice_geneff_diam_def
     else
       cell_ice_size(is:ie,js:je,:) = cell_ice_geneff_diam_def
     end if
      
  else if (do_input_cell_ice_size) then
     
     if (do_average) then       
       cell_ice_size(is:ie,js:je,:) = &
       cell_ice_size(is:ie,js:je,:) + cell_ice_geneff_diam_input
     else
       cell_ice_size(is:ie,js:je,:) = cell_ice_geneff_diam_input
     end if
     
  endif


!--------------------------------------------

end subroutine donner_deep_sum





subroutine cell_liquid_size_comp(pfull, temp, Don_conv, land)

!---------------------------------------------------
! This subroutine calculates the effective radii of
! liquid cloud drops in cumulus clouds following
! the prescription of Bower et al.
!---------------------------------------------------

real, dimension(:,:,:), intent(in)       :: pfull, temp  
type(donner_conv_type), intent(inout)    :: Don_conv
real, dimension(:,:)  , intent(in)       :: land


   
!------------------------------------------------------------------
! local variables
  
integer      :: i, j, k, idim, jdim, kdim

real, dimension (size(pfull,1),size(pfull,2)) ::   &
                                      cell_pbase, temp_cell_pbase, &
                                cell_land_ref_delp, cell_ocean_ref_delp
                                       
real, dimension (size(pfull,1),size(pfull,2),size(pfull,3)) ::   &
                          cell_liquid_eff_diam_land, cell_delp, &
                          cell_liquid_eff_diam_ocean
     integer :: unit

!------------------------
! code

      idim = SIZE(Don_conv%cuql,1)
      jdim = SIZE(Don_conv%cuql,2)
      kdim = SIZE(Don_conv%cuql,3)

!       compute liquid effective diameter using formulation of
!       Bower et al (JAS, 1994)
!
!
!   stage 1:       compute the pressure (cell_pbase) of the base of
!              cell clouds and the difference (cell_delp) between
!              (cell_pbase) and the pressure at full levels above
!              cell_pbase)( so long as a cell cloud exists).
!                   also compute the pressure difference between
!              (cell_pbase) and a point 500 and 1500 meters above 
!               the cloud base pressure (cell_land_delp and 
!              cell_ocean_delp, respectively). obtain the weighted
!   stage 2  :      compute separately for the land and ocean portion
!              of the grid box:
!              if (cell_delp) is greater than (cell_land_delp(say))
!              use specified effective radius. otherwise derive 
!              effective radius using formula.
!   stage 3  : the effective diameter is twice the weighted (by land/
!              ocean fraction) effective radius.
!
!   stage 1

       cell_pbase = pfull(:,:,1)
       temp_cell_pbase = temp  (:,:,1)

       do i=1,idim
       do j=1,jdim
       do k=nlev,1,-1
         if (Don_conv%cuql(i,j,k) >= 1.0e-11 )  then
           cell_pbase(i,j) = pfull(i,j,k)
           temp_cell_pbase(i,j) = temp  (i,j,k)
           exit
         end if
       end do
       end do
       end do

       do k=1,nlev
!         cell_delp(:,:,k) = pfull(:,:,k) - cell_pbase(:,:)
          cell_delp(:,:,k) = cell_pbase(:,:) - pfull(:,:,k)
       enddo

!      cell_land_ref_delp = cell_pbase*                     &
!              EXP( -(delz_land*gravm/(rair*temp_cell_pbase)) ) 
       cell_land_ref_delp = cell_pbase*(1.0 -               &
               EXP( -(delz_land*gravm/(rair*temp_cell_pbase)) ))
!      cell_ocean_ref_delp = cell_pbase*                     &
!              EXP( -(delz_ocean*gravm/(rair*temp_cell_pbase)) )
       cell_ocean_ref_delp = cell_pbase*(1.0 -               &
               EXP( -(delz_ocean*gravm/(rair*temp_cell_pbase)) ))

!      unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' donner_deep,cell_liquid_size_comp,stage 1'
!        write (unit,*) idim,jdim,kdim
!        write (unit,*) ' Don_conv%cuql'
! write (unit,*) Don_conv%cuql
!        write (unit,*) ' cell_land_ref_delp'
! write (unit,*) cell_land_ref_delp
!        write (unit,*) ' cell_ocean_ref_delp'
! write (unit,*) cell_ocean_ref_delp
!        write (unit,*) ' land'
! write (unit,*) land
!        write (unit,*) ' cell_delp'
! write (unit,*) cell_delp
!        write (unit,*) ' cell_pbase'
! write (unit,*) cell_pbase
!        write (unit,*) ' pfull'
! write (unit,*) pfull
!      call close_file (unit)
!  stage 2
! set "default" land, ocean diameters to zero. these will be
! overwritten whenever the cell liquid conc (Don_conv%cuql)  >= 1.0E-11

       cell_liquid_eff_diam_land = 0.0
       cell_liquid_eff_diam_ocean = 0.0

       do k = 1,nlev
       do j=1,jdim
       do i=1,idim
       if (Don_conv%cuql(i,j,k) >= 1.0e-11) then
!      unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' Don_conv%cuql large enough', i,j,k
!      call close_file (unit)
         if (land(i,j) > 0.0) then
!   do land calculation only if land > 0
           if (cell_delp(i,j,k) >= cell_land_ref_delp(i,j)) then 
!      unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' land,cell_delp large enof',r_conv_land
!      call close_file (unit)
            cell_liquid_eff_diam_land(i,j,k) = 2.0*r_conv_land
           else
!      unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' land,cell_delp  too small'
!      call close_file (unit)
    cell_liquid_eff_diam_land(i,j,k) =  2.0 *  (1.0e6) *      &
     (3.0*(pfull(i,j,k)/(rair*temp  (i,j,k)))*Don_conv%cuql(i,j,k) / &
!      (4*pie*rho_water*N_land) ) ** (1./3.)      
      (4*pie*DENS_H2O *N_land) ) ** (1./3.)      
           endif
         endif
         if (land(i,j) < 1.0) then
!   do ocean calculation only if land < 1.0
           if (cell_delp(i,j,k) >= cell_ocean_ref_delp(i,j)) then 
            cell_liquid_eff_diam_ocean(i,j,k) = 2.0*r_conv_ocean
           else
             cell_liquid_eff_diam_ocean(i,j,k) =  2.0 *  (1.0e6) *     &
       (3.0*(pfull(i,j,k)/(rair*temp  (i,j,k)))*Don_conv%cuql(i,j,k) / &
!      (4*pie*rho_water*N_ocean) ) ** (1./3.)    
          (4*pie*DENS_H2O *N_ocean) ) ** (1./3.)    
           endif
         endif
       endif
      enddo
      enddo
      enddo


!  stage 3
       do k = 1,nlev
       do i=1,idim
       do j=1,jdim
         Don_conv%cell_liquid_eff_diam(i,j,k) =                    &
              land(i,j)*        cell_liquid_eff_diam_land(i,j,k)   + &
              (1.0 - land(i,j))*cell_liquid_eff_diam_ocean(i,j,k) 
       enddo
       enddo
       enddo

!      Don_conv%cell_liquid_eff_diam = 10.0
       where (Don_conv%cuql == 0.0)
         Don_conv%cell_liquid_eff_diam = 10.0
       end where
       where (Don_conv%cell_liquid_eff_diam < 8.401E+00)
         Don_conv%cell_liquid_eff_diam = 8.401E+00
       end where
       where (Don_conv%cell_liquid_eff_diam > 33.199E+00)
         Don_conv%cell_liquid_eff_diam = 33.199E+00
       end where
       
!      unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' donner_deep,cell_liquid_size_comp'
!        write (unit,*) ' Don_conv%cell_liquid_eff_diam'
! write (unit,*) Don_conv%cell_liquid_eff_diam
!        write (unit,*) ' cell_liquid_eff_diam_land'
! write (unit,*) cell_liquid_eff_diam_land
!        write (unit,*) ' cell_liquid_eff_diam_ocean'
! write (unit,*) cell_liquid_eff_diam_ocean
!      call close_file (unit)
end subroutine cell_liquid_size_comp

!#######################################################################

 subroutine strat_cloud_donner_tend (is, ie, js, je, Dmeso, qlmeso,  &
                                     dt, qimeso, Mtot, phalf, ql,  &
                                     qi, cf)

!-----------------------------------------------------------------------
! input
!
!               vertical index 1 at model top
! Dmeso  mass detrainment rate from mesoscale region to large-scale
! region (sec-1)
! qlmeso cloud liquid specific humidity (kg condensate/kg air)
! qimeso cloud ice specific humidity (kg condensate/kg air)
! Mtot   total mass flux = mesoscale_mass_flux + convective_mass_flux
!               (kg /m2/sec) defined on level interfaces
!
!               NOTE: Regardless of what they contain, Mtot(:,:,1)
!                     Mtot(:,:,kdim+1) will be assumed to be zero.
!
! phalf pressure on model interfaces (Pa)
! ql large-scale cloud liquid specific humidity
! qi large-scale cloud ice specific humidity
! cf large-scale cloud fraction (0-1)
!
! output
!
! delta_ql large-scale cloud liquid increment (kg cond/kg air)
! delta_qi large-scale cloud ice increment (kg cond/kg air)
! delta_qa large-scale cloud increment tendency 
!
!-----------------------------------------------------------------------
   integer, intent(in) :: is, ie, js, je
   real, intent(in),  dimension(:,:,:) :: Dmeso, qlmeso, qimeso
   real, intent(in) :: dt
   real, intent(in),  dimension(:,:,:) :: ql, qi, cf
   real, intent(in),  dimension(:,:,:) :: Mtot,phalf
!-----------------------------------------------------------------------
   integer kdim
   integer unit
   real, dimension(size(Dmeso,1),size(Dmeso,2),size(Dmeso,3)) :: mass    ! kg
   real, dimension(size(cf,1),size(cf,2),size(cf,3)) :: cffix
   real :: cffmax
!  air/m2 of each level
   integer k
   integer kk
!-----------------------------------------------------------------------

   kdim = size(Dmeso,3)
   mass(:,:,1:kdim) = (phalf(:,:,2:kdim+1)-phalf(:,:,1:kdim))/gravm

    delta_ql (is:ie,js:je,1)=0.
    delta_qi (is:ie,js:je,1)=0.
    delta_qa (is:ie,js:je,1)=0.
!   qltend=0.
!   qitend=0.
!   cftend=0.
    
!       unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!         write (unit,*) ' donner_deep,strat_cloud_donner_tend_1'
!         write (unit,*) ' qltend'
!         write (unit,*) qltend
!         write (unit,*) ' qitend'
!         write (unit,*) qitend
!       call close_file (unit)
   !do incoming compensating subsidence fluxes from above
   delta_ql (is:ie,js:je,2:kdim)= Mtot(:,:,2:kdim) * &
                  0.5*(ql(:,:,1:kdim-1)+ql(:,:,2:kdim))/mass(:,:,2:kdim)
!   qitend(:,:,2:kdim)=qitend(:,:,2:kdim) +  Mtot(:,:,2:kdim) * &
   delta_qi (is:ie,js:je,2:kdim)= Mtot(:,:,2:kdim) * &
                  0.5*(qi(:,:,1:kdim-1)+qi(:,:,2:kdim))/mass(:,:,2:kdim)
!  cftend(:,:,2:kdim)=cftend(:,:,2:kdim) +  Mtot(:,:,2:kdim) * &
   delta_qa (is:ie,js:je,2:kdim)= Mtot(:,:,2:kdim) * &
                  0.5*(cf(:,:,1:kdim-1)+cf(:,:,2:kdim))/mass(:,:,2:kdim)
   !do outgoing compensating subsidence fluxes out the bottom
!  qltend(:,:,1:kdim-1)=qltend(:,:,1:kdim-1) -  Mtot(:,:,2:kdim) * &
   delta_ql (is:ie,js:je,1:kdim-1)=delta_ql(is:ie,js:je,1:kdim-1) -  &
                                   Mtot(:,:,2:kdim) * &
                                   0.5*(ql(:,:,1:kdim-1)+  &
                                   ql(:,:,2:kdim))/mass(:,:,1:kdim-1)
   delta_qi (is:ie,js:je,1:kdim-1)=delta_qi (is:ie,js:je,1:kdim-1) - &
                                    Mtot(:,:,2:kdim) * &
                                    0.5*(qi(:,:,1:kdim-1)+ &
                                    qi(:,:,2:kdim))/mass(:,:,1:kdim-1)
!  cftend(:,:,1:kdim-1)=cftend(:,:,1:kdim-1) -  Mtot(:,:,2:kdim) * &
   delta_qa( is:ie,js:je,1:kdim-1)=delta_qa (is:ie,js:je,1:kdim-1) - &
                                   Mtot(:,:,2:kdim) * &
                                   0.5*(cf(:,:,1:kdim-1)+ &
                                   cf(:,:,2:kdim))/mass(:,:,1:kdim-1)
   !do detrainment from meso region
!  qltend(:,:,:) = qltend(:,:,:) + Dmeso(:,:,:)*qlmeso(:,:,:)
!  qitend(:,:,:) = qitend(:,:,:) + Dmeso(:,:,:)*qimeso(:,:,:)
   delta_ql (is:ie,js:je,:) = delta_ql (is:ie,js:je,:) +  &
                              Dmeso(:,:,:)*qlmeso(:,:,:)
   delta_qi (is:ie,js:je,:) = delta_qi (is:ie,js:je,:) +  &
                              Dmeso(:,:,:)*qimeso(:,:,:)
   where ((qlmeso+qimeso) .ge. 1.e-10)
!  cftend(:,:,:) = cftend(:,:,:) + Dmeso(:,:,:)
     delta_qa (is:ie,js:je,:) = delta_qa (is:ie,js:je,:) + Dmeso(:,:,:)
   end where
      
!       unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!         write (unit,*) ' donner_deep,strat_cloud_donner_tend_2'
!         write (unit,*) ' ql'
!         write (unit,*) ql
!         write (unit,*) ' qi'
!         write (unit,*) qi
!         write (unit,*) ' cf'
!         write (unit,*) cf
!         write (unit,*) ' qltend'
!         write (unit,*) qltend
!         write (unit,*) ' qitend'
!         write (unit,*) qitend
!       call close_file (unit)
      where (cf > 1.00 .or. cf < 0.0)
         cffix = 100.
      else where
         cffix = 1.0
      end where
      cffmax = maxval(cffix)
!if (cffmax > 5.0) then
!       unit = open_file ('fort.152', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
! write (unit,*) ' cf out of range, cffmax = ', cffmax
!       call close_file (unit)
!       endif
!-----------------------------------------------------------------------

!    qltend = qltend*dt
!    qitend = qitend*dt
!    cftend = cftend*dt
     delta_ql (is:ie,js:je,:) = delta_ql(is:ie,js:je,:)*dt
     delta_qi (is:ie,js:je,:) = delta_qi (is:ie,js:je,:)*dt
     delta_qa (is:ie,js:je,:) = delta_qa (is:ie,js:je,:)*dt
!
!    mulsub allowed ice and liquid from convective system to evaporate
!    and/or sublimate as part of thermal and moisture forcing terms
!    remove those tendencies here. different assumptions used to
!    calculate these increments/tendencies here and in mulsub, so
!    some residual phase change will generally remain
!
     cememf(is:ie,js:je,:) = cememf(is:ie,js:je,:) -  &
     (delta_ql(is:ie,js:je,:)/dt) - (delta_qi(is:ie,js:je,:)/dt) 
     cemetf(is:ie,js:je,:) = cemetf(is:ie,js:je,:) + &
     (delta_ql(is:ie,js:je,:)*latvap/(cp_air*dt)) + &
     (delta_qi(is:ie,js:je,:)*(latvap+latice)/(cp_air*dt))


end subroutine strat_cloud_donner_tend

!##################################################################




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!          4. COLUMN DIAGNOSTICS SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!##################################################################

subroutine donner_column_init (pref, Time)

real, dimension(:), intent(in) :: pref
type(time_type), intent(in) :: Time


      integer :: n, k


        column_diagnostics_desired = .true.

!---------------------------------------------------------------------
!    check that array dimensions are sufficiently large for the number 
!    of columns requests.
!---------------------------------------------------------------------
        if (num_diag_pts > MAX_PTS) then
          call error_mesg ('donner_deep_mod', &
         'must reset MAX_PTS or reduce number of diagnostic points', &
                                                           FATAL)  
        endif

!---------------------------------------------------------------------
!    check that the specified time at which diagnostics are to be 
!    activated has been specified and is after the current time.
!---------------------------------------------------------------------
        do n=1,3
          if (diagnostics_start_time(n) == 0) then
            call error_mesg ('donner_deep_mod', &
             'year, month and/or day invalidly specified for column '//&
                  'diagnostics starting time', FATAL)
          endif
        end do
        Time_col_diagnostics = set_date (diagnostics_start_time(1), &
                                         diagnostics_start_time(2), &   
                                         diagnostics_start_time(3), &   
                                         diagnostics_start_time(4), &   
                                         diagnostics_start_time(5), &   
                                         diagnostics_start_time(6) )    
        if (Time_col_diagnostics < Time) then
          call error_mesg ( 'donner_deep_mod', &
            'specified time for column diagnostics is prior to start'//&
                                              ' of run', FATAL)
        endif

!---------------------------------------------------------------------
!    allocate space for the arrays used for column diagnostics.
!---------------------------------------------------------------------
        allocate (do_column_diagnostics   (jdf) )
        allocate (col_diag_unit             (num_diag_pts) )
        allocate (col_diag_lon                (num_diag_pts) )
        allocate (col_diag_lat                (num_diag_pts) )
        allocate (col_diag_i                  (num_diag_pts) )
        allocate (col_diag_j                  (num_diag_pts) )

!---------------------------------------------------------------------
!    call initialize_diagnostic_columns to determine the locations 
!    (i,j,lat and lon) of any diagnostic columns in this processor's
!    space and to open output files for the diagnostics.
!---------------------------------------------------------------------
        call initialize_diagnostic_columns   &
                     (mod_name, num_diag_pts_latlon, num_diag_pts_ij, &
                      i_coords_gl, j_coords_gl, lat_coords_gl, &
                      lon_coords_gl, do_column_diagnostics,  &
                      col_diag_lon, col_diag_lat, col_diag_i,  &
                      col_diag_j, col_diag_unit)

!---------------------------------------------------------------------
!    verify that requested pressure cutoff for column diagnostics output
!    is valid. define the model k index which corresponds.
!---------------------------------------------------------------------
        do k=1,nlev
          if (pref(k) >= diagnostics_pressure_cutoff) then
            kstart_diag = k
            exit
          endif
        end do
        if (kstart_diag == -99) then
          call error_mesg ( 'donner_deep_mod', &
           'diagnostics_pressure_cutoff is higher than pressure at '//&
                                     'any model level', FATAL)
        endif



end subroutine donner_column_init


!#####################################################################

subroutine donner_column_control (is, js, Time)

integer, intent(in) :: is, js
type(time_type), intent(in) :: Time

      integer :: nn, j, i

!--------------------------------------------------------------------
!    initialize debug control in this subdomain. if debug is desired, 
!    this is the requested debug timestep, and data from this window
!    is desired as output, in_diagnostics_window will be true and j_dc
!    is set to the debug column's j index. 
!-------------------------------------------------------------------
      in_diagnostics_window = .false.
      j_dc(:) = -99
      i_dc(:) = -99
      ncols_in_window = 0
      if (column_diagnostics_desired) then
        if (Time >= Time_col_diagnostics) then
          do nn =1, num_diag_pts
            do j=1,jsize      
              if (js + j - 1 == col_diag_j(nn)) then
                do i=1,isize       
                  if (is+i-1 == col_diag_i(nn)) then
                    ncols_in_window = ncols_in_window + 1
                    i_dc(ncols_in_window) = i
                    j_dc(ncols_in_window) = j
                    in_diagnostics_window = .true.
                    igl_dc(ncols_in_window) = col_diag_i(nn)
                    jgl_dc(ncols_in_window) = col_diag_j(nn)
                    unit_dc(ncols_in_window) = col_diag_unit(nn)
                    call column_diagnostics_header &
                     (mod_name, col_diag_unit(nn), Time, nn,  &
                      col_diag_lon, col_diag_lat, col_diag_i,  &
                      col_diag_j)
                  endif
                end do  ! (i loop)
              endif
            end do  ! (j loop) 
          end do  ! (num_diag_pts loop)
        endif  ! (Time >= starting time)
      endif ! (diagnostics desired)

end subroutine donner_column_control


!#####################################################################

subroutine donner_column_cape_call (tempbl, ratpbl, ttnd, qtnd)

real, dimension(:,:,:), intent(in) :: tempbl, ratpbl, ttnd, qtnd

      integer :: n, k
      if (in_diagnostics_window) then
        do n=1,ncols_in_window
          do k=1,model_levels_in_sfcbl
            write (unit_dc(n), '(a, i4, f20.14, e20.12)')  &
                        'in donner_deep: k, tempbl,ratpbl: ',  &
                         nlev-k+1, tempbl(igl_dc(n),jgl_dc(n),k), &
                         ratpbl(igl_dc(n),jgl_dc(n),k)
          end do
          do k=1,nlev
            if (ttnd(i_dc(n), j_dc(n),k) /= 0.0) then
              write (unit_dc(n), '(a, i4, f20.14, e20.12)') &
                                'in donner_deep: k,ttnd,qtnd', &
                       k,ttnd(i_dc(n),j_dc(n),k),qtnd(i_dc(n),j_dc(n),k)
            endif
          end do
        end do
      endif

end subroutine donner_column_cape_call 

!#####################################################################

subroutine donner_column_input_fields (dt, conv_calc_on_this_step, &
                                       temp, mixing_ratio, phalf, omega)

real, intent(in) :: dt
logical, intent(in) :: conv_calc_on_this_step
real, dimension(:,:,:), intent(in) :: temp, mixing_ratio, phalf
real, dimension(:,:,:), intent(in) :: omega                  

      integer :: n, k
      do n=1,ncols_in_window
        write (unit_dc(n), '(a,f8.1, 2i4)')  &
                ' physics timestep, window i, window j= ',  &
                  dt, i_dc(n), j_dc(n)
        write (unit_dc(n),'(a,l4 )' )  'conv_calc_on_this_step = ',   &
                  conv_calc_on_this_step
        do k=kstart_diag,nlev
          write (unit_dc(n), '(a, i4, f20.14, e20.12)')  &
                       'in donner_deep A: k,temp  ,mixing ratio', &
             k,temp  (i_dc(n),j_dc(n),k),mixing_ratio(i_dc(n),j_dc(n),k)
        end do
        write (unit_dc(n),'(a,f19.10,2e20.12)')  &
                             'sfcprs,  omega_btm= ',   &
                            phalf(i_dc(n),j_dc(n),nlev+1),   &
!                           omega_btm(i_dc(n),j_dc(n)) 
                            omega(i_dc(n),j_dc(n), nlev) 
        write (unit_dc(n),'(a,f19.10,2e20.12)')  ' omint= ',   &
                            omint_acc(igl_dc(n),jgl_dc(n))
      end do

end subroutine donner_column_input_fields 


!####################################################################

subroutine donner_column_end_of_step (Don_conv, Don_Cape)

type(donner_cape_type), intent(in) :: Don_cape
type(donner_conv_type), intent(in) :: Don_conv


       integer :: k, n, kcont



       do n=1,ncols_in_window
         write (unit_dc(n), '(a, e20.12)')  & 
                'in donner_deep: plcl ', Don_cape%plcl(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
                 'in donner_deep: plfc ', Don_cape%plfc(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
              'in donner_deep: plzb ', Don_cape%plzb(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
               'in donner_deep: xcape ', xcape_lag(igl_dc(n),jgl_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
               'in donner_deep: coin ', Don_cape%coin(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
              'in donner_deep: dcape ', Don_conv%dcape(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
              'in donner_deep: qint ', qint_lag(igl_dc(n),jgl_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
               'in donner_deep: a1   ', Don_conv%a1  (i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
            'in donner_deep: amax ', Don_conv%amax(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
            'in donner_deep: amos ', Don_conv%amos(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
            'in donner_deep: tprea1 ', tprea1(igl_dc(n),jgl_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
             'in donner_deep: ampta1 ', Don_conv%ampta1(i_dc(n),j_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
                'in donner_deep: omint', omint_acc(igl_dc(n),jgl_dc(n))
         write (unit_dc(n), '(a, e20.12)')  & 
             'in donner_deep: rcoa1 ', Don_conv%rcoa1(i_dc(n),j_dc(n))

         do k=kstart_diag,nlev
           write (unit_dc(n), '(a, i4)')'in donner_deep: k = ', k
           write (unit_dc(n), '(a, e20.12)')  &
                    'in donner_deep: cemetf',  &
                              cemetf (igl_dc(n),jgl_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                      'in donner_deep: ceefc ',     &
                               Don_conv%ceefc  (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                      'in donner_deep: cecon ',  &
                                Don_conv%cecon  (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                     'in donner_deep: cemfc ',   &
                               Don_conv%cemfc  (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                     'in donner_deep: cememf',  &
                                 cememf (igl_dc(n),jgl_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                      'in donner_deep: cememf_mod',  &
                                Don_conv%cememf_mod (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                       'in donner_deep: cual  ',  &
                                  Don_conv%cual(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                       'in donner_deep: fre   ',   &
                                  Don_conv%fre    (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                        'in donner_deep: elt   ',  &
                                     Don_conv%elt    (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                         'in donner_deep: cmus  ',    &
                                     Don_conv%cmus   (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                         'in donner_deep ',   &
                                    Don_conv%ecds   (i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                         'in donner_deep: eces  ', &
                                     Don_conv%eces(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                         'in donner_deep: emds  ',  &
                                     Don_conv%emds(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                          'in donner_deep: emes  ',  &
                                     Don_conv%emes(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                          'in donner_deep: qmes  ',  &
                                    Don_conv%qmes(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                          'in donner_deep: wmps  ', &
                                      Don_conv%wmps(i_dc(n),j_dc(n),k)
          write (unit_dc(n), '(a, e20.12)')  &
                          'in donner_deep: wmms  ',  &
                                       Don_conv%wmms(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                           'in donner_deep: tmes  ',  &
                                        Don_conv%tmes(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                             'in donner_deep: dmeml ',   &
                                       Don_conv%dmeml(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                              'in donner_deep: uceml ',  &
                                       Don_conv%uceml(i_dc(n),j_dc(n),k)
           write (unit_dc(n), '(a, e20.12)')  &
                               'in donner_deep: umeml ',   &
                                      Don_conv%umeml(i_dc(n),j_dc(n),k)
           do kcont=1,ncont
           write (unit_dc(n), '(a, e20.12)')  &
                               'in donner_deep: xgcm1 ',   &
                           Don_conv%xgcm1(i_dc(n),j_dc(n),k,kcont)
           write (unit_dc(n), '(a, e20.12)')  &
                               'in donner_deep: qtren1 ',  &
                              Don_conv%qtren1(i_dc(n),j_dc(n), k,kcont)
           write (unit_dc(n), '(a, e20.12)')  &
                                'in donner_deep: qtmes1 ',  &
                              Don_conv%qtmes1(i_dc(n),j_dc(n), k,kcont)
!          end do
!          do kcont=1,ncont
             write (unit_dc(n), '(a, e20.12)')  &
                                   'in donner_deep: qtceme ',   &
                               Don_conv%qtceme(i_dc(n),j_dc(n),k,kcont)
           write (unit_dc(n), '(a, e20.12)')  &
                                  'in donner_deep: wtp1 ',   &
                                Don_conv%wtp1(i_dc(n),j_dc(n),k,kcont)
           end do


         end do
       end do



end subroutine donner_column_end_of_step




!######################################################################


                 end module donner_deep_mod

