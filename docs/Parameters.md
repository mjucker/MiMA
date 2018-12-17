[back to top](https://mjucker.github.io/MiMA)

# Parameter settings

This section shows most of the parameters and their default and/or recommended values. For details about the physical meaning of these parameters, please refer to the main MiMA reference paper and comments in the source code.

Go directly to [default values](#default-values) or an [input file with all recommended values](../input/input.nml).

## Recommended values

Some of the default values are "safe choices", designed to have the least impact if the user is not aware of them.
But these values are not necessarily the best ones, so this section gives some ballpark values which have been working well so far.
A sample input file with all recommended values resides in the repository [here](../input/input.nml).

### General

Namelist `coupler_nml`

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 dt_atmos | 900 | integration time step in [s]
 
### Dynamics
 
 Namelist `spectral_dynamics_nml`
 
 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 damping_order | 4 | 8<sup>th</sup> order numerical diffusion
 water_correction_limit | 200.e2 | [Pa] turn off water correction in stratosphere to avoid systematic sink
 initial_sphum | 2.e-6 | [kg/kg] Start with some water vapor as the stratosphere takes very long to fill up dynamically
 robert_coeff            | .03 | Robertson coefficient for implicit time stepping
 num_levels | >= 40 | Number of vertical levels
 vert_coord_option | 'uneven_sigma' | use hybrid sigma/pressure coordinates
 surf_res   | 0.5 | parameter 1 to define vertical level distribution
 scale_heights | 9.0 | parameter 2 to define vertical level distribution
 exponent | 7.0 | parameter 3 to define vertical level distribution
 topography_option | 'gaussian' | if orographic forcing required
 specify_initial_conditions | .false. | If set to true, the model will look for initial conditions in INPUT/initial_conditions.nc


 Namelist `gaussian_topog_nml`, example for 3km wave-two in midlatitudes
 
 Variable | Recommended Value | Meaning
 :--- | :---: | :---
   height   | 3000., 3000., 
   olat     |   45.,   45.,
   olon     |   90.,  270.,
   wlat     |   20.,   20.,
   wlon     |   20.,   20.,
   rlat     |    0.,    0.,
   rlon     |    0.,    0., /

  

 
### Mixed layer

Namelist `simple_surface_nml`

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 Tm   | 285 | [K] some meridional profile already in initial conditions
 do_qflux | .true. | export tropical heat into extratropics. Mimic ocean circulation
 do_warmpool | .true. | mimic tropical warmpool wave-one forcing
 const_albedo | 0.22 | global SSTs reasonable
 [land_option | 'lonlat' | if some land-sea contrast desired ]
 [slandlon    |   0, | [deg] don't forget the ',' ]
 [elandlon    | 360, | [deg] zonally symmetric land]
 [slandlat    |  40, | [deg] way outside the tropics ]
 [elandlat    |  50, | [deg] midlats only]
 
Namelist `qflux_nml`

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 qflux_amp | 30 | [W/m<sup>2</sup>] gives a good compromise between too strong jets and double ITCZ 
 warmpool_amp | 30 | [W/m<sup>2</sup> gives realistic warmpool anomaly and good cold point

### Moisture

Following *Frierson (2007)* we use large scale condenstion together with the Betts-Miller convection scheme with the "shallower" shallow convection adjustment.

Namelist `moist_processes_nml`

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 do_mca | .false. | Do moist convective adjustment
 do_lsc | .true. | Do large scale condensation
 do_bm  | .true. | Do Betts-Miller convetion
 use_df_stuff | .true. | When true, specific humidity = (rdgas/rvgas)*esat/pressure

Namelist `betts_miller_nml`

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 rhbm   | 0.7  | relax to 70% relative humidity
 do_simp | .false. | don't adjust time scales to make precipitation always continuous
 do_shallower | .true. | shallow convection: choose smaller depth to make precip zero

Namelist `moist_conv_nml` is only touched to make sure moisture handling is consistent accross namelists.

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 use_df_stuff | .true. | Make everything consistent with above `use_df_stuff`

Namelist `lscale_cond_nml`: We want to re-evaporate outfalling precipitation if any of the layers below are sub-saturated.

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 do_evap | .true. | re-evaporate in below sub-saturated layers (if any)
 use_df_stuff | .true. | Make everything consistent
 
 
### Radiation

Namelist `rrtm_radiation_nml`

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
  do_read_ozone         |    .true.      |  Read ozone from file (only way to have non-zero ozone)
  ozone_file            |    ’ozone_1990’     | If so, filename without ’.nc’ extension - this file is in the repository
  co2_val               |    300.           |  Constant value for CO<sub>2</sub> [ppm]
  dt_rad                |    7200           |  Radiation time step [s]. Every time step if < `dt_atmos`

Namelist `astro_nml`.

Variable | Recommended Value | Meaning
 :--- | :---: | :---
 solr_cnst  | 1360           | solar constant [W/m2]
 [ solday     | 90           | if perpetual equinox desired and calendar='thirt_day' ]


### Boundary conditions

Namelist `damping_driver_nml`.

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 do_rayleigh | .true. | do simple Rayleigh friction at the top
 trayfric | -0.5 | Rayleigh friction time scale of 1/2 day
 sponge_bottom | 100 | Rayleigh friction above 0.5hPa
 do_conserve_energy | .true. | account for heat release due to momentum loss
 
Namelist `surface_flux_nml`

Variable | Recommended Value | Meaning
 :--- | :---: | :---
 use_virtual_temp | .false. | for consistency with df_stuff
 old_dtaudv       | .true. | use alternative d(stress)/d(wind component)
 use_df_stuff     | .true. | for consistency with df_stuff
 no_surface_momentum_flux  | .false. | use to turn off surface momentum fluxes
 no_surface_moisture_flux  | .false. | use to turn off surface moisture fluxes
 no_surface_heat_flux      | .false. | use to turn off surface heat fluxes
 no_surface_radiative_flux | .false. | use to turn off surface radiative fluxes


Namelist `diffusivity_nml`

Variable | Recommended Value | Meaning
 :--- | :---: | :---
do_entrain          | .false. | don't account for entrainment (which is off anyway)
use_df_stuff        | .true. | for consistency with df_stuff


### Vertical numerics

Namelist `vert_turb_driver_nml`

Variable | Recommended Value
 :--- | :---: 
use_tau          | .false.
constant_gust | 0.0
do_mellor_yamada | .false.
use_df_stuff  | .true.
do_diffusivity         | .true.

Namelist `vert_diff_driver_nml`

Variable | Recommended Value
:-- | :--:
do_conserve_energy         | .true.
use_virtual_temp_vert_diff | .false.



## Default values

### General

Most general parameters are set in `coupler/coupler_main.f90` and `coupler_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 current_date | 0001,1,1,0,0,0 | Start from year 1
 calendar | 'thirty_day' | Use 12 30-day months

### Dynamics

The dynamical core is set up in `atmos_spectral/model/spectral_dynamics.f90` and the namelist `spectral_dynamics_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
   damping_order           | 2 | 8<sup>th</sup> order numerical diffusion
   do_water_correction     | .true. | make sure water mass is conserved in the advection step
   water_correction_limit  | 200.e2 | correct water mass only below this level [Pa]. Introduces artificial sink in stratosphere if corrected there. 
   vert_advect_uv          | 'second_centered' | second order vertical advection scheme
   vert_advect_t           | 'second_centered' | second order vertical advection scheme
   robert_coeff            | .04 | Robertson coefficient for implicit time stepping
   lon_max                 | 128 | T42 resolution
   lat_max                 | 64  | T42 resolution
   num_fourier             | 42  | T42 resolution
   num_spherical           | 43  | T42 resolution
   fourier_inc             | 1   | T42 resolution
   num_levels              | 18  | number of vertical levels
   surf_res        | 0.1 | parameter 1 to define uneven_sigma levels
   scale_heights | 4 | parameter 2 to define uneven_sigma levels
   vert_coord_option       | 'even_sigma' | use sigma levels
   exponent | 7.0 | parameter 3 to define vertical level distribution

Initial conditions are set in `atmos_spectral/init/spectral_init_cond.f90` and the namelist `spectral_init_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 initial_temperature | 264 | Initial temperature of isothermal atmosphere in [K]

### Physics 

The namelist `physics_driver_nml` steers which physics components are used. It resides in `atmos_param/physics_driver/physics_driver.f90`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 do_moist_processes | .true. | call moist_processes routines
 tau_diff | 3600.    | time scale for smoothing diffusion coefficients
 do_radiation | .false. | calculating radiative fluxes and heating rates with AM2 radiation?
 do_grey_radiation | .false. | rather do grey radiation?
 do_rrtm_radiation | .true. | or RRTM radiation?
 do_damping | .true. | do any of the damping schemes?
 do_local_heating | .false. | add artificial local heating? If so, see `local_heating_nml` namelist
 diff_min | 1.e-3    | minimum value of a diffusion coefficient beneath which the coefficient is reset to zero
 diffusion_smooth | .true. | diffusion coefficients should be smoothed in time?
 do_netcdf_restart | .true. | make restart files netCDF format?

### Mixed layer

All changes to the mixed layer ocean as compared to FHZ06 are discussed
in the main text. Table \[Surface\_default\] gives the default values
for all newly introduced input parameters, and those that might have
slightly changed meaning.

These parameters are set in `coupler/simple_surface.f90`.

  Variable           | Default Value | Meaning
  :--- | :---: | :---
  Tm              | 265 | Initial surface temperature [K]. 1K warmer than isothermal atmosphere seems reasonable to get convection going right away.
  heat_capacity   | 4e8 | [J/K/m<sup>2</sup>] 100m mixed layer depth
  land_capacity   | -1 | same as `heat_capacity`
  trop_capacity   | -1 | same as `heat_capacity`
  trop_cap_limit  | 15 | [deg] Poleward boundary for `trop_capacity`
  heat_cap_limit  | 60 | [deg] Equatorward boundary for `heat_capacity`
  albedo_choice   | 1  | Single value globally
  const_albedo    | 0.30 | Value of surface albedo
  albedo_exp      | 2  | Exponent for meridional albedo profile if `albedo_choice = 4`
  do_qflux        | .false. | Add meridional heat flux to surface?
  qflux_width     | 16 | if so, half-width in [degrees lat]
  qflux_amp       | 30 | if so, amplitude in [W/m<sup>2</sup>]
  do_warmpool     | .false. | Add zonally asymmetric heat flux?
  warmpool_amp     | 30 | if so, amplitude in [W/m<sup>2</sup>]
  warmpool_phase   | 0  | if so, phase in zonal direction in degrees longitude
  warmpool_centr  | 0  | if so, center of heatflux in degrees latitude
  warmpool_k      | 1  | if so, wave number of heatflux
  land_option     | ’none’ | Add land-sea contrast?
  slandlon        | (0,0,0,0,0,0,0,0,0,0) | if so, give land patches start longitude [deg]
  elandlon        | -(1,1,1,1,1,1,1,1,1,1)| if so, give land patches end longitude [deg]
  slandlat        | (0,0,0,0,0,0,0,0,0,0) | if so, give land patches start latitude [deg]
  elandlat        | -(1,1,1,1,1,1,1,1,1,1)| if so, gie land patches end longitude [deg]
 

### Moisture

In contrast to the initial *Frierson et al (2006)* setup, we turn off moist convective adjustment by default, and turn on Betts-Miller convection and use the alternative definition of specific humidity.

#### Moist processes

The moist processes parameters are set in `atmos_param/moist_processes/moist_processes.f90` and the namelist `moist_processes_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 do_mca | .true. | Do moist convective adjustment?
 do_lsc | .true. | Do large scale condensation?
 do_ras | .false. | Do relaxed Arakawa-Schubert?
 do_strat | .false. | Do stratiform clouds?
 do_dryadj | .false. | Do dry adjustment?
 do_rh_clodus | .false. | Do relative humidity cloud scheme?
 do_diag_clouds | .false. | Do Gordon's diagnostic cloud scheme?
 do_donner_deep | .false. | Do Donner deep convection scheme?
 use_tau | .false. | Use current time values? (future time if .false.)
 do_gust_cv | .false. | Do convective gustiness?
 do_bm  | .false. | Do Betts-Miller convetion?
 do_bmmass | .false. | Do Betts-Miller mass flux scheme?
 do_bmomp | .false. | Do Pauluis version of Betts-Miller scheme?
 use_df_stuff | .false. | When true, specific humidity = (rdgas/rvgas)*esat/pressure
 
 
#### Betts-Miller 

The Betts-Miller parameters are set in `atmos_param/betts_miller/betts_miller.f90` and the namelist `betts_miller_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 tau_bm | 7200 | relaxation time scale in [s]
 rhbm   | 0.8  | relative humidity to relax to
 do_simp | .true. | adjust time scales to make precipitation always continuous?
 do_shallower | .false. | shallow convection: choose smaller depth to make precip zero?
 do_changeqref | .false. | shallow convection: change both q and T to make precip zero?
 do_envsat | .false. | reference rhbm wrt environment (.true.) or parcel (.false.)?
 do_taucape | .false. | make taubm proportional to 1/sqrt(CAPE)?
 capetaubm | 900 | [J/kg] value of CAPE for which tau = tau_bm if do_taucape=.true.
 tau_min | 2400 | [s] minimum tau_bm allowed if do_taucape=.true.
 
#### Moist convective adjustment

Moist convective adjustment parameters are in `atmos_param/moist_conv/moist_conv.f90` and the namelist `moist_conv_nml`. 

 Variable | Default Value | Meaning
 :--- | :---: | :---
 beta | 0.0 | fraction of condensate detrained into stratiform cloud
 use_df_stuff | .false. | 

#### Large Scale Condensation

We want to re-evaporate outfalling precipitation if any of the layers below are sub-saturated. Parameters are described in `atmos_param/lscale_cond/lscale_cond.f90` and the namelist `lscale_cond_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 hc | 1.0 | relative humidity at which condensation occurs
 do_evap | .false. | re-evaporate in below sub-saturated layers (if any)?
 use_df_stuff | .false. | For consistency with above

### Radiation

#### RRTM

The new radiative transfer routine, which eventually calls the RRTM
short and long wave modules, has a few parameters that do not affect
RRTM directly. Those parameters are given in the first table and can be found in `atmos_param/rrtm_radiation/rrtm_radiation.f90` and namelist `rrtm_radiation_nml`.


 Variable | Default Value | Meaning
 :--- | :---: | :---
  do_read_radiation     |    .false.     |  Read SW and LW radiation in the atmosphere from file?
  radiation_file        |    ’radiation’ |  If so, filename without ’.nc’ extension
  do_read_sw_flux       |    .false.     |  Read SW surface flux from file?
  sw_flux_file          |    ’sw_flux’  |  If so, filename without ’.nc’ extension
  do_read_lw_flux       |    .false.     |  Read LW surface flux from file?
  lw_flux_file          |    ’lw_flux’  |  If so, filename without ’.nc’ extension
  do_read_ozone         |    .false.     |  Read ozone from file (only way to have non-zero ozone)?
  ozone_file            |    ’ozone’     | If so, filename without ’.nc’ extension
  do_read_h2o           |    .false.     |  Use external water vapor distribution instead of active tracer?
  h2o_file              |    ’h2o’       |  If so, filename without ’.nc’ extension
  include_secondary_gases|   .false.     |  Set CH<sub>4</sub>, N<sub>2</sub>O, O<sub>2</sub>, CFC-11, CFC-12, CFC-22, CCl<sub>4</sub> to non-zero?
  ch4_val               |    0           |  If so, set value for CH<sub>4</sub>
  n2o_val               |    0           |  If so, set value for N<sub>2</sub>O
  o2_val                |    0           |  If so, set value for O<sub>2</sub>
  cfc{1112,22}_val      |    0           |  If so, set value for CFC-{11,12,22}
  ccl4_val              |    0            | If so, set value for CCl<sub>4</sub>
  h2o_lower_limit       |    0.2ppm      |  Never use specific humidity values smaller than this in radiative transfer
  do_fixed_water        |    .false.     |  Use fixed value for specific humidity in radiative transfer?
  fixed_water           |    2ppm        |  If so, set value
  fixed_water_pres      |    100e2Pa     |  If so, above which pressure level?
  fixed_water_lat       |    90          |  If so, equatorward of which latitude?
  do_zm_tracers         |    .false.     |  Feed only the zonal mean of all absorbers to the radiative transfer
  do_zm_rad             |    .false.     |  Only pass zonal mean radiative forcing to dynamics
  dt_rad                |    0           |  Radiation time step [s]. Every time step if < `dt_atmos`
  store_intermediate_rad |   .true.      |  Keep radiative forcing constant or set to zero between radiative time steps?
  do_rad_time_avg       |    .true.      |  Compute zenith angle average over radiative time step (`true`) or use instantaneous (`false`)?
  dt_rad_avg            |    86400.      |  Average zenith angle over this time interval [s]. Defaults to diurnal mean, e.g. incoming SW flux at the TOA will be zonally symmetric. Set equal to `dt_rad` if < 0
  lonstep               |    1           |  Only compute radiation at every nth longitudinal grid point
  do_precip_albedo      |    .false.     |  Link surface albedo to precipitation?
  precip_albedo         |    0.35        |  If so, set albedo for 100% precipitating grid boxes
  precip_lat            |    0           |  If so, only use poleward of this latitude [deg]
  precip_albedo_mode     |   ’full’      |  If so, use ’full’ (total), ’lscale’ (large scale), or ’conv’ (convective) precipitation for albedo calculation


The next table gives the default values for all needed RRTM
input parameters, with the parameter names as given by the original
source code. These are not part of any namelist, but we report them here for completeness.
We refer the reader to RRTM documentation and source code
to learn the exact meanings of each input parameter. The parameter type
is defined as follows: ‘fixed’ means the value is hard coded and cannot
be changed without recompiling; ‘fixed input’ means that the variable
can be changed within a namelist at run time, but will not change during
the simulation; ‘interactive’ means that the value is re-evaluated at
each radiation time step before calling RRTM.

  Variable   |  Type     |     Default
:--- | :--- | :---
  iaer      | fixed        | 0
  h2vmr     | interactive  | tracer
  o3vmr     | fixed input  | from file
  co2vmr    | fixed input  | 300e-6
  ch4vmr    | fixed input  | 0
  n2ovmr    | fixed input  | 0
  o2vmr     | fixed input  | 0
  cfc11vmr  | fixed input  | 0
  cfc12vmr  | fixed input  | 0
  cfc22vmr  | fixed input  | 0
  ccl4vmr   | fixed input  | 0
  asdir     | fixed input  | 0.27
  aldir     | fixed        | = asdir
  asdif     | fixed        | = asdir
  aldir     | fixed        | = asdir
  coszen    | interactive  | computed
  adjes     | fixed input  | 1
  dyofyr    | fixed input  | 0
  scon      | fixed input  | 1368.22
  emis      | fixed        | 1
  inflgsw   | fixed        | 0
  iceflgsw  | fixed        | 0
  liqflgsw  | fixed        | 0
  cldfr     | fixed        | 0
  taucld    | fixed        | 0
  ssacld    | fixed        | 0
  asmcld    | fixed        | 0
  asmcld    | fixed        | 0
  fsfcld    | fixed        | 0
  cicewp    | fixed        | 0
  cliqwp    | fixed        | 0
  reice     | fixed        | 10
  reliq     | fixed        | 10
  tauaer    | fixed        | 0
  ssaaer    | fixed        | 0
  asmaer    | fixed        | 0
  acaer     | fixed        | 0

#### Astronomy

All things which define Earth versus other planets/planetary systems are set in `atmos_param/rrtm_radiation/astro.f90` and `astro_nml`.

Variable | Default Value | Meaning
 :--- | :---: | :---
obliq      | 23.439             | Earth's obliquity in [degrees latitude]
use_dyofyr | .false.            | use day of the year to compute Earth-Sun distance? Note that this is done internally in RRTM, and assumes 365days/year.
solrad     | 1.0                | distance Earth-Sun [AU] if use_dyofyr=.false.
solr_cnst  | 1368.22            | solar constant [W/m2]
solday     | 0                  | if >0, do perpetual run corresponding to day of the year = solday in [0,days per year]
equinox_day | 0.25              | fraction of the year defining March equinox.
        

#### Local heating
      

If `do_local_heating = .true.` in `physics_driver_nml`, the namelist `local_heating_nml` can be used to set the form and position of the desired local heating.

Variable | Default Value | Meaning
:--- | :---: | :---
hamp | 0 | amplitude of Gaussian heating in [K/d], maximum 10 entries
loncenter| -1 | zonal center of the Gaussian in [degrees longitude]. Zonally symmetric if <0. maximum 10 entries
lonwidth | -1 | zonal width of the Gaussian, if loncenter >= 0. [degrees longitude], maximum 10 entries
lonmove  | 0  | zonal speed of moving heat source [deg/day], maximum 10 entries
latcenter|  0 | meridional center of the Gaussian in [degrees latitude], maximum 10 entries
latwidth | 15 | meridional width of the Gaussian in the [degrees latitude], maximum 10 entries
latmove  | 0  | meridonal speed of moving heat source [deg/day], maximum 10 entries
pcenter | -1 | vertical center of the Gaussian in the vertical [hPa] no vertical structure if <0, maximum 10 entries
pwidth  | 1  | vertical width of the Gaussian in orders of magnitude [log10(hPa)], constant with pressure if <0, maximum 10 entries
pmove   | 0  | vertical speed of moving heat source [hPa/day], maximum 10 entries
is_periodic | .false. | reset location periodically? Periodicity is unidirectional in longitude and pressure, back-and-forth in latitude, maximum 10 entries
twidth  | -1 | temporal width of Gaussian [days]. constant in time if <0, maximum 10 entries
tphase  | 0  | temporal phase of Gaussian heating [days], maximum 10 entries
tperiod | -1 | temporal period of Gaussian heating; [fraction of year] if <0, [days] if >0, maximum 10 entries

      


### Boundary conditions

The upper boundary conditions are defined by how to deal with gravity waves. The simplest choice is simple Rayleigh friction (damping towards zero momentum), non-orographic gravity wave parameterization, orographic gravity wave drag, or a constant "gravity wave" drag.
By default, no scheme is active, which will almost certainly cause the code to crash. See [recommended values](#recommended-values) above and source code comments for guidance.

Parameters are described in `atmos_param/damping_driver/damping_driver.f90` and the namelist `damping_driver_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 do_rayleigh | .false. | do simple Rayleigh friction at the top?
 trayfric | 0.0 | Rayleigh friction time scale, [s] if > 0, [-d] if < 0
 sponge_bottom | 50 | [Pa] bottom of Rayleigh friction layer (sponge)
 do_mg_drag | .false. | mountain gravity wave scheme (not tested!)
 do_cg_drag | .false. | non-orographic gravity wave scheme (under development)
 do_topo_drag | .false. | topographic drag scheme (not tested!)
 do_const_drag | .false. | constant "gravity wave" scheme (not tested!)
 do_conserve_energy | .false. | account for heat release due to momentum loss?

The lower boundary condition is set by the Monin-Obukhov boundary layer. It is used for atmospheric diffusivities in `atmos_param/diffusivity/diffusivity.f90` and surface fluxes in `coupler/surface_flux.f90`.  We didn't change anything in this part as compared to *Frierson et al (2006)*, and only report the default values here. Refer to [recommended values](#recommended-values) for the appropriate settings in MiMA.

Namelist `surface_flux_nml`

Variable | Default Value 
 :--- | :---: 
 no_neg_q | .false.  
 use_virtual_temp | .true. 
 alt_gustiness    | .false. 
 old_dtaudv       | .false.
 use_mixing_ratio | .false. 
 use_df_stuff     | .false. 
 gust_const       |  1.0 
 ncar_ocean_flux  | .false. 
 raoult_sat_vap   | .false. 

Namelist `diffusivity_nml`

Variable | Default Value
 :--- | :---: 
fixed_depth | .false. 
depth_0     | 5000.0 
frac_inner  | 0.1 
rich_crit_pbl | 1.0 
entr_ratio    |  0.2 
parcel_buoy         |  2.0 
znom                |  1000.0 
free_atm_diff       | .false. 
free_atm_skyhi_diff | .false. 
pbl_mcm             | .false. 
rich_crit_diff      |  0.25 
mix_len             | 30. 
rich_prandtl        |  1.00 
background_m        |  0.0 
background_t        |  0.0 
ampns               | .false. 
ampns_max           | 1.0E20  
do_entrain          | .true.
use_df_stuff        | .false.


### Vertical numerics

Again, nothing has been done to these schemes, but we report the default values. 

File `atmos_param/vert_turb_driver/vert_turb_driver.f90` and namelist `vert_turb_driver_nml`. Make sure to check the [recommended values](#recommended-values).

Variable | Default Value
 :--- | :---: 
do_shallow_conv  | .false.
do_mellor_yamada | .true.
do_diffusivity         | .false.
do_molecular_diffusion | .false.
do_edt                 | .false.
do_stable_bl     | .false.
use_tau          | .true.
do_entrain    | .false.
gust_scheme  | 'constant' 
constant_gust | 1.0
gust_factor   | 1.0
use_df_stuff  | .false.


File `atmos_param/vert_diff_driver/vert_diff_driver.f90` and namelist `vert_diff_driver_nml`. Make sure to check the [recommended values](#recommended-values).

Variable | Default Value
:-- | :--:
do_conserve_energy         | .false.
do_mcm_no_neg_q            | .false.
use_virtual_temp_vert_diff | .true.
do_mcm_plev                | .false.
do_mcm_vert_diff_tq        | .false.
