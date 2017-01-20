# Parameter settings

This section shows most of the parameters and their default and/or recommended values. For details about the physical meaning of these parameters, please refer to the main MiMA reference paper and comments in the source code.

## Recommended values

Some of the default values are "safe choices", designed to have the least impact if the user is not aware of them.
But these values are not necessarily the best ones, so this section gives some ballpark values which have been working well so far.

### General

Namelist `coupler_nml`

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 dt_atmos | 600 | integration time step in [s]
 
### Dynamics
 
 Namelist `spectral_dynamics_nml`
 
 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 water_correction_limit | 200.e2 | [Pa] turn off water correction in stratosphere to avoid systematic sink
 initial_sphum | 2.e6 | [kg/kg] Start with some water vapor as the stratosphere takes very long to fill up dynamically
 num_levels | >= 40 | Number of vertical levels
 vert_coord_option | 'uneven_sigma' | use hybrid sigma/pressure coordinates
 surf_res   | 0.4 | parameter 1 to define vertical level distribution
 scale_heights | 11.0 | parameter 2 to define vertical level distribution
 topography_option | 'gaussian' | if orographic forcing required
 
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

### Boundary conditions

Namelist `damping_driver_nml`.

 Variable | Recommended Value | Meaning
 :--- | :---: | :---
 do_rayleigh | .true. | do simple Rayleigh friction at the top
 trayfric | -0.5 | Rayleigh friction time scale of 1/2 day
 sponge_bottom | 50 | Rayleigh friction above 0.5hPa
 do_conserve_energy | .true. | account for heat release due to momentum loss
 
Namelist `surface_flux_nml`

Variable | Recommended Value | Meaning
 :--- | :---: | :---
 use_virtual_temp | .false. | for consistency with df_stuff
 old_dtaudv       | .true. | use alternative d(stress)/d(wind component)
 use_df_stuff     | .true. | for consistency with df_stuff

Namelist `diffusivity_nml`

Variable | Recommended Value | Meaning
 :--- | :---: | :---
do_entrain          | .false. | don't account for entrainment (which is off anyway)
use_df_stuff        | .true. | for consistency with df_stuff


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
   damping_order           | 4 | 8<sup>th</sup> order numerical diffusion
   do_water_correction     | .true. | make sure water mass is conserved in the advection step
   water_correction_limit  | 200.e2 | correct water mass only below this level [Pa]. Introduces artificial sink in stratosphere if corrected there. 
   vert_advect_uv          | 'second_centered' | second order vertical advection scheme
   vert_advect_t           | 'second_centered' | second order vertical advection scheme
   robert_coeff            | .03 | Robertson coefficient for implicit time stepping
   lon_max                 | 128 | T42 resolution
   lat_max                 | 64  | T42 resolution
   num_fourier             | 42  | T42 resolution
   num_spherical           | 43  | T42 resolution
   fourier_inc             | 1   | T42 resolution

Initial conditions are set in `atmos_spectral/init/spectral_init_cond.f90` and the namelist `spectral_init_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 initial_temperature | 264 | Initial temperature of isothermal atmosphere in [K]


### Mixed layer

All changes to the mixed layer ocean as compared to FHZ06 are discussed
in the main text. Table \[Surface\_default\] gives the default values
for all newly introduced input parameters, and those that might have
slightly changed meaning.

These parameters are set in `coupler/simple_surface.f90`.

  Variable           | Default Value | Meaning
  :--- | :---: | :---
  Tm              | 265 | Initial surface temperature [K]. 1K warmer than isothermal atmosphere seems reasonable to get convection going right away.
  heat_capacity   | 4e-8 | [J/K/m<sup>2</sup>] 100m mixed layer depth
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


### Radiation

The new radiative transfer routine, which eventually calls the RRTM
short and long wave modules, has a few parameters that do not affect
RRTM directly, and are given with their default values in Table
\[Rad\_default\]. All of these variables can be set within the
corresponding namelist in the input file.

  Variable                    Default       Meaning
  --------------------------- ------------- ----------------------------------------------------------------------------------------------------------------
  `do_read_radiation`         .false.       Read SW and LW radiation in the atmosphere from file?
  `radiation_file`            ’radiation’   If so, filename without ’.nc’ extension
  `do_read_sw_flux`           .false.       Read SW surface flux from file?
  `sw_flux_file`              ’sw\_flux’    If so, filename without ’.nc’ extension
  `do_read_lw_flux`           .false.       Read LW surface flux from file?
  `lw_flux_file`              ’lw\_flux’    If so, filename without ’.nc’ extension
  `do_read_ozone`             .false.       Read ozone from file (only way to have non-zero ozone)?
  `ozone_file`                ’ozone’       If so, filename without ’.nc’ extension
  `do_read_h2o`               .false.       Use external water vapor distribution instead of active tracer?
  `h2o_file`                  ’h2o’         If so, filename without ’.nc’ extension
  `include_secondary_gases`   .false.       Set CH$_4$, N$_2$O, O$_2$, CFC-11, CFC-12, CFC-22, CCl$_4$ to non-zero?
  `ch4_val`                   0             If so, set value for CH$_4$
  `n2o_val`                   0             If so, set value for N$_2$O
  `o2_val`                    0             If so, set value for O$_2$
  `cfc{1112,22}_val`          0             If so, set value for CFC-{11,12,22}
  `ccl4_val`                  0             If so, set value for CCl$_4$
  `h2o_lower_limit`           0.2ppm        Never use specific humidity values smaller than this in radiative transfer
  `do_fixed_water`            .false.       Use fixed value for specific humidity in radiative transfer?
  `fixed_water`               2ppm          If so, set value
  `fixed_water_pres`          100e2Pa       If so, above which pressure level?
  `fixed_water_lat`           90            If so, equatorward of which latitude?
  `do_zm_tracers`             .false.       Feed only the zonal mean of all absorbers to the radiative transfer
  `do_zm_rad`                 .false.       Only pass zonal mean radiative forcing to dynamics
  `solday`                    0             Sets day of the year to run perpetual simulations. Seasonal cycle if =0
  `equinox_day`               0.25          Fraction of the year defining ‘March’ equinox
  `dt_rad`                    0             Radiation time step \[s\]. Every time step if $<$ `dt_atmos`
  `store_intermediate_rad`    .true.        Keep radiative forcing constant or set to zero between radiative time steps?
  `do_rad_time_avg`           .true.        Compute zenith angle average over radiative time step or use instantaneous?
  `dt_rad_avg`                -1            Average zenith angle over this time step. =`dt_rad` if $<0$
  `lonstep`                   1             Only compute radiation at every nth longitudinal grid point
  `do_precip_albedo`          .false.       Link surface albedo to precipitation?
  `precip_albedo`             0.35          If so, set albedo for 100% precipitating grid boxes
  `precip_lat`                0             If so, only use poleward of this latitude \[deg\]
  `precip_albedo_mode`        ’full’        If so, use ’full’ (total), ’lscale’ (large scale), or ’conv’ (convective) precipitation for albedo calculation

Table \[RRTM\_default\] gives the default values for all needed RRTM
input parameters, with the parameter names as given by the original
source code. We refer the reader to RRTM documentation and source code
to learn the exact meanings of each input parameter. The parameter type
is defined as follows: ‘fixed’ means the value is hard coded and cannot
be changed without recompiling; ‘fixed input’ means that the variable
can be changed within a namelist at run time, but will not change during
the simulation; ‘interactive’ means that the value is re-evaluated at
each radiation time step before calling RRTM.

  Variable     Type          Default
  ------------ ------------- ------------
  `iaer`       fixed         0
  `h2vmr`      interactive   tracer
  `o3vmr`      fixed input   from file
  `co2vmr`     fixed input   300e-6
  `ch4vmr`     fixed input   0
  `n2ovmr`     fixed input   0
  `o2vmr`      fixed input   0
  `cfc11vmr`   fixed input   0
  `cfc12vmr`   fixed input   0
  `cfc22vmr`   fixed input   0
  `ccl4vmr`    fixed input   0
  `asdir`      fixed input   0.27
  `aldir`      fixed         $=$`asdir`
  `asdif`      fixed         $=$`asdir`
  `aldir`      fixed         $=$`asdir`
  `coszen`     interactive   computed
  `adjes`      fixed input   1
  `dyofyr`     fixed input   0
  `scon`       fixed input   1368.22
  `emis`       fixed         1
  `inflgsw`    fixed         0
  `iceflgsw`   fixed         0
  `liqflgsw`   fixed         0
  `cldfr`      fixed         0
  `taucld`     fixed         0
  `ssacld`     fixed         0
  `asmcld`     fixed         0
  `asmcld`     fixed         0
  `fsfcld`     fixed         0
  `cicewp`     fixed         0
  `cliqwp`     fixed         0
  `reice`      fixed         10
  `reliq`      fixed         10
  `tauaer`     fixed         0
  `ssaaer`     fixed         0
  `asmaer`     fixed         0
  `acaer`      fixed         0

  : Input parameter settings for both RRTM short and long wave modules.
  ‘fixed’ means hard coded, ‘fixed input’ is read from the input file
  but remains constant, and ‘interactive’ is updated before every
  radiation call.<span data-label="RRTM_default"></span>
