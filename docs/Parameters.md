# Parameter settings

This section shows (hopefully) all parameters and their default values. For details about the physical meaning of these parameters, please refer to the main MiMA reference paper.

General
-------

Most general parameters are set in `coupler/coupler_main.f90` and `coupler_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 current_date | 0001,1,1,0,0,0 | Start from year 1
 calendar | 'thirty_day' | Use 12 30-day months

Dynamics
--------

The dynamical core is set up in `atmos_spectral/model/spectral_dynamics.f90` and the namelist `spectral_dynamics_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
   damping_order           | 4 | 8<sup>th</sup> order numerical diffusion
   do_water_correction     | .true. | make sure water mass is conserved in the advection step
   water_correction_limit  | 200.e2 | correct water mass only below 200hPa. removes stratospheric sink
   vert_advect_uv          | 'second_centered' | second order vertical advection scheme
   vert_advect_t           | 'second_centered' | second order vertical advection scheme
   robert_coeff            | .03 | Robertson coefficient for implicit time stepping
   lon_max                 | 128 | T42 resolution
   lat_max                 | 64  | T42 resolution
   num_fourier             | 42  | T42 resolution
   num_spherical           | 43  | T42 resolution
   fourier_inc             | 1   | T42 resolution
   initial_sphum           | 2.e-6 | constant specific humidity initial condition

Initial conditions are set in `atmos_spectral/init/spectral_init_cond.f90` and the namelist `spectral_init_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 initial_temperature | 264 | Initial temperature of isothermal atmosphere in [K]


Mixed layer
-----------

All changes to the mixed layer ocean as compared to FHZ06 are discussed
in the main text. Table \[Surface\_default\] gives the default values
for all newly introduced input parameters, and those that might have
slightly changed meaning.

These parameters are set in `coupler/simple_surface.f90`.

  Variable           | Default Value | Meaning
  :--- | :---: | :---
  heat_capacity   | 4e-8 | 100m mixed layer depth
  land_capacity   | -1 | same as `heat_capacity`
  trop_capacity   | -1 | same as `heat_capacity`
  trop_cap_limit  | 15 | Poleward boundary for `trop_capacity`
  heat_cap_limit  | 60 | Equatorward boundary for `heat_capacity`
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
  slandlon        | (0,0,0,0,0,0,0,0,0,0) | if so, give land patches start longitude
  elandlon        | -(1,1,1,1,1,1,1,1,1,1)| if so, give land patches end longitude
  slandlat        | (0,0,0,0,0,0,0,0,0,0) | if so, give land patches start latitude
  elandlat        | -(1,1,1,1,1,1,1,1,1,1)| if so, gie land patches end longitude
 

Moist processes
---------------

In contrast to the initial *Frierson et al (2007)* setup, we turn off moist convective adjustment by default, and turn on Betts-Miller convection and use the alternative definition of specific humidity.

The moist processes parameters are set in `atmos_param/moist_processes/moist_processes.f90` and the namelist `moist_processes_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 do_mca | .false. | Do moist convective adjustment
 do_lsc | .true. | Do large scale condensation
 do_bm  | .true. | Do Betts-Miller convetion
 do_df_stuff | .true. | When true, specific humidity = (rdgas/rvgas)*esat/pressure
 
Betts-Miller
------------

We use a simple Betts-Miller convection
scheme ($\tau_\mathrm{SBM} = $2hrs, RH$_\mathrm{SBM}=0.7$) with the
‘shallower’ shallow convection as described by [@Frierson2007a].
Precipitation only falls out if all layers below the condensation level
are saturated, otherwise, the condensate is re-evaporated.

The Betts-Miller parameters are set in `atmos_param/betts_miller/betts_miller.f90` and the namelist `betts_miller_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 tau_bm | 7200 | relaxation time scale in [s]
 rhbm   | 0.7  | relative humidity to relax to
 do_sim | .false. | don't adjust time scales to make precipitation always continuous
 do_shallower | .true. | use shallower convection scheme
 
Moist convection
----------------

Moist convection parameters are in `atmos_param/moist_conv/moist_conv.f90` and the namelist `moist_conv_nml`. This is only touched to make sure moisture handling is consistent accross namelists.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 use_df_stuff | .true. | Make everything consistent with above `do_df_stuff`

Large Scale Condensation
------------------------

We want to re-evaporate outfalling precipitation if any of the layers below are sub-saturated. Parameters are described in `atmos_param/lscale_cond/lscale_cond.f90` and the namelist `lscale_cond_nml`.

 Variable | Default Value | Meaning
 :--- | :---: | :---
 do_evap | .true. | re-evaporate in below sub-saturated layers (if any)
 use_df_stuff | .true. | Make everything consistent


Boundary conditions
-------------------

The lower boundary conditions are unchanged with respect to FHZ06, i.e.
the same Monin-Obukhov theory is applied for determining surface drag,
setting surface fluxes and boundary layer height. The top of the
atmosphere, however, needed to be updated due to fully resolved
stratosphere, and in particular the polar vortex. For simplicity, a
simple Rayleigh drag is included as a sponge layer to parameterize
gravity wave drag. This is the exact same sponge layer as in the
Newtonian cooling codes of
[@Polvani2002; @Kushner2004; @Gerber2009; @Jucker2013; @Jucker2014], and
is applied above 0.5hPa with a time scale of 0.5day at the top layer.
The possibility to use a more accurate gravity wave drag scheme is
already implemented, but has not been thoroughly tested. This will be
subject of future code development.

Radiation
---------

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
