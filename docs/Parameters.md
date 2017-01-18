# Parameter settings

MiMA is designed with the idea of including as much detail and
sophistication as needed, but no more. This definition for the
complexity of the model depends on the exact experiment one wants to
undertake. To add more flexibility, a few somewhat more sophisticated
options are included in the model, which are not used for the here
presented studies of the cold point.

Mixed layer depth
-----------------

Instead of having a constant heat capacity $C_O$ in Equation
(\[mixedLayer\]) everywhere, we introduced the possibility to have a
different heat capacity in the tropics than in the extratropics, with a
linear interpolation in-between. This is motivated by the fact that the
real ocean mixed layer is shallower in the tropics than higher latitudes
\[SOME REFERENCE HERE?\], and the heat capacity is the only quantity in
the surface parametrization that is directly connected to mixed layer
depth. Mathematically, the ocean heat capacity is evaluated as follows:
$$\label{heatCapacitySea}
    C_O^s(\phi) =
        \begin{cases}   C_{t} & ,|\phi| < \phi_{t} \\
                        C_{t} + \Delta(\phi)(C_{p}-C_{t}) &,\phi_{t}\leq |\phi| \leq \phi_{p} \\
                        C_{p} & ,|\phi| > \phi_{p},
        \end{cases}$$ where
$\Delta(\phi) = (\phi-\phi_{t})/(\phi_{p}-\phi_{t})$. Typically,
$\phi_{t}=15^\circ$ and $\phi_{p}=60^\circ$. One of the effects of a
variable tropical mixed layer is that the tropical precipitation maximum
(i.e. the Intertropical Convergence Zone ITCZ) has a more (shallow
tropics) or less (deep tropics) pronounced seasonal cycle in terms of
meridional position.

A second variation to the heat capacity setup described in the main text
concerns the definition of the ‘land’ sections at the surface. Instead
of defining rectangular patches of land, it is possible to link the heat
capacity of the mixed layer to the geopotential height of the surface,
if there is surface topography. With this choice, all regions of the
surface where the geopotential is larger than 10m$^2$/s$^2$ will take
the value for the heat capacity defined for ‘land’.

Meridional albedo profile
-------------------------

Instead of using a constant surface albedo, in some situations it might
be helpful to define an albedo which varies with latitude. The most
obvious example is the positioning and strength of the eddy-driven jet,
which can be displaced and/or strengthened by changing the meridional
temperature gradients though the albedo effect. Several basic options
are present in MiMA:

-   Discontinuous step: The albedo can change abruptly between two
    different values at a given latitude. This can either be done
    symmetrically in both hemispheres or only in one hemisphere

-   Continuous profile: Vary the albedo according to the form
    $$\alpha(\phi) = \alpha_0 + (\alpha_1-\alpha_0)\left(|\phi|/90\right)^p,$$
    with the input albedos $\alpha_0 and \alpha_1$, the order $p$, and
    the latitude $\phi$ in degrees.

-   Continuous step: For ‘best of both worlds’, there is an option to
    set the profile according to
    $$\alpha(\phi) = \alpha_0 + (\alpha_1-\alpha_0)*0.5*[1+\tanh\left((\phi-\phi_0)/\Delta\phi\right)],$$
    with the input albedos $\alpha_0$ and $\alpha_1$, the center
    latitude $\phi_0$, and the width $\Delta\phi$ .

‘Cloud’ albedo
--------------

As mentioned by [@Merlis2013], in the real world the meridional
structure of cloud properties plays some part in setting up the
meridional temperature gradient through radiative effects. But rather
than prescribing a fixed cloud fraction and liquid water content as
these authors, we implemented a scheme where the surface albedo can be
linked to precipitation.

The scheme first counts the number of time steps when precipitation was
non-zero in each horizontal grid box, either total precipitation, large
scale condensation, or precipitation due to convection. Whenever the
radiative scheme is called, this number is divided by the total number
of time steps since the last radiation call. This gives the ratio of
precipitating time steps for each grid box at the surface,
$r_p(\lambda,\phi) \in [0,1]$. To avoid the risk of instabilities, the
albedo is changed symmetrically in the zonal direction, i.e. only the
zonal mean $\overline{r}_p(\phi)$ is used to compute the surface albedo
$\alpha(\phi)$ as $$\label{albedo}
    \alpha(\phi) = \alpha_c + (\alpha_p - \alpha_c)\overline{r}_p(\phi),$$
where $\alpha_c$ is the constant part of the albedo, related to the
absence of precipitation, and $\alpha_p$ is the maximum albedo linked to
precipitation (the ‘cloud’ albedo). Typically, $\alpha_c=0.1$ and
$\alpha_p=0.35$.

This way of including one of the most basic effects of clouds can have
an impact on the dynamics of the model, in particular in the tropics and
around the midlatitude jet streams. It also allows for a more dynamic
positioning of this cloud albedo effect, as compared to a fixed
meridional profile.

Radiative transfer
------------------

#### Perpetual simulations

Whereas all simulations presented in this work are based on full
seasonal cycle simulations, it is possible to run MiMA in perpetual
simulation mode. The optional input variable `solday` can be set to the
day of the year of the perpetual simulation, and it forces the zenith
angle of the sun to always correspond to that day of the year. We note
that this option has not been rigorously tested by the authors, and it
is possible that for some choices of perpetual simulations the
temperatures reach unphysical values, for instance in the polar night
due to constant cooling.

#### Time stepping

In addition to setting the radiative time step, it is possible to use an
average instead of an instantaneous zenith angle for the short wave
radiative forcing. By default, the zenith angle averaged over the
radiation time step is used, but the averaging time can be adjusted
independently, or turned off completely. If the radiation time step is
longer than the dynamics time step (usually the case), the radiative
forcing is kept constant over the intermediate dynamics time steps. This
can be turned off and radiative forcing is only added at every radiation
time step.

#### Spatial resolution

The radiative transfer calculations take a considerable amount of time.
In order to speed up simulations, it is possible to compute the
radiation only at every $n$th longitudinal grid point. The radiative
forcing is then interpolated linearly back onto the full longitudinal
grid. By default, this variable, `lonstep`, is set to 4, i.e. only every
fourth longitudinal grid point is used to compute the radiative
transfer. There is no such possibility for latitude, as the
parallelization of the code is done in meridional direction, and
typically only two meridional grid points are present in each compute
node.

\[modelparamsapp\]

All input files used for the here presented simulations are available
for download \[DATA REFERENCE\]. As described in the main text, the
surface drag, fluxes, mixed layer, and water vapor (thermo)dynamics are
as described in FHZ06 and [@Frierson2007a]. We will restrict ourselves
here to a description of the general setup, and give exact values for
important newly introduced variables.

Dynamics
--------

In the dynamical core, we use eighth order numerical diffusion, mass,
energy, and water correction schemes. The latter is only applied for
pressure surfaces below 200hPa, i.e. in the troposphere. Vertical
advection uses a second order centered scheme, and implicit time
stepping is used with $\alpha=0.5$ and a Robertson coefficient of 0.03.
The vertical coordinates default to an uneven sigma coordinate, spanning
11 scale heights, but in this work we used an input profile with highly
increased resolution between 90 and 130hPa. Spinup is achieved starting
from an isothermal atmosphere at 264K and at rest. Time stepping is
10min here (T42 resolution), but can be relaxed to 15min for more
regular vertical grids.

Mixed layer
-----------

All changes to the mixed layer ocean as compared to FHZ06 are discussed
in the main text. Table \[Surface\_default\] gives the default values
for all newly introduced input parameters, and those that might have
slightly changed meaning.

  Variable           Default
  ------------------ ------------------------------
  `heat_capacity`    4e-8
  `land_capacity`    -1 \[i.e. `=heat_capacity`\]
  `trop_capacity`    -1 \[i.e. `=heat_capacity`\]
  `trop_cap_limit`   15
  `heat_cap_limit`   60
  `albedo_choice`    1 \[i.e. constant\]
  `const_albedo`     0.30
  `albedo_exp`       2
  `do_qflux`         .false.
  `qflux_width`      16
  `qflux_amp`        30
  `do_warmpool`      .false.
  `warmpool_amp`     30
  `warmpool_phase`   0
  `warmpool_centr`   0
  `warmpool_k`       1
  `land_option`      ’none’
  `slandlon`         (0,0,0,0,0,0,0,0,0,0)
  `elandlon`         -(1,1,1,1,1,1,1,1,1,1)
  `slandlat`         (0,0,0,0,0,0,0,0,0,0)
  `elandlat`         -(1,1,1,1,1,1,1,1,1,1)

  : Default values for surface mixed layer. Only quantities that have
  been added or changed meaning compared to FHZ06 are included.<span
  data-label="Surface_default"></span>

Moist processes
---------------

We use large scale condensation and a simple Betts-Miller convection
scheme ($\tau_\mathrm{SBM} = $2hrs, RH$_\mathrm{SBM}=0.7$) with the
‘shallower’ shallow convection as described by [@Frierson2007a].
Precipitation only falls out if all layers below the condensation level
are saturated, otherwise, the condensate is re-evaporated.

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
