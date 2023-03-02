[back to top](https://mjucker.github.io/MiMA)

# Getting started with MiMA

This model is based on the gray radiation model of [Frierson, Held, and Zurita-Gotor, JAS (2006)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3753.1).
In fact, it even includes that exact model with a namelist switch flag. The major development step of MiMA is the replacement of the gray radiation scheme with a full radiative transfer code. For maximum portability and generality, that radiative transfer code is the Rapid Radiative Transfer Model [RRTM](http://rtweb.aer.com/rrtm_frame.html), developed by AER, and described in the references below.

## Downloading source
You can download MiMA for free. However, we ask that you cite all relevant references given on the [front page](https://mjucker.github.com/MiMA/) with any publications that might result from its use.

Get the latest version from [GitHub](https://github.com/mjucker/MiMA/releases/latest).

## Compiling

* Dependencies
  * The code in its present form will only compile with Intel `ifort` and `icc` compilers. This description will assume `ifort` and `icc` are available.
  * MiMA reads and writes to `netCDF`, so `netCDF` needs to be installed on the system.
  * Being parallel, `MPI` needs to be there too.

* Build Systems
  * There are two build systems that may be used to compile MiMA, CMake, and `mkmf`.
  * Instructions for building with both of these are provided below.

### CMake
Building using CMake should follow the same process, regardless of the platform MiMA is being built on.

* Dependencies
  * In addition to the above list, building with CMake requires `cmake` to be installed on the system.

* To build MiMA using CMake after cloning the repository and navigating into it (e.g. via `cd MiMA/`) run the following commands:
  ```
  mkdir build
  cd build
  cmake ..
  make
  ```
  This takes you into the MiMA directory, creates a build directory, runs the CMake script `CMakeLists.txt` to generate a makefile for the system and then builds using the makefile.

* The output executable will be at `build/mima.x`

### mkmf
Building using mkmf requires the user to amend the `compilescript.csh` scipt as appropriate for the platform they are building on, and possibly defining a mkmf template for their platform.

* MiMA can be built using mkmf via the following steps:

  * Select the appropriate mkmf template file for your platform from `bin/`.
    These are of the form `bin/mkmf.template.$PLATFORM`, where `$PLATFORM` is typically slightly different on each machine as libraries may not be found at the same locations.

  * Adjust the relevant compilation flags in the `bin/mkmf.template.$PLATFORM` file of choice as appropriate. 
    * These may or may not use environment variables. For instance, `netCDF` libraries or debug flags could be read from environment variables for more dynamic compilation. The first thing to do is to create an appropriate `mkmf.template.something`, which contains the relevant flags. Look at some of the template files that are already there to get an idea how to set the flags.

  * Compile script: A compilescript is provided in `exp/compilescript.csh`. Make sure to set the first variable, `platform`, to match the `$PLATFORM` of the mkmf template in the previous step. In our example, set it to `something`.

  * MiMA can now be built with the following commands:
    ```
    cd exp
    ./compilescript.csh
    ```

* The output executable will be in `exp/exec.$PLATFORM/mima.x`


### Adding files to the build process

* If you work on your own version of MiMA, make sure every extension is in a new file, so as to not disturb the main branch and any other fork that might exist.
* When adding a source file you should:
  * add the file to the `CMakeLists.txt` file in its local directory,
  * add the path to the file in `exp/path_names`,
  To ensure that it will be compiled the next time you build using CMake or run `./compilescript.csh`.


## Test run

A small test run is defined in the input/ directory.
* `input.nml`: This is the most important file, where all the input parameters within the various namelists of MiMA can be set. If a variable is not present in `input.nml`, it will be set to the (hard coded) default value. This file completely defines the simulation you are running.
* `diag_table`: A list of diagnostic outputs you would like to have in your output files. This does not change the simulation you are running, but simply lets you decide which variables you'd like to have in the output, how frequently you'd like the output, and whether the output should be averaged or instantaneous.
* `field_table`: A list of passive tracers you'd like to advect during the simulation. There are two types: grid or spectral tracers. To get the temporal evolution of the tracers (or the time average), add the name of the tracer as diagnostic output to diag_table.

MiMA will automatically look for `input.nml`, so it is run without any explicit input (i.e. ``./mima.x`` is fine - don't try ``./mima.x < input.nml``).

The test run will be one 360-day year with the following parameters:
* no surface topography
* ozone from the file `INPUT/ozone_1990.nc` (which is also in the input/ directory)
* 300ppm CO2
* solar constant of 1360W/m2
* circular Earth-Sun orbit with 1UA radius
* NH solstice on December 30 (day 360)
* mixed layer ocean depth of 100m
* constant surface albedo of 0.27
* meridional Q flux of 30W/m2
* sponge layer Rayleigh friction above 50Pa
* Betts-Miller convection and large scale condensation
* RRTM radiation scheme

To run the test simulation, do the following:
```
EXECDIR=/PATH/TO/RUN/DIRECTORY
cp -r input/* $EXECDIR/
cp exp/exec.$PLATFORM/mima.x $EXECDIR/
cd $EXECDIR
mkdir RESTART
mpiexec -n $N_PROCS ./mima.x
CCOMB=/PATH/TO/MiMA/REPOSITORY/bin/mppnccombine.$PLATFORM
$CCOMB -r atmos_daily.nc atmos_daily.nc.*
$CCOMB -r atmos_avg.nc atmos_avg.nc.*
```
The last three lines make sure the indiviudal diagnostics files from each CPU are combined into one daily and one average file.

## Radiation options

By default, MiMA uses the RRTM radiation code. This is set by `do_rrtm_radiation = .true.` (default). There are, however, two more options for radiation, described below.

MiMA includes the gray radiation scheme developed by Dargan Frierson ([Frierson, Held, Zurita-Gotor, JAS (2006)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3753.1) ). To switch between the radiation schemes, the flags `do_grey_radiation`, and `do_rrtm_radiation` in the namelist `physics_driver_nml` can be set accordingly (only one of them should be `.true.` of course). 

Theoretically, there is also the possibility of running the full AM2 radiation scheme, with the flag `do_radiation in physics_driver_nml`. However, this option will need a lot of input files for tracer concentration, which are not part of the MiMA repository. This option, although all the relevant files are present and being compiled, has never been tested, and should only be used with great caution.

## Life cycle calculations

MiMA can be used to run life cycle experiments, as explore in Yamada and Pauluis (2017).  To specify the initial conditions, activate this flag in "spectral_dynamics_nml":

```fortran
specify_initial_conditions = .true.
```

Then a netcdf file containing the initial conditions for zonal wind, meridional wind, temperature, specific humidity, and surface pressure (ucomp, vcomp, temp, sphum, and ps, respectively) must be provided.  It should be at the resolution of the model.  It should be named initial_conditions.nc and placed in the INPUT/ directory where the model is executed.  Note that if you do not include a slight zonal perturbation, the model will maintain a zonally symmetric state, stuck to the unstabled fixed point.  There are different strategies for exciting zonal asymmetries.  You can add random noise, or focus in on a particular wavenumber, as detailed below.

A traditional life cycle is run with no forcing.  To shut off all diabatic processes, you must make these adjustments to the name list.  To turn off radiation and damping (except for hyperdiffusion), add these options to "physics_driver_nml"

```fortran
do_grey_radiation = .false.,
do_rrtm_radiation = .false.,
do_damping = .false. 
```

Then, to turn off any diabatic forcings at the lower boundary, in "surface_flux_nml" add these options:

```fortran
no_surface_momentum_flux  = .true.,
no_surface_moisture_flux  = .true.,
no_surface_heat_flux      = .true.,
no_surface_radiative_flux = .true. 
```

Lastly, the spectral dynamical core allows one to focus in on a particular wavenumber, as was done by Yamada and Pauluis (2017).  For example, to run a T170 resolution model, but enforce 6 fold symmetry (i.e., only capture instabilities at wave 6 and harmonics, use these options in "spectral_dynamics_nml":

* `lon_max                 = 128,`     [an ideal number for the Fourier transforms, close to 512/6]
* `lat_max                 = 256,`     [this grid corresponds to T170 resolution]
* `num_fourier             = 29,`     [this is approximately 170/6]
* `num_spherical           = 171,`     [this is always the T-resolution + 1]
* `fourier_inc             = 6,`       [this allows zonal waves 0, 6, 12, ...]

This trick allows you to run a higher resolution integration about 6 times faster.

Lastly, note that hyperdiffusion is still required for stability.  For the Yamada and Pauluis life cycles, these options were selected in "spectral_dynamics_nml":

```fortran
damping_option          = 'resolution_dependent',
damping_order           = 3,
damping_coeff           = 6.94444444e-5,
damping_order_vor       = 3,
damping_order_div       = 3,
damping_coeff_vor       = 6.94444444e-5,
damping_coeff_div       = 6.94444444e-5,
```

### Reference:
[Yamada, R., and O. Pauluis, 2017: Wave-mean-flow interactions in moist baroclinic lifecycles. J. Atmo. Sci., 74, 2143-2162, doi:10.1175/JAS-D-16-0329.1](https://doi.org/10.1175/JAS-D-16-0329.1).


