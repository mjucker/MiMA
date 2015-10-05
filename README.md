# MiMA
Model of an idealized Moist Atmosphere: Intermediate-complexity General Circulation Model with interactive water vapor and full radiation.

This model is based on the gray radiation model of [Frierson, Held, and Zurita-Gotor, JAS (2006)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3753.1).
In fact, it even includes that exact model with a namelist switch flag. The major development step of MiMA is the replacement of the gray radiation scheme with a full radiative transfer code. For maximum portability and generality, that radiative transfer code is the Rapid Radiative Transfer Model [RRTM](http://rtweb.aer.com/rrtm_frame.html), developed by AER, and described in the references below.


## Compiling

* Dependencies
  * The code has only been tested with Intel ifort and icc compilers. There is no reason why it should not compile with any other Fortran and C compiler, but the compiler flags will have to be adjusted. This description will assume ifort and icc are available.
  * MiMA reads and writes to netCDF, so netcdf needs to be installed on the system.
  * Being parallel, MPI needs to be there too.
  * The flags need to be adjusted in the bin/mkmf.template.* file of choice.

* Compilation flags: The relevant flags are defined in bin/mkmf.template.*, and might or might not use environment variables. For instance, netCDF libraries or debug flags could be read from environment variables for more dynamic compilation. The first thing to do is to create an appropriate mkmf.template.something, which contains the relevant flags. Look at some of the template files that are already there to get an idea how to set the flags.

* Compile script: A compilescript is provided in exp/compilescript.csh. Make sure to set the first variable, platform, to whatever name you gave the mkmf template in the previous step. In our example, set it to something. The output executable will be in exp/exec.platform/mima.x

* Adding files: If you work on your own version of MiMA, make sure every extension is in a new file, so as to not disturb the main branch and any other fork that might exist. When adding a source file, add the path to the file in exp/path_names, and it will be compiled the next time you run ./compilescript.csh.



## Test run

A small test run is defined in the exp/ directory.
* input.nml: This is the most important file, where all the input parameters within the various namelists of MiMA can be set. If a variable is not present in input.nml, it will be set to the (hard coded) default value. This file completely defines the simulation you are running.
* diag_table: A list of diagnostic outputs you would like to have in your output files. This does not change the simulation you are running, but simply lets you decide which variables you'd like to have in the output, how frequently you'd like the output, and whether the output should be averaged or instantaneous.
* field_table: A list of passive tracers you'd like to advect during the simulation. There are two types: grid or spectral tracers. To get the temporal evolution of the tracers (or the time average), add the name of the tracer as diagnostic output to diag_table.

MiMA will automatically look for input.nml, so it is run without any explicit input (i.e. ./mima.x is fine - don't try ./mima.x < input.nml).

The test run will be one 360-day year with the following parameters:
* no surface topography
* ozone from the file input/ozone_1990.nc
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

## Radiation options

MiMA includes the gray radiation scheme developed by Dargan Frierson ([Frierson, Held, Zurita-Gotor, JAS (2006)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3753.1) ). To switch between the radiation schemes, the flags do_grey_radiation, and do_rrtm_radiation in the namelist physics_driver_nml can be set accordingly (only one of them should be .true. of course). 

Theoretically, there is also the possibility of running the full AM2 radiation scheme, with the flag do_radiation in physics_driver_nml. However, this option will need a lot of input files for tracer concentration, which are not part of the MiMA repository. This option, although all the relevant files are present and being compiled, has never been tested, and should only be used with great caution.

## References

### MiMA

* v1.0: TBD


### Gray radiation model

* [Frierson, Held, Zurita-Gotor, JAS (2006)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3753.1)
* [Frierson, JAS (2007)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3935.1)
* [Frierson, Held, Zurita-Gotor, JAS (2007)](http://journals.ametsoc.org/doi/abs/10.1175/JAS3913.1)


### RRTM

* [Mlawer et al., JGR (1997)](http://doi.wiley.com/10.1029/97JD00237)
* [Iacono et al., JGR (2000)](http://doi.wiley.com/10.1029/2000JD900091)
* [Iacono et al., JGR (2008)](http://onlinelibrary.wiley.com/doi/10.1029/2008JD009944/abstract)
* [Clough et al., JQSRT (2005)](http://www.sciencedirect.com/science/article/pii/S0022407304002158)


## License

MiMA is distributed under a GNU GPLv3 license. That means you have permission to use, modify, and distribute the code, even for commercial use. However, you must make your code publicly available under the same license. See LICENSE.txt for more details.

AM2 is distributed under a GNU GPLv2 license. That means you have permission to use, modify, and distribute the code, even for commercial use. However, you must make your code publicly available under the same license.

RRTM/RRTMG: Copyright © 2002-2010, Atmospheric and Environmental Research, Inc. (AER, Inc.). This software
may be used, copied, or redistributed as long as it is not sold and this copyright notice is reproduced
on each copy made. This model is provided as is without any express or implied warranties.
