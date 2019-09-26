#!/bin/csh -f
#Minimal runscript for atmospheric dynamical cores
set echo 
#--------------------------------------------------------------------------------------------------------
# define variables
set platform  = nci                                   # A unique identifier for your platform
set npes      = $PBS_NCPUS                            # number of processors
set template  = $cwd/../bin/mkmf.template.$platform   # path to template for your platform
set mkmf      = $cwd/../bin/mkmf                      # path to executable mkmf
set sourcedir = $cwd/../src                           # path to directory containing model source code
set mppnccombine = $cwd/../bin/mppnccombine.$platform # path to executable mppnccombine
#--------------------------------------------------------------------------------------------------------
set execdir   = $cwd/exec.$platform       # where code is compiled and executable is created
set workdir   = $cwd/workdir              # where model is run and model output is produced
set pathnames = $cwd/path_names           # path to file containing list of source paths
set namelist  = $cwd/namelists            # path to namelist file
set diagtable = $cwd/diag_table           # path to diagnositics table
set fieldtable = $cwd/field_table         # path to field table (specifies tracers)
#--------------------------------------------------------------------------------------------------------
# compile mppnccombine.c, will be used only if $npes > 1
set add_flag = 0
if (! $?NETCDF_INC || ! $?NETCDF_LIB) then
    set add_flag = 1
else
    if ("$NETCDF_INC" == "" || "$NETCDF_LIB" == "") then
        set add_flag = 1
    endif
endif	
if ($add_flag) then
    set NETCDF_INC = `nc-config --includedir`
    set NETCDF_LIB = `nc-config --libs`
    echo "NETCDF_INC: "$NETCDF_INC
    echo "NETCDF_LIB: "$NETCDF_LIB
endif
if ( ! -f $mppnccombine ) then
  icc -O -o $mppnccombine -I$NETCDF_INC $NETCDF_LIB $cwd/../postprocessing/mppnccombine.c 
endif
if ($add_flag) then
    set NETCDF_INC = `nc-config --fflags`
    set NETCDF_LIB = `nc-config --flibs`
endif
if ($add_flag) then
    set NETCDF_INC = `nc-config --fflags`
    set NETCDF_LIB = `nc-config --flibs`
endif
#--------------------------------------------------------------------------------------------------------
# setup directory structure
if ( ! -d $execdir ) mkdir $execdir
if ( -e $workdir ) then
  echo "ERROR: Existing workdir may contaminate run.  Move or remove $workdir and try again."
  exit 1
endif
#--------------------------------------------------------------------------------------------------------
# compile the model code and create executable
cd $execdir
set cppDefs = "-Duse_libMPI -Duse_netCDF"
$mkmf -p mima.x -t $template -c "$cppDefs" -a $sourcedir $pathnames /usr/local/include $NETCDF_INC $sourcedir/shared/mpp/include $sourcedir/shared/include
make -f Makefile -j $npes
