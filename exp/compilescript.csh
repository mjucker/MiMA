#!/bin/csh -f
#Minimal runscript for atmospheric dynamical cores
set echo 
#--------------------------------------------------------------------------------------------------------
# define variables
set platform  = nci                                      # A unique identifier for your platform
set npes      = $PBS_NCPUS                               # number of processors
set template  = $cwd/../bin/mkmf.template.$platform   # path to template for your platform
set mkmf      = $cwd/../bin/mkmf                      # path to executable mkmf
set sourcedir = $cwd/../src                           # path to directory containing model source code
set mppnccombine = $cwd/../bin/mppnccombine.$platform # path to executable mppnccombine
#--------------------------------------------------------------------------------------------------------
set execdir   = $cwd/exec.$platform       # where code is compiled and executable is created
set pathnames = $cwd/path_names           # path to file containing list of source paths
#--------------------------------------------------------------------------------------------------------
set define_netcdf = 0
if ( ! $?NETCDF_INC ) then
  set define_netcdf = 1
else
  if ( "$NETCDF_INC" == "" ) then
    set define_netcdf = 1
  endif
endif
if ( $define_netcdf == 1 ) then
  set NETCDF_INC = `nc-config --fflags`
  set NETCDF_LIB = `nc-config --flibs`
  echo $NETCDF_INC
  echo $NETCDF_LIB
endif	
# compile mppnccombine.c, will be used only if $npes > 1
if ( ! -f $mppnccombine ) then
  icc -O -o $mppnccombine -I$NETCDF_INC -L$NETCDF_LIB $cwd/../postprocessing/mppnccombine.c -lnetcdf
endif
#--------------------------------------------------------------------------------------------------------
# setup directory structure
if ( ! -d $execdir ) mkdir $execdir
#--------------------------------------------------------------------------------------------------------
# compile the model code and create executable
cd $execdir
set cppDefs = "-Duse_libMPI -Duse_netCDF"
$mkmf -p mima.x -t $template -c "$cppDefs" -a $sourcedir $pathnames /usr/local/include $NETCDF_INC $sourcedir/shared/mpp/include $sourcedir/shared/include
make -f Makefile -j $npes
