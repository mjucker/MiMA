#!/bin/bash

path_to_inputfiles=$HOME/path/to/input/files

rundir=/path/to/run/directory

npes=32
num_executions=2

cp -r $path_to_inputfiles/* $rundir/

cd $rundir
mkdir RESTART

n=0
while [ $n -lt $num_executions ]
do
  index=`printf %04d ${n%*} ${n##*}`
  mpiexec -n $npes ./mima.x > out.${index}.txt
  mppnccombine -r ${index}.atmos_daily.nc atmos_daily.nc.????
  mppnccombine -r ${index}.atmos_avg.nc atmos_avg.nc.????
  cp RESTART/*res* INPUT/
  let n=$n+1
done
