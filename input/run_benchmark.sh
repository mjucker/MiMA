ulimit -s unlimited

platform=nci
plevel=plevel.sh
NCPUS=32

mkdir ../benchmark
cp -r * ../benchmark/
cd ../benchmark
cp input_benchmark.nml input.nml
mkdir RESTART
cp ../exp/exec.${platform}/mima.x .
cp ../bin/mppnccombine.${platform} .
cp ../postprocessing/plevel_interpolation/scripts/plevel.sh .


#run 30 years
for y in {0..29}
do
    if [[ $y -lt 10 ]]
    then
	yy=0$y
    else
	yy=$y
    fi
    # run mima for one year
    mpiexec -np $NCPUS ./mima.x > ${yy}.out.txt
    # postproc
    mkdir ${yy}
    for file in *.nc.0000
    do
	baseName="${file%%.*}"
	fileName=${yy}.${baseName}.nc
	if [ -f ${fileName} ]
	then
            rm ${fileName}
	fi
	./mppnccombine.${platform} ${fileName} ${baseName}.nc.????
	if [ $? -eq 0 ]
	then
            rm ${baseName}.nc.????
	fi
    done
    for fileName in *.nc
    do
        sh ${plevel} -a -i ${fileName} slp hght
        if [ $? -ne 0 ]
        then
	    sh plevel.sh -a -i ${fileName}
        fi
        mv plevel.nc ${fileName/.nc/.plev.nc}
    done
    tar cvf ${yy}.restart.tar RESTART/*
    mv ${yy}.* ${yy}/
    cp RESTART/* INPUT/
done

#python benchmark_analysis.py
