FROM robcking/eddy-builder
ENV FFLAGS="-O2 -i4 -r8 -g -I/opt/netcdf-fortran/include -I/opt/netcdf-c/include -I/opt/hdf5/include -I/include -cpp"
ENV CFLAGS="-g -I/opt/netcdf-fortran/include -I/opt/netcdf-c/include -I/opt/hdf5/include/ -I/include"
ENV LDFLAGS="-shared-intel -L/opt/netcdf-fortran/lib -lnetcdff -L/opt/netcdf-c/lib -lnetcdf"
ENV cppDefs="-Duse_libMPI -Duse_netCDF -DgFortran"
WORKDIR /opt/MiMA
COPY ./src ./src
COPY ./bin ./bin
COPY ./exp ./exp
COPY ./postprocessing ./postprocessing
WORKDIR /opt/MiMA/build
RUN icx -g -O2 -o mppnccombine.docker -I/opt/netcdf-c/include -I/opt/netcdf-fortran/include -I/opt/hdf5/include -L/opt/netcdf-c/lib -L/opt/hdf5/lib -L/opt/netcdf-fortran/lib -lnetcdf -lnetcdff /opt/MiMA/postprocessing/mppnccombine.c
RUN chmod +x mppnccombine.docker

RUN /opt/MiMA/bin/mkmf -p mima.x -t "/opt/MiMA/bin/mkmf.template.singularity" -c "${cppDefs}" -a /opt/MiMA/src /opt/MiMA/exp/path_names ${NETCDF_INC} ${NETCDF_LIB} ${NETCDF_FORTRAN_INC} ${NETCDF_FORTRAN_LIB} ${HDF5_INC} ${HDF5_LIB} /opt/MiMA/src/shared/mpp/include /opt/MiMA/src/shared/include
RUN make -f Makefile clean
RUN make -f Makefile
RUN chmod +x /opt/MiMA/build/mima.x
WORKDIR /root
CMD ulimit -s unlimited && mpiexec -n 1 /opt/MiMA/build/mima.x
