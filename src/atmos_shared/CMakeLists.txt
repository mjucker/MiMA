set ( ATMOS_SHARED_SOURCES
tracer_driver/atmos_convection_tracer.f90
tracer_driver/atmos_tracer_driver.f90
tracer_driver/atmos_radon.f90
tracer_driver/atmos_carbon_aerosol.f90
tracer_driver/atmos_sulfur_hex.f90
tracer_driver/atmos_tracer_utilities.f90
interpolator/interpolator.F90
vert_advection/vert_advection.f90
)

target_sources ( mima.x PRIVATE ${ATMOS_SHARED_SOURCES} )
