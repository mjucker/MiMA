set ( COUPLER_SOURCES
simple_surface.f90
coupler_main.f90
# No longer seems to be required in new MiMA?
#flux_exchange.f90
surface_flux.f90
)

# NB this is where the executable is defined.
target_sources ( mima.x PRIVATE ${COUPLER_SOURCES} )
