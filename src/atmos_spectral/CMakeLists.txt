set ( ATMOS_SPECTRAL_SOURCES
# driver/solo/atmosphere.f90
driver/coupled/mcm_mca_lsc.f90
driver/coupled/mcm_moist_processes.f90
driver/coupled/spectral_physics.f90
driver/coupled/atmosphere.f90
model/leapfrog.f90
model/global_integral.f90
model/water_borrowing.f90
model/every_step_diagnostics.f90
model/tracer_type.f90
model/fv_advection.f90
model/press_and_geopot.f90
model/implicit.f90
model/spectral_dynamics.f90
model/spectral_damping.f90
model/matrix_invert.f90
# init/spectral_init_cond.f90.axi
init/vert_coordinate.f90
init/spectral_initialize_fields.f90
# init/spectral_init_cond.f90.nonaxi
init/spectral_init_cond.f90
init/topog_regularization.f90
tools/spec_mpp.f90
tools/gauss_and_legendre.f90
tools/transforms.f90
tools/spherical.f90
tools/spherical_fourier.f90
tools/grid_fourier.f90
)

target_sources( mima.x PRIVATE ${ATMOS_SPECTRAL_SOURCES} )
