set ( ICE_PARAM_SOURCES
ocean_rough.f90
)

target_sources( mima.x PRIVATE ${ICE_PARAM_SOURCES} )
