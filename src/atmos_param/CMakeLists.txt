set ( ATMOS_PARAM_SOURCES
cloud_zonal/null/cloud_zonal.f90
vert_diff_driver/vert_diff_driver.f90
moist_conv/moist_conv.f90
# moist_conv/null/moist_conv.f90
moist_processes/moist_processes.f90
damping_driver/damping_driver.f90
cloud_rad/cloud_rad.f90
cg_drag/cg_drag.f90
ras/ras.f90
diag_integral/diag_integral.f90
edt/edt.f90
# edt/null/edt.f90
# my25_turb/null/my25_turb.f90
my25_turb/my25_turb.f90
diag_cloud/null/diag_cloud.f90
vert_diff/vert_diff.f90
# shallow_conv/shallow_conv.f90
shallow_conv/null/shallow_conv.f90
grey_radiation/grey_radiation.f90
vert_turb_driver/vert_turb_driver.f90
qflux/qflux.f90
# donner_deep/null/donner_deep.f90
donner_deep/donner_deep.f90
entrain/entrain.f90
cu_mo_trans/cu_mo_trans.f90
mg_drag/mg_drag.f90
cloud_generator/betaDistribution.f90
cloud_generator/null/cloud_generator.f90
topo_drag/null/topo_drag.f90
# diffusivity/null/diffusivity.f90
diffusivity/diffusivity.f90
physics_driver/physics_driver.f90
rh_clouds/null/rh_clouds.f90
betts_miller/bm_omp.f90
betts_miller/betts_miller.f90
betts_miller/bm_massflux.f90
monin_obukhov/monin_obukhov.f90
lscale_cond/lscale_cond.f90
# lscale_cond/null/lscale_cond.f90
radiation_driver/radiation_driver.f90
local_heating/local_heating.f90
astronomy/astronomy.f90
sea_esf_rad/sea_esf_rad.f90
sea_esf_rad/rad_output_file.f90
sea_esf_rad/microphys_rad.f90
sea_esf_rad/microphys_cloud.f90
# FIXME this is the overlay thing, this file is replaced by the above e.g. sea_esf_rad/null/specified_clouds_W.f90
# sea_esf_rad/standalone_clouds.f90
sea_esf_rad/bulkphys_rad.f90
sea_esf_rad/longwave_fluxes.f90
sea_esf_rad/optical_path.f90
sea_esf_rad/aerosolrad_package.f90
sea_esf_rad/longwave_clouds.f90
sea_esf_rad/gas_tf.f90
sea_esf_rad/sealw99.f90
sea_esf_rad/rh_based_clouds.f90
sea_esf_rad/diag_clouds_W.f90
sea_esf_rad/isccp_clouds.f90
# sea_esf_rad/specified_clouds_W.f90
sea_esf_rad/donner_deep_clouds_W.f90
sea_esf_rad/radiation_diag.f90
sea_esf_rad/longwave_tables.f90
sea_esf_rad/rad_utilities.f90
sea_esf_rad/lhsw_driver.f90
sea_esf_rad/ozone.f90
sea_esf_rad/cloudrad_diagnostics.f90
sea_esf_rad/radiative_gases.f90
sea_esf_rad/esfsw_driver.f90
sea_esf_rad/zetac_clouds_W.f90
sea_esf_rad/longwave_params.f90
sea_esf_rad/aerosol.F90
# sea_esf_rad/original_fms_rad.f90
sea_esf_rad/cloud_spec.f90
sea_esf_rad/null/standalone_clouds.f90
# sea_esf_rad/null/bulkphys_rad.f90
# sea_esf_rad/null/rh_based_clouds.f90
# sea_esf_rad/null/diag_clouds_W.f90
# FIXME this is the overlay thing, this file is replaced by the above e.g. sea_esf_rad/null/specified_clouds_W.f90
sea_esf_rad/null/specified_clouds_W.f90
# sea_esf_rad/null/donner_deep_clouds_W.f90
# sea_esf_rad/null/lhsw_driver.f90
sea_esf_rad/null/original_fms_rad.f90
# sea_esf_rad/null/mgrp_prscr_clds.f90
sea_esf_rad/shortwave_driver.f90
sea_esf_rad/strat_clouds_W.f90
sea_esf_rad/lw_gases_stdtf.f90
sea_esf_rad/mgrp_prscr_clds.f90
sea_esf_rad/longwave_driver.f90
sea_esf_rad/cloudrad_package.f90
sea_esf_rad/esfsw_parameters.f90
stable_bl_turb/stable_bl_turb.f90
dry_adj/dry_adj.f90
# dry_adj/null/dry_adj.f90
rrtm_radiation/astro.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/mcica_subcol_gen_sw.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_spcvmc.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_cldprop.f90
# rrtm_radiation/rrtmg_sw/gcm_model/src/mcica_random_numbers.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_rad.nomcica.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_spcvrt.f90
# rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_rad.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_cldprmc.f90
# rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_read_nc.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_setcoef.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_init.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_reftra.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_vrtqdr.f90
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_taumol.f90
# Long Time to Compile so special flags in main CMakeLists
rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_k_g.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg22.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg20.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg26.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg23.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_vsn.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg24.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_ref.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg19.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_cld.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg25.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg28.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_wvn.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg21.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg27.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_con.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/parrrsw.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_aer.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg16.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_ncpar.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_tbl.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg17.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg18.f90
# rrtm_radiation/rrtmg_sw/gcm_model/modules/parkind.f90
rrtm_radiation/rrtmg_sw/gcm_model/modules/rrsw_kg29.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_cldprmc.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_setcoef.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/mcica_subcol_gen_lw.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_init.f90
# Not required if rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_k_g.f90 used
# rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_read_nc.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_taumol.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_rtrnmc.f90
# Long Time to Compile so special flags in main CMakeLists
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_k_g.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/mcica_random_numbers.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_cldprop.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_rtrn.f90
# rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_rad.f90
rrtm_radiation/rrtm_radiation.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_rtrnmr.f90
rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_rad.nomcica.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg03.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_con.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg14.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg11.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg04.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_wvn.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_ncpar.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_vsn.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg15.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/parrrtm.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg05.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_ref.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg16.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg09.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg06.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg10.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg12.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg01.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/parkind.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg08.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_cld.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg13.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_tbl.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg02.f90
rrtm_radiation/rrtmg_lw/gcm_model/modules/rrlw_kg07.f90
diag_cloud_rad/null/diag_cloud_rad.f90
strat_cloud/strat_cloud.f90
)

target_sources ( mima.x PRIVATE ${ATMOS_PARAM_SOURCES} )
