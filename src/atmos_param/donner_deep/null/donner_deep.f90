                       module donner_deep_mod

use time_manager_mod,       only: time_manager_init, time_type, &
                                  set_date, get_time,   &
                                  get_calendar_type, &
                                  operator(>=), operator (<)
use diag_manager_mod,       only: register_diag_field, send_data, &
                                  diag_manager_init
use sat_vapor_pres_mod,     only: lookup_es, sat_vapor_pres_init
use fms_mod,                only: fms_init, mpp_pe, mpp_root_pe,  &
                                  file_exist,  check_nml_error,  &
                                  error_mesg, FATAL, WARNING, NOTE,  &
                                  close_file, open_namelist_file,    &
                                  stdlog, write_version_number,  &
                                  read_data, write_data,    &
                                  open_restart_file
use constants_mod,          only: constants_init, DENS_H2O, RDGAS,   &
                                  CP_AIR, pie=>PI
use column_diagnostics_mod, only: column_diagnostics_init, &
                                  initialize_diagnostic_columns, &
                                  column_diagnostics_header, &
                                  close_column_diagnostics_units


implicit none
private

!--------------------------------------------------------------------
!         module to compute the effects of deep convection
!
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: donner_deep.f90,v 11.0 2004/09/28 19:16:10 fms Exp $'
character(len=128)  :: tagname =  '$Name: lima $'


!--------------------------------------------------------------------
!---interfaces------

public   &
        donner_deep_init, donner_deep, donner_deep_avg, donner_deep_end


logical :: module_is_initialized =  &
                              .false.   ! has this module been 
                                        ! initialized ?
!---------------------------------------------------------------------



                        contains




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################

subroutine donner_deep_init (lonb, latb, pref, axes, Time, &
                             tracers_in_donner)

!---------------------------------------------------------------------
!    donner_deep_init is the constructor for donner_deep_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real,    dimension(:), intent(in)      :: lonb, latb, pref
integer, dimension(4), intent(in)      :: axes
type(time_type),       intent(in)      :: Time
logical, dimension(:), intent(in), optional  :: tracers_in_donner

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      lonb     array of model longitudes on cell boundaries [ radians ]
!      latb     array of model latitudes on cell boundaries [ radians ]
!      pref     array of reference pressures at full levels (plus 
!               surface value at nlev+1), based on 1013.25 hPa pstar
!               [ Pa ]
!      axes     data axes for diagnostics
!      Time     current time [ time_type ]
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!--------------------------------------------------------------------
!    set flag to indicate that donner_deep_mod has been initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('donner_deep_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_init



!###################################################################

subroutine donner_deep (is, ie, js, je, temp, mixing_ratio, pfull,  &
                        phalf, omega, dt, land, Time, &
                        ttnd, qtnd, precip, ahuco, qrat,  &
                        kbot, cf, qlin, qiin, &
                        delta_qa, delta_ql, delta_qi, mtot, &
                        tracers, qtrceme)
                        
!-------------------------------------------------------------------
!   donner_deep is the prognostic subroutine of donner_deep_mod. it 
!   returns the tendencies of temperature and mixing ratio (and changes
!   in cloudwater, cloudice and cloud area when run with a prognostic 
!   cloud scheme), and the precipitation and the mass flux (if needed)
!   produced by deep convection, as determined by this parameterization.
!-------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                 intent(in)            :: is, ie, js, je
real, dimension(:,:,:),  intent(in)            :: temp, mixing_ratio
real, dimension(:,:,:),  intent(in)            :: pfull, phalf
real, dimension(:,:,:),  intent(in)            :: omega
real,                    intent(in)            :: dt
real, dimension(:,:),    intent(in)            :: land
type(time_type),         intent(in)            :: Time
real, dimension(:,:,:),  intent(out)           :: ttnd, qtnd
real, dimension(:,:),    intent(out)           :: precip      
real, dimension(:,:,:),  intent(inout)         :: ahuco, qrat
real, dimension(:,:,:,:), intent(in), optional :: tracers
real, dimension(:,:,:,:), intent(out),optional :: qtrceme
integer, dimension(:,:), intent(in),  optional :: kbot
real, dimension(:,:,:),  intent(in),  optional :: cf
real, dimension(:,:,:),  intent(in),  optional :: qlin,qiin
real, dimension(:,:,:),  intent(out), optional :: delta_qa, delta_ql, &
                                                  delta_qi
real, dimension(:,:,:),  intent(out), optional :: mtot
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   mixing ratio field at model levels [ kg/kg ]
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field at half-levels 1:nlev+1  [ Pa ]
!     omega          model omega field [ Pa/sec ]
!     dt             physics time step [ sec ]
!     land           fraction of surface (grid box) covered by land
!                    [ nondimensional ]
!     Time           current time (time_type)
!
!   intent(out) variables:
!
!     ttnd           time tendency of temperature due to deep convection
!                    [ deg K / sec ]
!     qtnd           time tendency of mixing ratio due to deep
!                    convection [ (kg/kg) / sec ]
!     precip         precipitation generated by deep convection
!                    [ kg / m**2 ]
!     ahuco          cell + meso fraction for specific humidity
!                    (mesoscale downdraft does not contain cloud
!                    water, so this fraction should not be used
!                    for radiation)
!                    index 1 nearest ground
!     qrat           ratio of large-scale specific humidity to
!                    specific humidity in environment outside convective
!                    system
!                    index 1 nearest ground
!
!   intent(in), optional variables:
!
!     kbot           index of lowest model level (only needed for eta
!                    coordinate case)

!     These variables are present when a prognostic cloud scheme is 
!     employed:

!     cf             large-scale cloud fraction  
!                    [ nondimensional ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(liquid) / kg ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(ice) / kg ]
!
!   intent(out), optional variables:
!
!     These variables are present when a prognostic cloud scheme is 
!     employed:
!
!     delta_qa       change in cloud area due to deep convection
!                    during the time step [ dimensionless ]
!     delta_ql       change in cloud water due to deep convection 
!                    during the time step [ (kg/kg) ]
!     delta_qi       change in cloud ice due to deep convection 
!                    during the time step [ (kg/kg) ]
!     mtot           mass flux at model levels, convective plus meso-
!                    scale, due to donner_deep_mod [ (kg/m**2) / sec ]
!
!
!---------------------------------------------------------------------

      call error_mesg('donner_deep', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep



!####################################################################

subroutine donner_deep_avg (is, ie, js, je,   &
                            cell_cloud_frac_out,  cell_liquid_amt_out,&
                            cell_liquid_size_out, cell_ice_amt_out,   &
                            cell_ice_size_out,    meso_cloud_frac_out,&
                            meso_liquid_amt_out,  meso_liquid_size_out,&
                            meso_ice_amt_out,     meso_ice_size_out)

!------------------------------------------------------------------
!   this subroutine provides the cloud microphysical quantities assoc-
!   iated with donner_deep convection to the radiation package. these
!   include the cloud liquid and ice amounts, liquid and ice sizes,
!   and fractional coverage, for both the convective cell and mesoscale
!   components of the convective system.
!------------------------------------------------------------------

integer, intent(in)                    :: is, ie, js, je
real,    intent(out), dimension(:,:,:) :: cell_cloud_frac_out,        &
                                          cell_liquid_amt_out,   &
                                          cell_liquid_size_out,  &
                                          cell_ice_amt_out, &
                                          cell_ice_size_out,   &
                                          meso_cloud_frac_out,   &
                                          meso_liquid_amt_out, &
                                          meso_liquid_size_out, &
                                          meso_ice_amt_out,  &
                                          meso_ice_size_out

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!
!  intent(out) variables:
!
!     cell_cloud_frac_out  fractional coverage of convective cells in 
!                          grid box [ dimensionless ]
!     cell_liquid_amt_out  liquid water content of convective cells 
!                          [ kg(h2o)/kg(air) ]
!     cell_liquid_size_out assumed effective size of cell liquid drops
!                          [ micrometers ]
!     cell_ice_amt_out     ice water content of cells 
!                          [ kg(h2o)/kg(air) ]
!     cell_ice_size_out    generalized effective diameter for ice in
!                          convective cells [ micrometers ]
!     meso_cloud_frac_out  fractional area of mesoscale clouds in grid 
!                          box [ dimensionless ]
!     meso_liquid_amt_out  liquid water content in mesoscale clouds
!                          currently set to 0.0
!                          [ kg(h2o)/kg(air) ]
!     meso_liquid_size_out assumed effective size of mesoscale drops
!                          currently set to 0.0 [ micrometers ]
!     meso_ice_amt_out     ice water content of mesoscale elements 
!                          [ kg(h2o)/kg(air) ]
!     meso_ice_size_out    generalized ice effective size for anvil ice
!                          [ micrometers ]
!
!---------------------------------------------------------------------

      call error_mesg('donner_deep_avg', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_avg




!##################################################################


subroutine donner_deep_end

!---------------------------------------------------------------------
!   this is the destructor for donner_deep_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------

      module_is_initialized = .false. 

!---------------------------------------------------------------------

      call error_mesg('donner_deep_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_end

!######################################################################

  end module donner_deep_mod

