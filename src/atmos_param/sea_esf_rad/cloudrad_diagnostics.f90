                 module cloudrad_diagnostics_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    cloudrad_diagnostics_mod generates any desired netcdf output
!    fields involving the cloud fields seen by the radiation package
!    or the cloud radiation interaction variables.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

! shared modules:

use fms_mod,                 only: fms_init, open_namelist_file, &
                                   write_version_number, mpp_pe, &
                                   mpp_root_pe, stdlog, file_exist,  &
                                   check_nml_error, error_mesg,   &
                                   FATAL, NOTE, WARNING, close_file,  &
                                   read_data, write_data
use time_manager_mod,        only: time_type, time_manager_init
use diag_manager_mod,        only: register_diag_field, send_data, &
                                   diag_manager_init

! shared radiation package modules:

use rad_utilities_mod,       only: rad_utilities_init, &
                                   cldrad_properties_type, &
                                   cld_specification_type, &
                                   solar_spectrum_type, &
                                   longwave_control_type, Lw_control, &
                                   microrad_properties_type, &
                                   microphysics_type, atmos_input_type,&
                                   Cldrad_control

use esfsw_parameters_mod,    only: Solar_spect, esfsw_parameters_init

use microphys_rad_mod,       only : isccp_microphys_sw_driver, isccp_microphys_lw_driver
!  other cloud diagnostics modules

use isccp_clouds_mod,        only: isccp_clouds_init, isccp_clouds_end,&
                                   isccp_output, isccp_cloudtypes, isccp_cloudtypes_stochastic

use constants_mod,           only: diffac

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloudrad_diagnostics_mod generates any desired netcdf output
!    fields involving the cloud fields seen by the radiation package
!    or the cloud radiation interaction variables.
!
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloudrad_diagnostics.f90,v 12.0 2005/04/14 15:44:36 fms Exp $'
character(len=128)  :: tagname =  '$Name: lima $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloudrad_diagnostics_init, cloudrad_netcdf, &
         cloudrad_diagnostics_end

private          &
!   called from cloudrad_diagnostics_init:
         diag_field_init, &
!   called from cloudrad_netcdf:
         isccp_diag, isccp_diag_stochastic, compute_isccp_clds,  &
!   called from isccp_diag:  
         cloud_optical_properties_diag


!---------------------------------------------------------------------
!-------- namelist  ---------
!
! do_isccp                 should isccp_cloudtypes processing be done?
!
! isccp_actual_radprops    should the GCMs radiative properties be 
!                          used in the isccp_cloudtypes processing?
!                          If false, then use properties diagnosed
!                          locally from cloud_optical_properties_diag.
!       
! isccp_scale_factor       This scale factor is here to remove the
!                          scaling of liquid water and ice water 
!                          paths in the cloud_rad to account for the
!                          plane-parallel homogenous cloud bias.
!
!                          NOTE THAT THIS SCALE FACTOR SHOULD BE
!                          SET IDENTICAL TO THE ONE SET IN THE
!                          NAMELIST TO CLOUD_RAD.f90
!
!                          It is put here because the diagnostics
!                          are on the clouds themselves, not the
!                          radiative fluxes.  The scale factor
!                          only exists to compute radiative transfer
!                          more accurately.    

logical :: do_isccp = .false.
logical :: isccp_actual_radprops = .true.
real :: isccp_scale_factor = 0.85

 namelist /cloudrad_diagnostics_nml /  do_isccp, isccp_actual_radprops, &
                                       isccp_scale_factor


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

real, parameter     :: taumin = 1.E-06  ! minimum value allowed for 
                                        ! optical depth 
                                        ! [ dimensionless ]

real,  parameter ::   mid_btm  = 6.8e4  ! isccp boundaries
real, parameter  ::   high_btm = 4.4e4  ! isccp boundaries
      
!----------------------------------------------------------------------
!    diagnostics variables.     
!----------------------------------------------------------------------
character(len=8)    :: mod_name = 'cloudrad'
real                :: missing_value = -999.

integer :: id_tot_cld_amt, id_cld_amt, id_em_cld_lw, id_em_cld_10u, & 
           id_abs_lsc_cld_10u, id_abs_lsc_cld_lw,  &
           id_abs_cell_cld_10u, id_abs_cell_cld_lw,  &
           id_abs_meso_cld_10u, id_abs_meso_cld_lw,  &
           id_abs_cld_10u, id_abs_cld_lw,  &
           id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt,  &
           id_lsc_cld_amt, id_cell_cld_amt, id_meso_cld_amt,  &
           id_lsc_cld_ext_uv, id_lsc_cld_ext_vis, id_lsc_cld_ext_nir, &
           id_lsc_cld_sct_uv, id_lsc_cld_sct_vis, id_lsc_cld_sct_nir, &
           id_lsc_cld_asymm_uv, id_lsc_cld_asymm_vis,    &
           id_lsc_cld_asymm_nir, &
           id_cell_cld_ext_uv, id_cell_cld_ext_vis,    &
           id_cell_cld_ext_nir, &
           id_cell_cld_sct_uv, id_cell_cld_sct_vis,    &
           id_cell_cld_sct_nir, &
           id_cell_cld_asymm_uv, id_cell_cld_asymm_vis,    &
           id_cell_cld_asymm_nir, &
           id_meso_cld_ext_uv, id_meso_cld_ext_vis,   &
           id_meso_cld_ext_nir, &
           id_meso_cld_sct_uv, id_meso_cld_sct_vis,   &
           id_meso_cld_sct_nir, &
           id_meso_cld_asymm_uv, id_meso_cld_asymm_vis,    &
           id_meso_cld_asymm_nir, &
           id_ext_cld_uv,   id_sct_cld_uv,  id_asymm_cld_uv, &
           id_ext_cld_vis,  id_sct_cld_vis, id_asymm_cld_vis, &
           id_ext_cld_nir,  id_sct_cld_nir, id_asymm_cld_nir, &
           id_alb_uv_cld, id_alb_nir_cld, id_abs_uv_cld, id_abs_nir_cld
   
! Diagnostics strat cloud microphysical properties
integer::  id_strat_area_liq, id_strat_conc_drop, id_strat_size_drop,&
           id_strat_area_ice, id_strat_conc_ice, id_strat_size_ice

! Diagnostics for stochastic clouds
integer :: id_cldfrac_ave, id_cldfrac_tot,    &
           id_ice_conc_ave, id_drop_conc_ave

integer, dimension(:), allocatable ::    &
                                   id_cldfrac_cols,                &
                                   id_ice_conc_cols, id_ice_size_cols, &
                                   id_drop_conc_cols, id_drop_size_cols

logical :: module_is_initialized =                            &
                         .false.    ! module  initialized ?


!----------------------------------------------------------------------
!----------------------------------------------------------------------



                        contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################
! <SUBROUTINE NAME="cloudrad_diagnostics_init">
!  <OVERVIEW>
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_diagnostics_init (axes, Time)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_diagnostics_init (axes, Time)

!---------------------------------------------------------------------
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!------------------------------------------------------------------

integer, dimension(4),   intent(in)    ::   axes
type(time_type),         intent(in)    ::   Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer         :: unit, io, ierr

!---------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      io       error status returned from io operation  
!      ierr     error code
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call time_manager_init
      call esfsw_parameters_init
      call diag_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=cloudrad_diagnostics_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'cloudrad_diagnostics_nml')
        enddo
10      call close_file (unit)
      endif
 
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() )    &
                       write (stdlog(), nml=cloudrad_diagnostics_nml)
 
!---------------------------------------------------------------------
!    allocate the arrays holding the band-dependent cloud fractions 
!    for use when stochastic clouds are being used.
!---------------------------------------------------------------------
      allocate (id_cldfrac_cols  (Cldrad_control%nlwcldb + Solar_spect%nbands))
      allocate (id_ice_conc_cols (Cldrad_control%nlwcldb + Solar_spect%nbands))
      allocate (id_drop_conc_cols(Cldrad_control%nlwcldb + Solar_spect%nbands))
      allocate (id_ice_size_cols (Cldrad_control%nlwcldb + Solar_spect%nbands))
      allocate (id_drop_size_cols(Cldrad_control%nlwcldb + Solar_spect%nbands))
 
!-------------------------------------------------------------------
!    initialize the netcdf diagnostics provided with this module.
!-------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then
        call diag_field_init (Time, axes)
      endif

!---------------------------------------------------------------------
!    initialize isccp_clouds_init 
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then
          if (do_isccp) call isccp_clouds_init (axes, Time)
      endif 

!--------------------------------------------------------------------
!    mark the module initialized.
!--------------------------------------------------------------------
      module_is_initialized= .true.

!--------------------------------------------------------------------



end subroutine cloudrad_diagnostics_init


!###################################################################
! <SUBROUTINE NAME="cloudrad_netcdf">
!  <OVERVIEW>
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_netcdf (is, js, Time_diag, Atmos_input, cosz, &
!                            Lsc_microphys, Meso_microphys, &
!                            Cell_microphys, Lscrad_props,   &
!                            Mesorad_props, Cellrad_props, Cldrad_props,&
!                            Cld_spec, mask)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="cosz" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </IN>
!  <IN NAME="Lscrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the large-scale 
!                      clouds   
!  </IN>
!  <IN NAME="Mesorad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the meso-scale
!                      clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cellrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the convective cell
!                      clouds associated with donner convection 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties on model grid
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_netcdf (is, js, Time_diag, Atmos_input, cosz, &
                            Lsc_microphys, Meso_microphys, &
                            Cell_microphys, Lscrad_props,   &
                            Mesorad_props, Cellrad_props, Cldrad_props,&
                            Cld_spec, mask)

!---------------------------------------------------------------------
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!---------------------------------------------------------------------

integer,                        intent(in)           :: is, js
type(time_type),                intent(in)           :: Time_diag
type(atmos_input_type),         intent(in)           :: Atmos_input
real, dimension(:,:),           intent(in)           :: cosz        
type(microphysics_type),        intent(in)           :: Lsc_microphys, &
                                                        Meso_microphys,&
                                                        Cell_microphys
type(microrad_properties_type), intent(in)           :: Lscrad_props, &
                                                        Mesorad_props, &
                                                        Cellrad_props
type(cldrad_properties_type),   intent(in)           :: Cldrad_props
type(cld_specification_type),   intent(in)           :: Cld_spec       
real, dimension(:,:,:),         intent(in), optional :: mask

!-------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Time_diag       time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ]
!      cosz            cosine of zenith angle --  mean value over
!                      appropriate averaging interval
!                      [ non-dimensional ]
!      Lsc_microphys   microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!      Meso_microphys  microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!      Cell_microphys  microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!      Lscrad_props    cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!      Mesorad_props   cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!      Cellrad_props   cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!      Cldrad_props    total-cloud radiative properties,
!                      [ cldrad_properties_type ]
!      Cld_spec        variables on the model grid which define the 
!                      cloud location and amount     
!                      [ cld_specification_type ]
!
!   intent(in), optional variables:
!
!      mask              present when running eta vertical coordinate,
!                        mask to remove points below ground
!
!------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2),                       &
                      size(Atmos_input%rh2o,3))     :: cloud, tmpmask

      logical, dimension(size(Atmos_input%rh2o,1),                    &
                         size(Atmos_input%rh2o,2),                  &
                         size(Atmos_input%rh2o,3))   :: hi_cld_lvl,   &
                                                        md_cld_lvl, &
                                                        lo_cld_lvl

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2))     :: tca     

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2), 3)  :: hml_ca        

      real,    dimension(:, :, :, :), allocatable :: Tau, LwEm
      
      logical    :: used
      integer    :: kx, kCol
      integer    :: i, j, k, n, iband, isccpSwBand, isccpLwBand, nswbands
      integer    :: iuv, ivis, inir
      logical, dimension(size(Atmos_input%rh2o,1),                   &
                         size(Atmos_input%rh2o,2))  :: flagij
      integer, dimension(size(Atmos_input%rh2o,1),                   &
                         size(Atmos_input%rh2o,2))  :: count5, count6, &
                                                       count7        
      logical,    &
        dimension (Cldrad_control%nlwcldb + Solar_spect%nbands) ::   &
                                      looking_for_hi, looking_for_md, &
                                      looking_for_lo
      
!---------------------------------------------------------------------
!  local variables:
!
!      cloud                array used to hold the various netcdf 
!                           output fields as they are sent off to 
!                           diag_manager_mod
!      tca                  total column cloud amount [ dimensionless ]
!      hml_ca               total column cloud amount for isccp high, 
!                           middle and low clouds, individually
!      used                 flag returned from send_data indicating
!                           whether diag_manager_mod has received 
!                           data that was sent
!      kx                   number of model layers, number of stochastic columns
!      i,j,k                do-loop indices
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kx =  size(Cld_spec%camtsw,3)

!--------------------------------------------------------------------
!    define the number of shortwave bands and set integer correspond-
!    ance for diagnostics output
!
!    The understanding used in this code is that there are 2 
!    resolutions to the shortwave spectrum.  A high resolution with
!    25 bands and a low resolution with 18 bands.  The low resolution
!    is used conventional for AM2.      Here are the bands in the 
!    high and low res used for the UV, VIS, and NIR prescriptions
!    below.
!
!
!    For Low Resolution (nswbands = 18) :
!
!    Region   iband     Wavenumbers (cm-1)         Wavelength (microns)
!    ------   -----     ------------------         --------------------
!
!     UV       15          35300-36500                    0.274-0.283
!     VIS       7          16700-20000                      0.5-0.6
!     NIR       3           4200-8200                      1.22-2.38
!
!
!    For High Resolution (nswbands = 25) :
!
!    Region   iband     Wavenumbers (cm-1)         Wavelength (microns)
!    ------   -----     ------------------         --------------------
!
!     UV       22          35300-36500                    0.274-0.283
!     VIS      12          16700-20000                      0.5-0.6
!     NIR       8           6200-8200                      1.22-1.61
!
!---------------------------------------------------------------------

      nswbands = size(Lscrad_props%cldext,4)
      if (nswbands .eq. 25) then      
          iuv=22
          ivis=12
          inir=8
      else if (nswbands .eq. 18) then
          iuv=15
          ivis=7
          inir=3
      else
          iuv=15
          ivis=7
          inir=3
      end if        
      
!---------------------------------------------------------------------
!
!
!
!                   ISCCP SIMULATOR SECTION
!
!
!
!
!---------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    when running strat clouds, call isccp_diag to generate isccp-
!    relevant diagnostics.
!
!    Note that cloud optical thickness in the visible band is sent
!    to isccp diag.  Band 6 corresponds to 14600-16700 cm-1 or 
!    0.6-0.685 microns, from 18 band structure.
!
!    If the multi-band lw emissivity is active, longwave emissivity 
!    is taken from the band closest to 10 microns (900-990 cm-1 band, 
!    10.1-11.1 microns, band 4 of 8). If the multi-band lw cloud 
!    emissivity formulation is not active, longwave emissivity is 
!    taken from band 1 (0-2200 cm-1).
!---------------------------------------------------------------------

      if (do_isccp) then
          if (Cldrad_control%do_strat_clouds) then
!
! Which bands to use for ISCCP cloud detection?
!
            select case(nswbands)
              case (25) 
                isccpSwBand = 11
              case (18) 
                isccpSwBand = 6
              case default
                isccpSwBand = 6
            end select
            if (Cldrad_control%do_lw_micro) then
              isccpLwBand = 4
            else
              isccpLwBand = 1    
            end if
      
            if (Cldrad_control%do_stochastic_clouds) then
            ! Cloud fraction, averaged across bands
            !
              allocate( Tau(size(Lsc_microphys%stoch_cldamt, 1),  &
                            size(Lsc_microphys%stoch_cldamt, 2),  &
                            size(Lsc_microphys%stoch_cldamt, 3),  &
                            size(Lsc_microphys%stoch_cldamt, 4)), &
                       LwEm(size(Lsc_microphys%stoch_cldamt, 1),  &
                            size(Lsc_microphys%stoch_cldamt, 2),  &
                            size(Lsc_microphys%stoch_cldamt, 3),  &
                            size(Lsc_microphys%stoch_cldamt, 4)) )
 
              kCol = size(Lsc_microphys%stoch_cldamt, 4)
      
! After this call the Tau array is actually extinction
!
              call isccp_microphys_sw_driver (is, js, isccpSwBand, &
                                           Lsc_microphys, cldext = Tau) 
! And to get optical thickness...
! 
              Tau(:,:,:,:) = (Tau(:, :, :, :) *                  &
                     spread(Atmos_input%deltaz(:,:,:)/1000.,   &
                      dim = 4, nCopies = size(Tau, 4)) / &
                     isccp_scale_factor)
     !
     ! At first the LwEm array holds the absorption coefficient...
     !
      call isccp_microphys_lw_driver (is, js, isccpLwBand, &
                                              Lsc_microphys,abscoeff = LwEm)
              !
      ! ... and then the emissivity 
      !
              LwEm(:, :, :, :) = 1. - exp( -1. * diffac *            &
                                          (LwEm(:, :, :, :) *        &
               spread(Atmos_input%deltaz(:,:,:)/1000.,    &
                 dim = 4, nCopies = size(Tau, 4))) / &
                isccp_scale_factor)  
                        
              call isccp_diag_stochastic (is, js, Atmos_input, cosz, &
                                          Tau, LwEm,   &
                                          Lsc_microphys%stoch_cldamt, &
                                          Time_diag)
            else
              allocate( Tau(size(Atmos_input%rh2o,1),     &
                            size(Atmos_input%rh2o,2),     &
                            size(Atmos_input%rh2o,3), 1), &
                       LwEm(size(Atmos_input%rh2o,1),     &
                            size(Atmos_input%rh2o,2),     &
                            size(Atmos_input%rh2o,3), 1))
              Tau(:,:,:, 1) = (Lscrad_props%cldext(:,:,:,isccpSwBand)* &
                                Atmos_input%deltaz(:,:,:)/1000.) / &
                                isccp_scale_factor
              LwEm(:,:,:, 1) =  1. - exp( -1. * diffac *        &
                                (Lscrad_props%abscoeff(:,:,:,isccpLwBand)* &
                                Atmos_input%deltaz(:,:,:)/1000.)/ &
                                isccp_scale_factor)           
              call isccp_diag (is, js, Cld_spec, Atmos_input, cosz, &
                               Tau(:, :, :, 1), LwEm(:, :, :, 1), &
                               Time_diag)
            end if
            deallocate(Tau, LwEm)
         endif
      endif

!---------------------------------------------------------------------
!
!
!
!          COMPUTE ACTUAL HIGH, MIDDLE, AND LOW CLOUD AMOUNTS
!
!
!
!
!---------------------------------------------------------------------

if (Cldrad_control%do_stochastic_clouds) then
      kCol = size(Lsc_microphys%stoch_cldamt, 4)

!---------------------------------------------------------------------
!    if desired as a diagnostic, define the total cloud amount. send to
!    diag_manager_mod for netcdf output.
!---------------------------------------------------------------------
 
      if ( id_tot_cld_amt > 0  .or. id_cldfrac_tot > 0) then
        count5 = 0
        do n=1, kCol
          do j=1, size(tca,2)
            do i=1, size(tca,1)
              do k=1, size(Lsc_microphys%stoch_cldamt, 3)
                if(Lsc_microphys%stoch_cldamt(i,j,k,n) > 0.0) then
                  count5(i,j) = count5(i,j) + 1
                  exit
                endif
              end do
            end do
          end do
        end do

        do j=1, size(tca,2)
          do i=1, size(tca,1)
            tca (i,j) = 100.0*(real(count5(i,j))/real(kCol))
          end do
        end do

        used = send_data (id_tot_cld_amt, tca, Time_diag, is, js)
      end if 

!---------------------------------------------------------------------
!    if high, mid or low cloud diagnostics are desired, call 
!    compute_isccp_clds to define the amount of each. send to 
!    diag_manager_mod.
!---------------------------------------------------------------------
      if (id_high_cld_amt > 0 .or. id_mid_cld_amt > 0 .or. &
          id_low_cld_amt > 0) then
                         
        do k = 1,size(Atmos_input%pflux,3)-1
          do j=1,size(Atmos_input%pflux,2)
            do i=1,size(Atmos_input%pflux,1)
              hi_cld_lvl(i,j,k) = Atmos_input%pflux(i,j,k)  <=  high_btm
              md_cld_lvl(i,j,k) = (Atmos_input%pflux(i,j,k) <= mid_btm &
                                       .and. &
                                   Atmos_input%pflux(i,j,k) >  high_btm)
              lo_cld_lvl(i,j,k) = Atmos_input%pflux(i,j,k ) >  mid_btm
            end do
          end do
        end do

        count5 = 0
        count6 = 0
        count7 = 0
        do j=1,size(Atmos_input%pflux,2)
          do i=1,size(Atmos_input%pflux,1)
            looking_for_hi = .true.
            looking_for_md = .true.
            looking_for_lo = .true.
            do k = 1,size(Atmos_input%pflux,3)-1
              if (hi_cld_lvl(i,j,k)) then
                do n=1,kCol
                  if (looking_for_hi(n)) then
                    if (Lsc_microphys%stoch_cldamt(i,j,k,n) > 0. ) then
                      count5(i,j) = count5(i,j) + 1
                      looking_for_hi(n) =  .false.
                    endif
                  endif
                end do
              else if (md_cld_lvl(i,j,k)) then
                do n=1,kCol
                  if (looking_for_md(n)) then
                    if (Lsc_microphys%stoch_cldamt(i,j,k,n) > 0. ) then
                      count6(i,j) = count6(i,j) + 1
                      looking_for_md(n) = .false.
                    endif
                  endif
                end do
              else if (lo_cld_lvl(i,j,k)) then
                do n=1,kCol
                  if (looking_for_lo(n)) then
                    if (Lsc_microphys%stoch_cldamt(i,j,k,n) > 0. ) then
                      count7(i,j) = count7(i,j) + 1
                      looking_for_lo(n) = .false.
                    endif
                  endif
                end do
              endif
            end do
            hml_ca(i,j,1) = real(count5(i,j))/real(kCol)
            hml_ca(i,j,2) = real(count6(i,j))/real(kCol)
            hml_ca(i,j,3) = real(count7(i,j))/real(kCol)
          end do
        end do

        hml_ca = 100.*hml_ca

        if (id_high_cld_amt > 0)  used =    &
           send_data (id_high_cld_amt, hml_ca(:,:,1), Time_diag, is, js)
        if (id_mid_cld_amt > 0)  used =     &
           send_data (id_mid_cld_amt,  hml_ca(:,:,2), Time_diag, is, js)
        if (id_low_cld_amt > 0)  used =    &
           send_data (id_low_cld_amt,  hml_ca(:,:,3), Time_diag, is, js)
      
      endif   !do high, middle or low cloud amount

else       !branched on do_stochastic clouds

!---------------------------------------------------------------------
!    if desired as a diagnostic, define the total cloud amount. send to
!    diag_manager_mod for netcdf output.
!---------------------------------------------------------------------
      if ( id_tot_cld_amt > 0 ) then
        tca = 1.0
        do k=1,kx        
          tca(:,:) = tca(:,:)*(1.0 - Cld_spec%camtsw(:,:,k))
        end do
        tca = 100.*(1. - tca)
        used = send_data (id_tot_cld_amt, tca, Time_diag, is, js)
      endif

!---------------------------------------------------------------------
!    if high, mid or low cloud diagnostics are desired, call 
!    compute_isccp_clds to define the amount of each. send to 
!    diag_manager_mod.
!---------------------------------------------------------------------
      if (id_high_cld_amt > 0 .or. id_mid_cld_amt > 0 .or. &
          id_low_cld_amt > 0) then
        call compute_isccp_clds (Atmos_input%pflux, Cld_spec%camtsw, &
                                 Cld_spec%camtsw_band, hml_ca)
   
        if (id_high_cld_amt > 0)  used =    &
           send_data (id_high_cld_amt, hml_ca(:,:,1), Time_diag, is, js)
        if (id_mid_cld_amt > 0)  used =     &
           send_data (id_mid_cld_amt, hml_ca(:,:,2), Time_diag, is, js)
        if (id_low_cld_amt > 0)  used =    &
           send_data (id_low_cld_amt, hml_ca(:,:,3), Time_diag, is, js)
      endif

end if

!---------------------------------------------------------------------
!
!
!
!                   3 DIMENSION CLOUD AMOUNT
!
!
!
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    send the 3D cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
      if (id_cld_amt > 0)  used =        &
        send_data (id_cld_amt, Cld_spec%camtsw, Time_diag, is, js, 1, &
                   rmask=mask)

!---------------------------------------------------------------------
!
!
!
!              SHORTWAVE RADIATIVE PROPERTIES OF STRATIFORM CLOUDS
!
!
!
!
!---------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    the following diagnostics are meaningful only when strat_clouds
!    is active:
!----------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then

!----------------------------------------------------------------------
!    send the 3D large-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
 
        if (id_lsc_cld_amt > 0)  used =      &
           send_data (id_lsc_cld_amt, Lsc_microphys%cldamt, Time_diag, &
                      is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    send various large-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_lsc_cld_ext_uv > 0) then
          cloud(:,:,:) = Lscrad_props%cldext(:,:,:,iuv)
          used = send_data (id_lsc_cld_ext_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_ext_vis > 0) then
          cloud(:,:,:) = Lscrad_props%cldext(:,:,:,ivis)
          used = send_data (id_lsc_cld_ext_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_ext_nir > 0) then
          cloud(:,:,:) = Lscrad_props%cldext(:,:,:,inir)
          used = send_data (id_lsc_cld_ext_nir, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_sct_uv > 0) then
          cloud(:,:,:) = Lscrad_props%cldsct(:,:,:,iuv)
          used = send_data (id_lsc_cld_sct_uv, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_sct_vis > 0) then
          cloud(:,:,:) = Lscrad_props%cldsct(:,:,:,ivis)
          used = send_data (id_lsc_cld_sct_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_sct_nir > 0) then
          cloud(:,:,:) = Lscrad_props%cldsct(:,:,:,inir)
          used = send_data (id_lsc_cld_sct_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_asymm_uv > 0) then
          cloud(:,:,:) = Lscrad_props%cldasymm(:,:,:,iuv)
          used = send_data (id_lsc_cld_asymm_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_asymm_vis > 0) then
          cloud(:,:,:) = Lscrad_props%cldasymm(:,:,:,ivis)
          used = send_data (id_lsc_cld_asymm_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_asymm_nir > 0) then
          cloud(:,:,:) = Lscrad_props%cldasymm(:,:,:,inir)
          used = send_data (id_lsc_cld_asymm_nir, cloud, Time_diag,  &
                            is, js, 1, rmask=mask)
        endif
      endif ! (do_strat_clouds)
 
!---------------------------------------------------------------------
!
!
!
!             SHORTWAVE RADIATIVE PROPERTIES OF DONNER CLOUDS
!
!
!
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    the following diagnostics are meaningful only when 
!    donner_deep_clouds is active:
!----------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then

!----------------------------------------------------------------------
!    send the 3D cell-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_cell_cld_amt > 0 ) then
          used = send_data (id_cell_cld_amt, Cell_microphys%cldamt, &
                            Time_diag, is, js, 1, rmask=mask)
        endif

!----------------------------------------------------------------------
!    send various cell-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_cell_cld_ext_uv > 0) then
          cloud(:,:,:) = Cellrad_props%cldext(:,:,:,iuv)
          used = send_data (id_cell_cld_ext_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_ext_vis > 0) then
          cloud(:,:,:) = Cellrad_props%cldext(:,:,:,ivis)
          used = send_data (id_cell_cld_ext_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_ext_nir > 0) then
          cloud(:,:,:) = Cellrad_props%cldext(:,:,:,inir)
          used = send_data (id_cell_cld_ext_nir, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_sct_uv > 0) then
          cloud(:,:,:) = Cellrad_props%cldsct(:,:,:,iuv)
          used = send_data (id_cell_cld_sct_uv, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_sct_vis > 0) then
          cloud(:,:,:) = Cellrad_props%cldsct(:,:,:,ivis)
          used = send_data (id_cell_cld_sct_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_sct_nir > 0) then
          cloud(:,:,:) = Cellrad_props%cldsct(:,:,:,inir)
          used = send_data (id_cell_cld_sct_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_asymm_uv > 0) then
          cloud(:,:,:) = Cellrad_props%cldasymm(:,:,:,iuv)
          used = send_data (id_cell_cld_asymm_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if ( id_cell_cld_asymm_vis > 0 ) then
          cloud(:,:,:) = Cellrad_props%cldasymm(:,:,:,ivis)
          used = send_data (id_cell_cld_asymm_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if ( id_cell_cld_asymm_nir > 0 ) then
          cloud(:,:,:) = Cellrad_props%cldasymm(:,:,:,inir)
          used = send_data (id_cell_cld_asymm_nir, cloud, Time_diag,  &
                            is, js, 1, rmask=mask )
        endif

!----------------------------------------------------------------------
!    send the 3D meso-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_meso_cld_amt > 0) then
          used = send_data (id_meso_cld_amt, Meso_microphys%cldamt,  &
                            Time_diag, is, js, 1, rmask=mask)
        endif

!----------------------------------------------------------------------
!    send various meso-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_meso_cld_ext_uv > 0) then
          cloud(:,:,:) = Mesorad_props%cldext(:,:,:,iuv)
          used = send_data (id_meso_cld_ext_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_ext_vis > 0) then
          cloud(:,:,:) = Mesorad_props%cldext(:,:,:,ivis)
          used = send_data (id_meso_cld_ext_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_ext_nir > 0) then
          cloud(:,:,:) = Mesorad_props%cldext(:,:,:,inir)
          used = send_data (id_meso_cld_ext_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_sct_uv > 0) then
          cloud(:,:,:) = Mesorad_props%cldsct(:,:,:,iuv)
          used = send_data (id_meso_cld_sct_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_sct_vis > 0) then
          cloud(:,:,:) = Mesorad_props%cldsct(:,:,:,ivis)
          used = send_data (id_meso_cld_sct_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_sct_nir > 0) then
          cloud(:,:,:) = Mesorad_props%cldsct(:,:,:,inir)
          used = send_data (id_meso_cld_sct_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_asymm_uv > 0) then
          cloud(:,:,:) = Mesorad_props%cldasymm(:,:,:,iuv)
          used = send_data (id_meso_cld_asymm_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_asymm_vis > 0) then
          cloud(:,:,:) = Mesorad_props%cldasymm(:,:,:,ivis)
          used = send_data (id_meso_cld_asymm_vis, cloud, Time_diag,  &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_asymm_nir > 0) then
          cloud(:,:,:) = Mesorad_props%cldasymm(:,:,:,inir)
          used = send_data (id_meso_cld_asymm_nir, cloud, Time_diag,  &
                            is, js, 1, rmask=mask)
        endif
      endif ! (do_donner_deep_clouds)

!---------------------------------------------------------------------
!
!
!             STRATIFORM PHYSICAL PROPERTIES
!
!
!
!
!---------------------------------------------------------------------
   if (Cldrad_control%do_strat_clouds) then

      if (id_strat_area_ice > 0) then
           where(Lsc_microphys%conc_ice > 0)
                 tmpmask = Lsc_microphys%cldamt
           elsewhere
                 tmpmask = 0.
           endwhere      
           used = send_data (id_strat_area_ice, tmpmask, Time_diag,  &
                            is, js, 1, rmask=mask)
      end if
      if (id_strat_size_ice > 0) then
           where(Lsc_microphys%conc_ice > 0)
                 tmpmask = Lsc_microphys%cldamt*Lsc_microphys%size_ice
           elsewhere
                 tmpmask = 0.
           endwhere      
           used = send_data (id_strat_size_ice, tmpmask, Time_diag,  &
                            is, js, 1, rmask=mask)
      end if
      if (id_strat_conc_ice > 0) then
           used = send_data (id_strat_conc_ice, Lsc_microphys%conc_ice, &
                             Time_diag, is, js, 1, rmask=mask)
      end if
      
      if (id_strat_area_liq > 0) then
           where(Lsc_microphys%conc_drop > 0)
                 tmpmask = Lsc_microphys%cldamt
           elsewhere
                 tmpmask = 0.
           endwhere      
           used = send_data (id_strat_area_liq, tmpmask, Time_diag,  &
                            is, js, 1, rmask=mask)
      end if
      if (id_strat_size_drop > 0) then
           where(Lsc_microphys%conc_drop > 0)
                 tmpmask = Lsc_microphys%cldamt*Lsc_microphys%size_drop
           elsewhere
                 tmpmask = 0.
           endwhere      
           used = send_data (id_strat_size_drop, tmpmask, Time_diag,  &
                            is, js, 1, rmask=mask)
      end if
      if (id_strat_conc_drop > 0) then
           used = send_data (id_strat_conc_drop, Lsc_microphys%conc_drop, &
                             Time_diag, is, js, 1, rmask=mask)
      end if
                    

   end if

!---------------------------------------------------------------------
!
!
!
!             LONGWAVE RADIATIVE PROPERTIES OF CLOUDS
!
!
!
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a total cloud field emissivity that is the weighted average
!    of the random and max overlap emissivities, over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_em_cld_10u > 0) then
          cloud(:,:,:) =    &
               (Cld_spec%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,5,1) + &
                Cld_spec%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,5,1))/ &
               (Cld_spec%crndlw(:,:,:) + Cld_spec%cmxolw(:,:,:) +      &
                                                               1.0E-10)
          used = send_data (id_em_cld_10u, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define a total cloud field emissivity that is the weighted average
!    of the random and max overlap emissivities, over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        if (id_em_cld_lw > 0   ) then
          cloud(:,:,:) =      &
              (Cld_spec%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,1,1) +  &
               Cld_spec%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,1,1))/ &
              (Cld_spec%crndlw(:,:,:) + Cld_spec%cmxolw(:,:,:) +     &
                                                                1.0E-10)
          used = send_data (id_em_cld_lw, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a large scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_lsc_cld_10u > 0) then
          used = send_data (id_abs_lsc_cld_10u,    &
                            Lscrad_props%abscoeff(:,:,:,5), Time_diag, &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the large scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        if (id_abs_lsc_cld_lw > 0) then
          used = send_data (id_abs_lsc_cld_lw,      &
                            Lscrad_props%abscoeff(:,:,:,1), Time_diag, &
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    the cell scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_cell_cld_10u > 0) then
          used = send_data (id_abs_cell_cld_10u,     &
                            Cellrad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      else

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the cell scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
        if (id_abs_cell_cld_lw > 0) then
          used = send_data (id_abs_cell_cld_lw,     &
                            Cellrad_props%abscoeff(:,:,:,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a meso-scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_meso_cld_10u > 0) then
          used = send_data (id_abs_meso_cld_10u,    &
                            Mesorad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      else
 
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the meso-scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
        if ( id_abs_meso_cld_lw > 0  ) then
          used = send_data (id_abs_meso_cld_lw,    &
                            Mesorad_props%abscoeff(:,:,:,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    the total-cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_cld_10u > 0) then
          used = send_data (id_abs_cld_10u,      &
                            Cldrad_props%abscoeff(:,:,:,5,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the total-cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        if (id_abs_cld_lw > 0) then
          used = send_data (id_abs_cld_lw,    &
                            Cldrad_props%abscoeff(:,:,:,1,1), Time_diag, &
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!
!
!
!             SHORTWAVE RADIATIVE PROPERTIES OF ALL CLOUDS COMBINED
!
!
!
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define total-cloud diagnostics that are associated with the micro-
!    physically-based cloud shortwave radiative properties.
!---------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then
        if (id_ext_cld_uv > 0) then
          cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,iuv,1)
          used = send_data (id_ext_cld_uv, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_sct_cld_uv > 0) then
          cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,iuv,1)
          used = send_data (id_sct_cld_uv, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_asymm_cld_uv > 0) then
          cloud(:,:,:) = 100.0*Cldrad_props%cldasymm(:,:,:,iuv,1)
          used = send_data (id_asymm_cld_uv, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_ext_cld_vis > 0) then
          cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,ivis,1)
          used = send_data (id_ext_cld_vis, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_sct_cld_vis > 0) then
          cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,ivis,1)
          used = send_data (id_sct_cld_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_asymm_cld_vis > 0) then
          cloud(:,:,:) = 100.0*Cldrad_props%cldasymm(:,:,:,ivis,1)
          used = send_data (id_asymm_cld_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_ext_cld_nir > 0) then
          cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,inir,1)
          used = send_data (id_ext_cld_nir, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_sct_cld_nir > 0) then
          cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,inir,1)
          used = send_data (id_sct_cld_nir, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_asymm_cld_nir > 0) then
          cloud(:,:,:) = 100.0*Cldrad_props%cldasymm(:,:,:,inir,1)
          used = send_data (id_asymm_cld_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    define total-cloud diagnostics that are associated with the bulk 
!    cloud shortwave radiative properties.
!---------------------------------------------------------------------
      else

!---------------------------------------------------------------------
!    define the reflected ultra-violet.
!---------------------------------------------------------------------
        if (id_alb_uv_cld > 0) then
          cloud(:,:,:) = Cldrad_props%cvisrfsw(:,:,:)
          used = send_data (id_alb_uv_cld, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    define the reflected infra-red. 
!---------------------------------------------------------------------
        if (id_alb_nir_cld > 0) then
          cloud(:,:,:) =  Cldrad_props%cirrfsw(:,:,:)
          used = send_data (id_alb_nir_cld, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    define the absorbed  ultra-violet (not implemented).
!---------------------------------------------------------------------
!       if ( id_abs_uv_cld > 0 ) then
!         cloud = 0.0
!         used = send_data (id_abs_uv_cld, cloud, Time_diag,    &
!                           is, js, 1, rmask=mask)
!       endif

!---------------------------------------------------------------------
!    define the absorbed  infra-red.
!---------------------------------------------------------------------
        if (id_abs_nir_cld > 0) then
          cloud(:,:,:) =  Cldrad_props%cirabsw(:,:,:)
          used = send_data (id_abs_nir_cld, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
      endif 

!---------------------------------------------------------------------
!
!
!
!             STOCHASTIC CLOUD PROPERTIES
!
!
!
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    output the cloud properties when the stochastic cloud option is active.
!------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
!
! Cloud fraction, averaged across bands
!
        if (id_cldfrac_ave > 0) then
          cloud(:,:,:) =   &
                  sum(Lsc_microphys%stoch_cldamt(:,:,:,:), dim = 4) / &
                  size(Lsc_microphys%stoch_cldamt(:,:,:,:), 4)
          used = send_data (id_cldfrac_ave, cloud, Time_diag, &
                            is, js, 1, rmask=mask)
        endif
  
! Total projected cloud fraction
! RSH : NOTE THAT THIS WAS PREVIOUSLY CALCULATED AS id_tot_cld_amt (tca)

        if (id_cldfrac_tot > 0) then
          used = send_data (id_cldfrac_tot, 0.01*tca, Time_diag, is, js)
        endif

! Average water and ice contents
!
        if (id_ice_conc_ave > 0) then
          cloud(:,:,:) =    &
                sum(Lsc_microphys%stoch_conc_ice(:,:,:,:), dim = 4) / &
                size(Lsc_microphys%stoch_conc_ice(:,:,:,:), 4)
          used = send_data (id_ice_conc_ave, cloud, Time_diag, &
                            is, js, 1, rmask=mask)
        endif
        if (id_drop_conc_ave > 0) then
          cloud(:,:,:) =   &
               sum(Lsc_microphys%stoch_conc_drop(:,:,:,:), dim = 4) / &
               size(Lsc_microphys%stoch_conc_drop(:,:,:,:), 4)
          used = send_data (id_drop_conc_ave, cloud, Time_diag, &
                            is, js, 1, rmask=mask)
        endif


! Cloud fraction, ice and water contents, ice and drop sizes, band by band
!
        do n=1,size(Lsc_microphys%stoch_cldamt(:,:,:,:), 4)
  
          if (id_cldfrac_cols(n) > 0) then
            cloud(:,:,:) = Lsc_microphys%stoch_cldamt(:,:,:,n)
            used = send_data (id_cldfrac_cols(n), cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif
          if (id_ice_conc_cols(n) > 0) then
            cloud(:,:,:) = Lsc_microphys%stoch_conc_ice(:,:,:,n)
            used = send_data (id_ice_conc_cols(n), cloud, Time_diag, &
                              is, js, 1, rmask=mask)
          endif
          if (id_drop_conc_cols(n) > 0) then
            cloud(:,:,:) = Lsc_microphys%stoch_conc_drop(:,:,:,n)
            used = send_data (id_drop_conc_cols(n), cloud, Time_diag, &
                              is, js, 1, rmask=mask)
          endif
          if (id_ice_size_cols(n) > 0) then
            cloud(:,:,:) = Lsc_microphys%stoch_size_ice(:,:,:,n)
            used = send_data (id_ice_size_cols(n), cloud, Time_diag,&
                              is, js, 1, rmask=mask)
          endif
          if (id_drop_size_cols(n) > 0) then
            cloud(:,:,:) = Lsc_microphys%stoch_size_drop(:,:,:,n)
            used = send_data (id_drop_size_cols(n), cloud, Time_diag, &
                              is, js, 1, rmask=mask)
          endif
        end do
      endif
       
!------------------------------------------------------------------



end subroutine cloudrad_netcdf


!####################################################################
! <SUBROUTINE NAME="cloudrad_diagnostics_end">
!  <OVERVIEW>
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_diagnostics_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cloudrad_diagnostics_end

!-------------------------------------------------------------------
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!    close out the component modules.
!--------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds) then
          if (do_isccp) call isccp_clouds_end
        endif

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine cloudrad_diagnostics_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init (axes, Time)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine diag_field_init (Time, axes )

!---------------------------------------------------------------------
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       Time      initialization time for the netcdf output fields
!       axes      diagnostic variable axes
!
!---------------------------------------------------------------------
      character(len=8) :: chvers
      integer          :: n

!---------------------------------------------------------------------
!    register the total-cloud diagnostic fields in this module.
!---------------------------------------------------------------------
      id_tot_cld_amt = register_diag_field    &
                       (mod_name, 'tot_cld_amt', axes(1:2), Time, &
                        'total cloud amount', 'percent'            )

      id_high_cld_amt = register_diag_field   &
                        (mod_name, 'high_cld_amt', axes(1:2), Time, &
                         'high cloud amount', 'percent'            )

      id_mid_cld_amt =  register_diag_field     &
                        (mod_name, 'mid_cld_amt', axes(1:2), Time, &
                         'mid cloud amount', 'percent'            )
  
      id_low_cld_amt = register_diag_field    &
                       (mod_name, 'low_cld_amt', axes(1:2), Time, &
                        'low cloud amount', 'percent'            )

      id_cld_amt =  register_diag_field     &
                    (mod_name, 'cld_amt', axes(1:3), Time,      &
                     'cloud amount', 'percent',     &
                     missing_value=missing_value            )

      id_em_cld_lw =  register_diag_field    &
                      (mod_name, 'em_cld_lw', axes(1:3), Time, &
                       'lw cloud emissivity', 'percent',        &
                       missing_value=missing_value          )
 
      id_em_cld_10u = register_diag_field    &
                      (mod_name, 'em_cld_10u', axes(1:3), Time, &
                       'cloud emissivity 10 um band', 'percent',    &
                       missing_value=missing_value          )

      id_abs_cld_lw = register_diag_field    &
                      (mod_name, 'abs_lw', axes(1:3), Time, &
                       'cloud abs coeff lw', 'percent',        &
                       missing_value=missing_value          )

     id_abs_cld_10u = register_diag_field     &
                      (mod_name, 'abs_10u', axes(1:3), Time, &
                       'cloud abs coeff 10um band', 'percent',    &
                       missing_value=missing_value          )

!---------------------------------------------------------------------
!    register diagnostic fields associated with the bulk shortwave
!    parameterization.
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_sw_micro) then
        id_alb_uv_cld = register_diag_field     &
                        (mod_name, 'alb_uv_cld', axes(1:3), Time, &
                         'UV reflected by cloud', 'percent',       &
                         missing_value=missing_value              )

        id_alb_nir_cld = register_diag_field      &
                         (mod_name, 'alb_nir_cld', axes(1:3), Time, &
                          'IR reflected by cloud', 'percent',        &
                          missing_value=missing_value               )

!   --- do not output this field ---
!       id_abs_uv_cld =  register_diag_field    &
!                        (mod_name, 'abs_uv_cld', axes(1:3), Time, &
!                         'UV absorbed by cloud', 'percent',        &
!                         missing_value=missing_value              )

        id_abs_nir_cld = register_diag_field     &
                         (mod_name, 'abs_nir_cld', axes(1:3), Time, &
                          'IR absorbed by cloud', 'percent',         &
                          missing_value=missing_value               )

!---------------------------------------------------------------------
!    register the microphysically-based total-cloud diagnostic fields.
!---------------------------------------------------------------------
      else 
        id_ext_cld_uv = register_diag_field       &
                        (mod_name, 'ext_cld_uv', axes(1:3), Time, &
                         '.27um cloud extinction coeff', 'km-1',  &
                         missing_value=missing_value          )

        id_sct_cld_uv = register_diag_field     &
                        (mod_name, 'sct_cld_uv', axes(1:3), Time, &
                         '.27um cloud scattering coeff', 'km-1', &
                         missing_value=missing_value          )

        id_asymm_cld_uv = register_diag_field    &
                          (mod_name, 'asymm_cld_uv', axes(1:3), Time, &
                           '.27um cloud asymmetry parameter',   &
                           'percent', missing_value=missing_value  ) 

        id_ext_cld_vis =  register_diag_field     &
                          (mod_name, 'ext_cld_vis', axes(1:3), Time, &
                           '.55um cloud extinction coeff', 'km-1', &
                           missing_value=missing_value          )

        id_sct_cld_vis = register_diag_field    &
                         (mod_name, 'sct_cld_vis', axes(1:3), Time, &
                          '.55um cloud scattering coeff', 'km-1', &
                          missing_value=missing_value          )

        id_asymm_cld_vis = register_diag_field      &
                           (mod_name, 'asymm_cld_vis', axes(1:3), Time,&
                            '.55um cloud asymmetry parameter',   &
                            'percent', missing_value=missing_value  )
        id_ext_cld_nir = register_diag_field    &
                         (mod_name, 'ext_cld_nir', axes(1:3), Time, &
                          '1.4um cloud extinction coeff', 'km-1', &
                          missing_value=missing_value          )

        id_sct_cld_nir = register_diag_field    &
                         (mod_name, 'sct_cld_nir', axes(1:3), Time, &
                          '1.4um cloud scattering coeff', 'km-1', &
                          missing_value=missing_value          )
 
        id_asymm_cld_nir = register_diag_field   &
                           (mod_name, 'asymm_cld_nir', axes(1:3), Time,&
                            '1.4um cloud asymmetry parameter',   &
                            'percent', missing_value=missing_value )

!---------------------------------------------------------------------
!    register the microphysically-based large-scale cloud diagnostic 
!    fields.
!---------------------------------------------------------------------
        id_lsc_cld_ext_uv = register_diag_field    &
                            (mod_name, 'lsc_cld_ext_uv', axes(1:3),   &
                             Time, '.27um lsc cloud ext coeff', 'km-1',&
                             missing_value=missing_value            )

        id_lsc_cld_ext_vis = register_diag_field   &
                             (mod_name, 'lsc_cld_ext_vis', axes(1:3), &
                              Time, '.55um lsc cloud ext coeff',   &
                              'km-1', missing_value=missing_value  )
 
        id_lsc_cld_ext_nir = register_diag_field    &
                             (mod_name, 'lsc_cld_ext_nir', axes(1:3),  &
                              Time, '1.4um lsc cloud ext coeff',   &
                              'km-1', missing_value=missing_value    )

        id_lsc_cld_sct_uv = register_diag_field    &
                            (mod_name, 'lsc_cld_sct_uv', axes(1:3),  &
                             Time, '.27um lsc cloud sct coeff', 'km-1',&
                             missing_value=missing_value            )

        id_lsc_cld_sct_vis = register_diag_field    &
                             (mod_name, 'lsc_cld_sct_vis', axes(1:3),  &
                              Time, '.55um lsc cloud sct coeff',  &
                              'km-1', missing_value=missing_value  )

        id_lsc_cld_sct_nir = register_diag_field    &
                             (mod_name, 'lsc_cld_sct_nir', axes(1:3), &
                              Time, '1.4um lsc cloud sct coeff',  &
                              'km-1', missing_value=missing_value  )
 
        id_lsc_cld_asymm_uv = register_diag_field   &
                              (mod_name, 'lsc_cld_asymm_uv', axes(1:3),&
                               Time, '.27um lsc cloud asymm coeff',  &
                               'percent', missing_value=missing_value  )

        id_lsc_cld_asymm_vis = register_diag_field  &
                               (mod_name, 'lsc_cld_asymm_vis',  &
                                axes(1:3), Time,    &
                                '.55um lsc cloud asymm coeff',   &
                                'percent', missing_value=missing_value)

        id_lsc_cld_asymm_nir = register_diag_field   &
                               (mod_name, 'lsc_cld_asymm_nir',  &
                                axes(1:3), Time,   &
                                '1.4um lsc cloud asymm coeff', &
                                'percent', missing_value=missing_value)

      endif

!---------------------------------------------------------------------
!    register the microphysically-based large-scale cloud amount and
!    lw properties diagnostic fields.
!---------------------------------------------------------------------
      id_lsc_cld_amt = register_diag_field    &
                       (mod_name, 'lsc_cld_amt', axes(1:3), Time, &
                        'lsc cloud amount', 'percent',             &
                        missing_value=missing_value            )

      id_abs_lsc_cld_lw = register_diag_field    &
                          (mod_name, 'lsc_abs_lw', axes(1:3), Time, &
                           'lsc cloud abs coeff lw', 'percent',     &
                           missing_value=missing_value          )

      id_abs_lsc_cld_10u = register_diag_field   &
                           (mod_name, 'lsc_abs_10u', axes(1:3), Time,&
                            'lsc cloud abs coeff 10um band',   &
                            'percent', missing_value=missing_value   )

!---------------------------------------------------------------------
!    register the cell-scale cloud diagnostic fields.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
        id_cell_cld_amt = register_diag_field    &
                          (mod_name, 'cell_cld_amt', axes(1:3), Time, &
                           'cell cloud amount', 'percent',             &
                           missing_value=missing_value            )

        id_cell_cld_ext_uv = register_diag_field   &
                             (mod_name, 'cell_cld_ext_uv', axes(1:3), &
                              Time, '.27um cell cloud ext coeff',  &
                              'km-1', missing_value=missing_value  )

        id_cell_cld_ext_vis = register_diag_field   &
                              (mod_name, 'cell_cld_ext_vis', axes(1:3),&
                               Time, '.55um cell cloud ext coeff',  &
                               'km-1', missing_value=missing_value  )

        id_cell_cld_ext_nir = register_diag_field   &
                              (mod_name, 'cell_cld_ext_nir', axes(1:3),&
                               Time, '1.4um cell cloud ext coeff',  &
                               'km-1', missing_value=missing_value  )

        id_cell_cld_sct_uv = register_diag_field    &
                             (mod_name, 'cell_cld_sct_uv', axes(1:3),  &
                              Time, '.27um cell cloud sct coeff',   &
                              'km-1', missing_value=missing_value    )

        id_cell_cld_sct_vis = register_diag_field    &
                              (mod_name, 'cell_cld_sct_vis', axes(1:3),&
                               Time, '.55um cell cloud sct coeff',  &
                               'km-1', missing_value=missing_value   )

        id_cell_cld_sct_nir = register_diag_field    &
                              (mod_name, 'cell_cld_sct_nir', axes(1:3),&
                               Time, '1.4um cell cloud sct coeff', &
                               'km-1', missing_value=missing_value )

        id_cell_cld_asymm_uv = register_diag_field    &
                               (mod_name, 'cell_cld_asymm_uv',   &
                                axes(1:3), Time,     &
                                '.27um cell cloud asymm coeff',   &
                                'percent', missing_value=missing_value )

        id_cell_cld_asymm_vis = register_diag_field     &
                                (mod_name, 'cell_cld_asymm_vis',   &
                                 axes(1:3), Time,     &
                                 '.55um cell cloud asymm coeff',   &
                                 'percent', missing_value=missing_value)

        id_cell_cld_asymm_nir = register_diag_field    &
                                (mod_name, 'cell_cld_sct_nir',    &
                                 axes(1:3), Time,&
                                 '1.4um cell cloud sct coeff',    &
                                 'percent', missing_value=missing_value)

        id_abs_cell_cld_lw = register_diag_field    &
                             (mod_name, 'cell_abs_lw', axes(1:3), Time,&
                              'cell cloud abs coeff lw', 'percent',   &
                              missing_value=missing_value          )

        id_abs_cell_cld_10u = register_diag_field    &
                              (mod_name, 'cell_abs_10u', axes(1:3), &
                               Time, 'cell cloud abs coeff 10um band', &
                               'percent',  missing_value=missing_value )

!---------------------------------------------------------------------
!    register the meso-scale cloud diagnostic fields.
!---------------------------------------------------------------------
        id_meso_cld_amt = register_diag_field     &
                          (mod_name, 'meso_cld_amt', axes(1:3), Time, &
                           'meso cloud amount', 'percent',             &
                           missing_value=missing_value            )

        id_meso_cld_ext_uv = register_diag_field    &
                             (mod_name, 'meso_cld_ext_uv', axes(1:3), &
                              Time, '.27um meso cloud ext coeff',   &
                              'km-1', missing_value=missing_value)

        id_meso_cld_ext_vis = register_diag_field   &
                              (mod_name, 'meso_cld_ext_vis', axes(1:3),&
                               Time, '.55um meso cloud ext coeff',  &
                               'km-1', missing_value=missing_value )

        id_meso_cld_ext_nir = register_diag_field   &
                              (mod_name, 'meso_cld_ext_nir', axes(1:3),&
                               Time, '1.4um meso cloud ext coeff',  &
                               'km-1', missing_value=missing_value    )

        id_meso_cld_sct_uv = register_diag_field   &
                             (mod_name, 'meso_cld_sct_uv', axes(1:3), &
                              Time, '.27um meso cloud sct coeff',  &
                              'km-1', missing_value=missing_value )

        id_meso_cld_sct_vis = register_diag_field  &
                              (mod_name, 'meso_cld_sct_vis', axes(1:3),&
                               Time, '.55um meso cloud sct coeff',  &
                               'km-1', missing_value=missing_value   )

        id_meso_cld_sct_nir = register_diag_field  &
                              (mod_name, 'meso_cld_sct_nir', axes(1:3),&
                               Time, '1.4um meso cloud sct coeff',  &
                               'km-1', missing_value=missing_value )

        id_meso_cld_asymm_uv = register_diag_field  &
                               (mod_name, 'meso_cld_asymm_uv',   &
                                axes(1:3), Time,     &
                                '.27um meso cloud asymm coeff',   &
                                'percent', missing_value=missing_value)

        id_meso_cld_asymm_vis = register_diag_field   &
                                (mod_name, 'meso_cld_asymm_vis',  &
                                 axes(1:3), Time,      &
                                 '.55um meso cloud asymm coeff',    &
                                 'percent', missing_value=missing_value)

        id_meso_cld_asymm_nir = register_diag_field    &
                                (mod_name, 'meso_cld_sct_nir',   &
                                 axes(1:3), Time,&
                                 '1.4um meso cloud sct coeff',   &
                                 'percent', missing_value=missing_value)

        id_abs_meso_cld_lw = register_diag_field    &
                             (mod_name, 'meso_abs_lw', axes(1:3), Time,&
                              'meso cloud abs coeff lw', 'percent',   &
                              missing_value=missing_value          )

        id_abs_meso_cld_10u = register_diag_field   &
                              (mod_name, 'meso_abs_10u', axes(1:3),  &
                               Time, 'meso cloud abs coeff 10um band', &
                               'percent', missing_value=missing_value)
      endif

!--------------------------------------------------------------------
!     register stratiform microphysical properties
!--------------------------------------------------------------------

     id_strat_area_liq = register_diag_field     &
                      (mod_name, 'strat_area_liq', axes(1:3), Time, &
                       'Area of stratiform liquid clouds', &
                       'fraction',    &
                       missing_value=missing_value          )
     id_strat_conc_drop = register_diag_field     &
                      (mod_name, 'strat_conc_drop', axes(1:3), Time, &
                       'In-cloud liquid water content of stratifrom clouds', &
                       'grams/m3',    &
                       missing_value=missing_value          )
     id_strat_size_drop = register_diag_field     &
                      (mod_name, 'strat_size_drop', axes(1:3), Time, &
                       'Effective diameter for liquid clouds', &
                       'microns',    &
                       missing_value=missing_value          )
     id_strat_area_ice = register_diag_field     &
                      (mod_name, 'strat_area_ice', axes(1:3), Time, &
                       'Area of stratiform ice clouds', &
                       'fraction',    &
                       missing_value=missing_value          )
     id_strat_conc_ice = register_diag_field     &
                      (mod_name, 'strat_conc_ice', axes(1:3), Time, &
                       'In-cloud ice water content of stratiform clouds', &
                       'grams/m3',    &
                       missing_value=missing_value          )
     id_strat_size_ice = register_diag_field     &
                      (mod_name, 'strat_size_ice', axes(1:3), Time, &
                       'Effective diameter for stratiform ice clouds', &
                       'microns',    &
                       missing_value=missing_value          )
     
!--------------------------------------------------------------------
!    register the stochastic cloud fraction arrays for each band.
!    We could do this only if Cldrad_control%do_stochastic_clouds_iz 
!    and then Cldrad_control%do_stochastic_clouds are .true., but then 
!    we'd have to be sure this module was initalized after Cldrad_control
!    was initialized. 
!--------------------------------------------------------------------
     !
     ! Cloud fraction layer-by-layer and projected
     !
     id_cldfrac_ave = register_diag_field  &
            (mod_name, 'stoch_cld_ave', axes(1:3), Time, &
             'ave cloud fraction seen by stochastic clouds', &
             'fraction', missing_value=missing_value)
     id_cldfrac_tot = register_diag_field  &
            (mod_name, 'stoch_cld_tot', axes(1:2), Time, &
             'total projected cloud fraction seen by stochastic clouds', &
             'fraction', missing_value=missing_value)
     !
     ! Mean ice and water contents
     !
     id_ice_conc_ave = register_diag_field  &
            (mod_name, 'stoch_ice_conc_ave', axes(1:3), Time, &
             'Grid-box mean ice water content seen by stochastic clouds', &
             'g/m3', missing_value=missing_value)
     id_drop_conc_ave = register_diag_field  &
            (mod_name, 'stoch_drop_conc_ave', axes(1:3), Time, &
             'Grid box mean liquid water content seen by stochastic clouds', &
             'g/m3', missing_value=missing_value)
     !
     ! Cloud properties in each sub-column
     !     
     do n=1, Cldrad_control%nlwcldb + Solar_spect%nbands
        if (n < 10) then
          write (chvers,'(i1)') n 
        else if (n <100) then
          write (chvers,'(i2)') n 
        else
          call error_mesg ('cloudrad_diagnostics_mod', &
             'must modify code to allow writing of more than 99 columns', FATAL)
        endif
!
! Cloud fraction
       ! 
        id_cldfrac_cols(n) = register_diag_field  &
            (mod_name, 'stoch_cld_col_'//trim(chvers), axes(1:3), Time, &
             'cloud fraction in stochastic column '//trim(chvers), &
             'fraction', missing_value=missing_value)
        ! 
! Ice and water contents
!
        id_ice_conc_cols(n) = register_diag_field  &
            (mod_name, 'stoch_ice_conc_col_'//trim(chvers), axes(1:3), Time, &
             'ice concentration in stochastic column '//trim(chvers), &
             'g/m3', missing_value=missing_value)
        id_drop_conc_cols(n) = register_diag_field  &
            (mod_name, 'stoch_drop_conc_col_'//trim(chvers), axes(1:3), Time, &
             'water concentration in stochastic column '//trim(chvers), &
             'g/m3', missing_value=missing_value)
! 
! Ice and water sizes
!
        id_ice_size_cols(n) = register_diag_field  &
            (mod_name, 'stoch_ice_size_col_'//trim(chvers), axes(1:3), Time, &
             'ice particle dimension in stochastic column '//trim(chvers), &
             'microns', missing_value=missing_value)
        id_drop_size_cols(n) = register_diag_field  &
            (mod_name, 'stoch_drop_size_col_'//trim(chvers), axes(1:3), Time, &
             'drop radius in stochastic column '//trim(chvers), &
             'microns', missing_value=missing_value)
     end do


!---------------------------------------------------------------------

 
end subroutine diag_field_init



!#####################################################################
! <SUBROUTINE NAME="isccp_diag">
!  <OVERVIEW>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_diag (is, js, Cld_spec, Atmos_input, coszen, Time)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <IN NAME="Lsctau" TYPE="real">
!   cloud optical thickness in the visible
!  </IN>
!  <IN NAME="Lsclwem" TYPE="real">
!   10 micron cloud emissivity
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_diag (is, js, Cld_spec, Atmos_input, coszen,       &
                       Lsctau, Lsclwem, Time)

!--------------------------------------------------------------------
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!---------------------------------------------------------------------
 
integer,                      intent(in)   :: is,js
type(cld_specification_type), intent(in)   :: Cld_spec
type(atmos_input_type),       intent(in)   :: Atmos_input
real, dimension(:,:),         intent(in)   :: coszen
real, dimension(:,:,:),       intent(in)   :: Lsctau, Lsclwem
type(time_type),              intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Cld_spec        cloud specification properties on model grid,
!                      [ cld_specification_type ]
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ] 
!      coszen          cosine of zenith angle [ dimensionless ]
!      Lsctau          0.6-0.685 micron cloud optical thickness 
!                      [ dimensionless ]
!      Lsclwem         Longwave emissivity [ dimensionless ]
!                      This is from 10.1-11.1 microns if the multiband 
!                      longwave emissivity is active.
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3)) :: &
                                    tau_local, em_local, cldamt_local

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3) ) ::  qv, em_lw_local, rh2o

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3)+1 ) ::  temp   

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3), 4 ) ::  tau

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                      7, 7) ::       fq_isccp

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2)) :: &
                          npoints, ninhomog, inhomogeneity_parameter

      integer      :: kdim
      integer      :: max_cld
      integer      :: i, j, k
      integer, dimension(size(Cld_spec%lwp,1),size(Cld_spec%lwp,2)):: &
                        sunlit
      
!---------------------------------------------------------------------
!   local variables:
!
!      tau_local        optical depth in band 1 in the current column
!                       [ dimensionless ]
!      em_local         lw cloud emissivity in current column
!                       [ dimensionless ]
!      cldamt_local     cloud fraction in current column 
!                       [ dimensionless ]
!      qv               water vapor specific humidity
!                       [ kg vapor / kg air ]
!      em_lw_local      lw cloud emissivity [ dimensionless ]
!      rh2o             mixing ratio of water vapor at model full levels
!                       [ non-dimensional ]
!      temp             temperature at model levels (1:nlev), surface
!                       temperature is stored at value nlev+1; if eta
!                       coordinates, surface value stored in below 
!                       ground points
!                       [ deg K ]
!      tau              optical depth in 4 bands [ dimensionless ]
!      fq_isccp         matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges
!      npoints          flag indicating whether isccp cloud is present
!                       in column (cloud + daylight needed)
!      kdim             number of model layers
!      max_cld          greatest number of clouds in any column in the
!                       current physics window
!      i,j,k            do-loop indices
!     
!      sunlit           is the given i,j point sunlit?
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define number of model layers.
!----------------------------------------------------------------------
      kdim = size (Cld_spec%lwp,3)
        
!---------------------------------------------------------------------
!    If optical properties are needed and no clouds exist in the 
!    window, call cloud_optical_properties_diag to define the cloud 
!    optical depth, the optical depth due to cloud ice and the 
!    longwave emissivity. If no clouds exist in the window, all the 
!    optical depths and emissivities are left are their initial
!    values of zero.to zero.
!---------------------------------------------------------------------
     
     em_lw_local = 0.
     tau = 0.

     max_cld = MAXVAL(Cld_spec%ncldsw(:,:))
     if (max_cld >= 1 .and. .not.isccp_actual_radprops)          &
       call cloud_optical_properties_diag (Cld_spec, tau, em_lw_local)

!---------------------------------------------------------------------
!    Initialize fields
!---------------------------------------------------------------------

           npoints(:,:) = 0.
           fq_isccp(:,:,:,:) = 0.
           ninhomog(:,:) = 0.
           inhomogeneity_parameter(:,:) = 0.

!---------------------------------------------------------------------
!    Compute sunlit integer flag
!---------------------------------------------------------------------
           sunlit(:,:) = 0
           where(coszen(:,:) > 1.E-06) sunlit(:,:) = 1

!--------------------------------------------------------------------
!    define the specific humidity from the mixing ratio which has been
!    input.
!--------------------------------------------------------------------
                  qv(:,:,:) = Atmos_input%cloudvapor(:,:,:)/   &
                              (1. + Atmos_input%cloudvapor(:,:,:))
                
!---------------------------------------------------------------------
!    define the column values of cloud fraction, cloud optical depth, 
!    and lw cloud emissivity. if cloud is not present, set these var-
!    iables to clear sky values.
!---------------------------------------------------------------------
           do j=1,size(Cld_spec%lwp,2)
            do i=1,size(Cld_spec%lwp,1)
             do k=1,kdim                      

                  if (Cld_spec%camtsw(i,j,k) > 0.0) then
                    cldamt_local(i,j,k) = Cld_spec%camtsw(i,j,k) 
                    if (isccp_actual_radprops) then
                      tau_local(i,j,k) = Lsctau(i,j,k)
                      em_local(i,j,k) = Lsclwem(i,j,k)
                    else 
                      tau_local(i,j,k) = tau(i,j,k,1)/ &
                                 real(Cld_spec%cld_thickness(i,j,k))
                      em_local(i,j,k) = 1.-((1.-em_lw_local(i,j,k))** &
                             (1./real(Cld_spec%cld_thickness(i,j,k))) )
                    end if          
                  else
                    cldamt_local(i,j,k) = 0.
                    tau_local(i,j,k) = 0.
                    em_local(i,j,k) = 0.
                  endif
                  
                end do
               end do
              end do
               
!---------------------------------------------------------------------
!    call isccp_cloudtypes to map each model cloud to an isccp cloud
!    type, based on its optical depth and height above the surface.
!    set a flag to indicate the presence of isccp cloud in this column.
!---------------------------------------------------------------------
                call isccp_cloudtypes (sunlit(:,:), &
                                       Atmos_input%press(:,:,1:kdim), &
                                       Atmos_input%pflux(:,:,:),&
                                       qv(:,:,:),       &
                                       Atmos_input%cloudtemp(:,:,:),  &
                                       Atmos_input%temp(:,:,kdim+1),  &
                                       cldamt_local, tau_local,   &
                                       em_local, fq_isccp(:,:,:,:),   &
                                       npoints(:,:), &
                                       inhomogeneity_parameter(:,:), &
                                       ninhomog(:,:))
 
!----------------------------------------------------------------------
!    send any desired diagnostics to the diag_manager_mod.
!----------------------------------------------------------------------
              
                
          call isccp_output (is, js, fq_isccp, npoints, &
                             inhomogeneity_parameter, ninhomog, Time)
   
!---------------------------------------------------------------------
    
    
end subroutine isccp_diag        

!#####################################################################
! <SUBROUTINE NAME="isccp_diag_stochastic">
!  <OVERVIEW>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_diag (is, js, Cld_spec, Atmos_input, coszen, Time)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Lsctau" TYPE="real">
!   cloud optical thickness in the visible
!  </IN>
!  <IN NAME="Lsclwem" TYPE="real">
!   10 micron cloud emissivity
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_diag_stochastic (is, js, Atmos_input, coszen,       &
                                  Lsctau, Lsclwem, LscCldAmt, Time)

!--------------------------------------------------------------------
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!---------------------------------------------------------------------
 
integer,                     intent(in)   :: is,js
type(atmos_input_type),      intent(in)   :: Atmos_input
real, dimension(:,:),        intent(in)   :: coszen
real, dimension(:,:,:,:),    intent(in)   :: Lsctau, Lsclwem, LscCldAmt
type(time_type),             intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ] 
!      coszen          cosine of zenith angle [ dimensionless ]
!      Lsctau          0.6-0.685 micron cloud optical thickness 
!                      [ dimensionless ]
!      Lsclwem         Longwave emissivity [ dimensionless ]
!                      This is from 10.1-11.1 microns if the multiband 
!                      longwave emissivity is active.
!      LsCldAmt        Cloud fraction [ dimensionless ]
!                      Values should be identically 0 or 1. 
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real, dimension (size(Lsctau,1), size(Lsctau,2), size(Lsctau,3) ) ::  qv

      !
      ! Isccp histogram variables
      !
      real, dimension (size(Lsctau,1), size(Lsctau,2), 7, 7) &
                                                       :: fq_isccp

      real, dimension (size(Lsctau,1), size(Lsctau,2)) :: npoints, &
                       ninhomog, inhomogeneity_parameter
      

      integer      :: kdim
      integer      :: i, j, k
      integer, dimension(size(Lsctau,1),size(Lsctau,2)):: sunlit
      
!---------------------------------------------------------------------
!   local variables:
!
!      qv               water vapor specific humidity
!                       [ kg vapor / kg air ]
!      fq_isccp         matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges
!      npoints          flag indicating whether isccp cloud is present
!                       in column (cloud + daylight needed)
!      kdim             number of model layers
!      max_cld          greatest number of clouds in any column in the
!                       current physics window
!      i,j,k            do-loop indices
!     
!      sunlit           is the given i,j point sunlit?
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define number of model layers.
!----------------------------------------------------------------------
      kdim = size (Lsctau,3)
        
!---------------------------------------------------------------------
!    If optical properties are needed and no clouds exist in the 
!    window, call cloud_optical_properties_diag to define the cloud 
!    optical depth, the optical depth due to cloud ice and the 
!    longwave emissivity. If no clouds exist in the window, all the 
!    optical depths and emissivities are left are their initial
!    values of zero.to zero.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    Initialize ISCCP histograms
!---------------------------------------------------------------------

     npoints(:,:) = 0.
     fq_isccp(:,:,:,:) = 0.
     ninhomog(:,:) = 0.
     inhomogeneity_parameter(:,:) = 0.

!---------------------------------------------------------------------
!    Compute sunlit integer flag
!---------------------------------------------------------------------
     sunlit(:,:) = 0
     where(coszen(:,:) > 1.E-06) sunlit(:,:) = 1

!--------------------------------------------------------------------
!    define the specific humidity from the mixing ratio which has been
!    input.
!--------------------------------------------------------------------
     qv(:,:,:) = Atmos_input%cloudvapor(:,:,:)/ (1. + Atmos_input%cloudvapor(:,:,:))
               
!---------------------------------------------------------------------
!    call isccp_cloudtypes to map each model cloud to an isccp cloud
!    type, based on its optical depth and height above the surface.
!    set a flag to indicate the presence of isccp cloud in this column.
!---------------------------------------------------------------------
     call isccp_cloudtypes_stochastic (sunlit(:,:),        &
                            Atmos_input%press(:,:,1:kdim), &
                            Atmos_input%pflux(:,:,:),      &
                            qv(:,:,:),                     &
                            Atmos_input%cloudtemp(:,:,:),  &
                            Atmos_input%temp(:,:,kdim+1),  &
                            LscCldAmt(:, :, :, :),         &
                            LscTau(:, :, :, :),            &
                            Lsclwem(:, :, :, :),           &
                            fq_isccp(:,:,:,:), npoints(:,:), &
                            inhomogeneity_parameter(:,:), &
                            ninhomog(:,:))
                                                
!----------------------------------------------------------------------
!    send any desired diagnostics to the diag_manager_mod.
!----------------------------------------------------------------------
     call isccp_output (is, js, fq_isccp, npoints, &
                             inhomogeneity_parameter, ninhomog, Time)
   
!---------------------------------------------------------------------
    
    
end subroutine isccp_diag_stochastic        


!#####################################################################
! <SUBROUTINE NAME="compute_isccp_clds">
!  <OVERVIEW>
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call compute_isccp_clds (pflux, camtsw, hml_ca)
!  </TEMPLATE>
!  <IN NAME="pflux" TYPE="real">
!   average of pressure at adjacent model levels
!  </IN>
!  <IN NAME="camtsw" TYPE="real">
!   total cloud amount [ nondimensional ]
!  </IN>
!  <OUT NAME="hml_ca" TYPE="real">
!   cloud fraction for the 3 isccp cloud types
!  </OUT>
! </SUBROUTINE>
!
subroutine compute_isccp_clds (pflux, camtsw, camtsw_band, hml_ca)

!---------------------------------------------------------------------
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!--------------------------------------------------------------------- 

real,  dimension(:,:,:),   intent(in)  :: pflux, camtsw
real,  dimension(:,:,:,:), intent(in)  :: camtsw_band
real,  dimension(:,:,:),   intent(out) :: hml_ca

!---------------------------------------------------------------------
!  intent(in) variables:
!
!        pflux           average of pressure at adjacent model levels
!                        [ (kg /( m s^2) ]
!        camtsw          total cloud amount [ nondimensional ]
!        camtsw_band     total cloud amount in each sw band 
!                        [ nondimensional ]
!
!  intent(out) variable:
!
!        hml_ca          cloud fraction for the 3 isccp cloud types
!                        [ nondimensional ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer          ::   i, j, k

!---------------------------------------------------------------------
!  local variables:
!
!         mid_btm     pressure boundary between isccp middle and 
!                     isccp low clouds [ Pa ]
!         high_btm    pressure boundary between isccp middle and
!                     isccp high clouds [ Pa ]
!         i,j,k       do-loop indices
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    initialize a column integrated hi-mid-low cloud-free area array.
!--------------------------------------------------------------------
      hml_ca = 1.0
 
!---------------------------------------------------------------------
!    define arrays giving the cloud-free area in each column within the
!    pressure regionscorresponding to the ISCCP definitions of high 
!    (10-440 hPa), middle (440-680 hPa) and low (680-1000 hPa) clouds. 
!    compute high, middle and low cloud amounts assuming that independ-
!    ent clouds overlap randomly. note that model clouds above 10 hPa 
!    and below 1000 hPa are included in the totals.
!----------------------------------------------------------------------
      do k = 1,size(pflux,3)-1
        do j=1,size(pflux,2)
          do i=1,size(pflux,1)
            if (pflux(i,j,k)  <=  high_btm) then
              hml_ca(i,j,1) = hml_ca(i,j,1) * (1. - camtsw(i,j,k))
            else if ( (pflux(i,j,k) >  high_btm) .and.  &
                      (pflux(i,j,k) <=  mid_btm) ) then
              hml_ca(i,j,2) = hml_ca(i,j,2) * (1. - camtsw(i,j,k))
            else  if ( pflux(i,j,k) > mid_btm ) then
              hml_ca(i,j,3) = hml_ca(i,j,3) * (1. - camtsw(i,j,k))
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!    convert the cloud-free area to an integrated cloud fraction in 
!    the column. express the cloud area in percent.
!--------------------------------------------------------------------
      hml_ca = 1. - hml_ca
      hml_ca = 100.*hml_ca
  
!-------------------------------------------------------------------


end subroutine compute_isccp_clds



!####################################################################
! <SUBROUTINE NAME="cloud_optical_properties_diag">
!  <OVERVIEW>
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_optical_properties_diag (Cld_spec, tau, em_lw)
!  </TEMPLATE>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid
!  </IN>
!  <OUT NAME="tau" TYPE="real">
!   cloud optical depth in each of the
!                     num_slingo_bands
!  </OUT>
!  <OUT NAME="em_lw" TYPE="real">
!   longwave cloud emissivity
!  </OUT>
! </SUBROUTINE>
!
subroutine cloud_optical_properties_diag (Cld_spec, tau, em_lw)

!---------------------------------------------------------------------
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!---------------------------------------------------------------------
                              
type(cld_specification_type), intent(in)   :: Cld_spec
real, dimension(:,:,:,:),     intent(out)  :: tau
real, dimension(:,:,:),       intent(out)  :: em_lw       

!--------------------------------------------------------------------
!   intent(in) variable:
!
!      Cld_spec       cloud specification properties on model grid,
!                     [ cld_specification_type ]
!
!   intent(out) variables:
!
!      tau            cloud optical depth in each of the
!                     num_slingo_bands [ dimensionless ]
!      em_lw          longwave cloud emissivity [ dimensionless ]  
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:


      real, dimension (size(Cld_spec%lwp,1),size(Cld_spec%lwp,2),  &
                       size(Cld_spec%lwp,3), 4) ::     &
                                                tau_liq, tau_ice

      real, dimension (size(Cld_spec%lwp,1),size(Cld_spec%lwp,2),   &
                       size(Cld_spec%lwp,3)) ::     &
                                                k_liq, k_ice

!--------------------------------------------------------------------
!   local variables:
!
!       tau_liq    liquid cloud optical depth [ dimensionless ]
!       tau_ice    ice    cloud optical depth [ dimensionless ]
!       k_liq      liquid cloud mass absorption coefficient for longwave
!                  portion of the spectrum [ m**2 / kg condensate ]
!       k_ice      ice cloud mass absorption coefficient for longwave
!                  portion of the spectrum [ m**2 / kg condensate ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    compute uv cloud optical depths due to liquid. the formula for 
!    optical depth comes from: 
!    Slingo (1989), J. Atmos. Sci., vol. 46, pp. 1419-1427
!---------------------------------------------------------------------
      tau_liq(:,:,:,1) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02817 + (1.305/Cld_spec%reff_liq(:,:,:)))
      tau_liq(:,:,:,2) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02682 + (1.346/Cld_spec%reff_liq(:,:,:)))
      tau_liq(:,:,:,3) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02264 + (1.454/Cld_spec%reff_liq(:,:,:)))
      tau_liq(:,:,:,4) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.01281 + (1.641/Cld_spec%reff_liq(:,:,:)))
        
!---------------------------------------------------------------------
!    compute uv cloud optical depths due to ice. the formula for 
!    optical depth comes from:
!    Ebert and Curry (1992), J. Geophys. Res., vol. 97, pp. 3831-3836
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    IMPORTANT:  WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE 
!                BAND MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL 
!                OF SLINGO. THIS IS DONE BY COMBINING BANDS 3 and 4 OF 
!                EBERT AND CURRY. EVEN SO THE EXACT BAND LIMITS DO NOT 
!                MATCH.  FOR COMPLETENESS HERE ARE THE BAND LIMITS IN 
!                MICRONS:
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!---------------------------------------------------------------------
      tau_ice(:,:,:,1) = Cld_spec%iwp(:,:,:) * 1000. * &
                         (0.003448 + (2.431/Cld_spec%reff_ice(:,:,:)))
      tau_ice(:,:,:,2) = tau_ice(:,:,:,1)
      tau_ice(:,:,:,3) = tau_ice(:,:,:,1)
      tau_ice(:,:,:,4) = tau_ice(:,:,:,1)

!---------------------------------------------------------------------
!     back out scaling factor

      tau_liq = tau_liq / isccp_scale_factor
      tau_ice = tau_ice / isccp_scale_factor
        
!---------------------------------------------------------------------
!    compute total cloud optical depth. the mixed phase optical prop-
!    erties are based upon equation 14 of Rockel et al. 1991, 
!    Contributions to Atmospheric Physics, volume 64, pp.1-12. 
!    thus:

!          tau = tau_liq + tau_ice
!
!---------------------------------------------------------------------
      tau(:,:,:,:) = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)
        
!----------------------------------------------------------------------
!    place a minimum value on tau.
!----------------------------------------------------------------------
      where (tau(:,:,:,:) .lt. taumin)
        tau(:,:,:,:) = taumin
      end where   
        
!----------------------------------------------------------------------
!    define the  mass absorption coefficient for longwave radiation 
!    for cloud ice and cloud liquid.
!
!    NOTE THAT THE NUMBERS HERE ALREADY INCLUDE THE DIFFUSIVITY
!    FACTOR!
!----------------------------------------------------------------------
      k_liq(:,:,:) = 140.
      k_ice(:,:,:) = 4.83591 + 1758.511/Cld_spec%reff_ice(:,:,:)
        
!----------------------------------------------------------------------
!    compute combined lw emmisivity. the mixed phase optical properties
!    are based upon equation 14 of Rockel et al. 1991, Contributions to
!    Atmospheric Physics,  volume 64, pp.1-12.  thus:
!
!    transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    which can also be written as:
!
!    em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!    which is what is solved here.
!----------------------------------------------------------------------
      em_lw(:,:,:) =  1. - exp(-1.*(k_liq(:,:,:)*Cld_spec%lwp(:,:,:) + &
                                    k_ice(:,:,:)*Cld_spec%iwp(:,:,:))/ &
                                    isccp_scale_factor)

!---------------------------------------------------------------------


end subroutine cloud_optical_properties_diag


!######################################################################



                end module cloudrad_diagnostics_mod

