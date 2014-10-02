                       module rad_output_file_mod

! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Module that provides subroutines to write radiation output to
!  history file
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

!   shared modules:

use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, WARNING, close_file
use time_manager_mod,  only: time_manager_init, time_type
use diag_manager_mod,  only: register_diag_field, diag_manager_init, &
                             send_data
use constants_mod,     only: constants_init, GRAV

!  radiation package shared modules:

use rad_utilities_mod, only:  rad_utilities_init, radiative_gases_type,&
                              rad_output_type, cldrad_properties_type, &
                              cld_specification_type, atmos_input_type,&
                              Sw_control, aerosol_diagnostics_type, &
                              aerosol_type, aerosol_properties_type, &
                              surface_type, sw_output_type,  &
                              Lw_control, &
                              lw_output_type, Rad_control
use esfsw_parameters_mod, only : esfsw_parameters_init, Solar_spect

!--------------------------------------------------------------------

implicit none
private

!------------------------------------------------------------------
!    rad_output_file_mod writes an output file containing an assort-
!    ment of variables related to the sea_esf_rad radiation package.
!    this is an optionally-generated file, which may be used to sup-
!    plement the standard diagnostic model output files. NOTE THAT
!    THIS FILE IS GENERATED ONLY ON RADIATION TIMESTEPS, SO THAT WHEN 
!    SW FLUXES ARE BEING RENORMALIZED ON EACH PHYSICS STEP, VARIABLES 
!    IN THIS FILE RELATED TO SW RADIATION WILL NOT REFLECT THE EFFECTS 
!    OF THE RENORMALIZATION.
!------------------------------------------------------------------


!-------------------------------------------------------------------
!----------- version number for this module ------------------------

character(len=128)  :: version = &
'$Id: rad_output_file.f90,v 12.0 2005/04/14 15:47:24 fms Exp $'
character(len=128)  :: tagname =  '$Name: lima $'


!---------------------------------------------------------------------
!-------  interfaces --------

public   &
         rad_output_file_init, write_rad_output_file,   &
         rad_output_file_end

private   &

!   called from rad_output_file_init
        register_fields


!-------------------------------------------------------------------
!-------- namelist  ---------

logical :: write_data_file=.false.  ! data file to be written  ?


namelist  / rad_output_file_nml /  &
                                  write_data_file
     
!------------------------------------------------------------------
!----public data------


!------------------------------------------------------------------
!----private data------

!--------------------------------------------------------------------
!    DU_factor is Dobson units per (kg/m2). the first term is 
!    (avogadro's number/loeschmidt's number at STP = vol./kmol of an 
!    ideal gas at STP). second term = mol wt o3 . third term is units 
!    conversion. values are chosen from US Standard Atmospheres, 1976.
!--------------------------------------------------------------------
real, parameter  :: DU_factor =    &
                            (6.022169e26/2.68684e25)/(47.9982)*1.0e5
real, parameter  :: DU_factor2 = DU_factor/GRAV
                                   ! Dobson units per (kg/kg * dyn/cm2) 
                                   ! Dobson units per (kg/kg * N  /m2) 

!--------------------------------------------------------------------
! netcdf diagnostics field variables
!---------------------------------------------------------------------
character(len=16), parameter       :: mod_name='radiation'
real                               :: missing_value = -999.
integer, dimension(:), allocatable :: id_aerosol, id_aerosol_column
!integer, dimension(:), allocatable :: id_aerosol, id_aerosol_column, &
!                                     id_absopdep, id_absopdep_column, &
integer, dimension(:,:), allocatable :: id_absopdep,  &
                                        id_absopdep_column, &
                                      id_extopdep, id_extopdep_column
integer, dimension(2)              :: id_lw_absopdep_vlcno_column, &
                                      id_lw_extopdep_vlcno_column, &
                                      id_lwext_vlcno, id_lwssa_vlcno, &
                                      id_lwasy_vlcno, id_lw_xcoeff_vlcno
integer, dimension(2)              :: id_absopdep_vlcno_column, &
                                      id_extopdep_vlcno_column, &
                                      id_swext_vlcno, id_swssa_vlcno, &
                                      id_swasy_vlcno, id_sw_xcoeff_vlcno
integer, dimension(4)              :: id_lw_bdyflx_clr, id_lw_bdyflx, &
                                      id_sw_bdyflx_clr, id_sw_bdyflx
integer                            :: id_swheat_vlcno
integer, dimension(:), allocatable :: id_aerosol_fam, &
                                      id_aerosol_fam_column    
integer, dimension(:,:), allocatable :: id_absopdep_fam,  &
                                        id_absopdep_fam_column, &
                                        id_extopdep_fam,  &
                                        id_extopdep_fam_column
integer                            :: id_radswp, id_radp, id_temp, &
                                      id_rh2o, id_qo3, id_qo3_col,  &
                                      id_cmxolw, id_crndlw, id_flxnet, &
                                      id_fsw, id_ufsw, id_psj, &
                                      id_tmpsfc, id_cvisrfgd_dir,  &
                                      id_cirrfgd_dir, &
                                      id_cvisrfgd_dif, id_cirrfgd_dif, &
                                      id_radswpcf,  &
                                      id_cldwater, id_cldice,  &
                                      id_cldarea, &
                                      id_radpcf, id_flxnetcf, &
                                      id_fswcf, id_ufswcf, id_pressm,  &
                                      id_phalfm, id_pfluxm, &
                                      id_dphalf, id_dpflux, &
                                      id_ptop



!---------------------------------------------------------------------
!    miscellaneous variables
!---------------------------------------------------------------------
integer :: naerosol=0                      ! number of active aerosols
logical :: module_is_initialized= .false.  ! module initialized ?
integer, parameter              :: N_DIAG_BANDS = 5
character(len=16), dimension(N_DIAG_BANDS) ::   &
                     band_suffix = (/ '_vis', '_nir', '_con',  &
                                      '_bd5', '_bd6' /)


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!#####################################################################
! <SUBROUTINE NAME="rad_output_file_init">
!  <OVERVIEW>
!   Constructor of rad_output_file module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize and set up rad_output_file module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  rad_output_file_init (axes, Time, names)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="integer">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
!  <IN NAME="names" TYPE="character">
!   aerosol names
!  </IN>
! </SUBROUTINE>
!
subroutine rad_output_file_init (axes, Time, names, family_names)

!--------------------------------------------------------------------
!    rad_output_file_init is the constructor for rad_output_file_mod.
!--------------------------------------------------------------------

integer, dimension(4),           intent(in)    :: axes
type(time_type),                 intent(in)    :: Time
!character(len=64), dimension(:), intent(in)    :: names
character(len=*), dimension(:), intent(in)    :: names
character(len=*), dimension(:), intent(in)    :: family_names

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    these variables are present when running the gcm, not present
!    when running the standalone code.
!  
!       axes      diagnostic variable axes for netcdf files
!       Time      current time [ time_type(days, seconds) ]
!       names     names of active aerosols
!       family_names 
!                 names of active aerosol families
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: unit, io, ierr
      integer   :: nfields

!---------------------------------------------------------------------
!   local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!        nfields         number of active aerosol fields
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call constants_init
      call rad_utilities_init
      call esfsw_parameters_init
      call diag_manager_init  
      call time_manager_init 

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=rad_output_file_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'rad_output_file_nml')
        end do
10      call close_file (unit)
      endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                         write (stdlog(), nml=rad_output_file_nml)

!--------------------------------------------------------------------
!    if running gcm, continue on if data file is to be written. 
!--------------------------------------------------------------------
        if (write_data_file) then

!---------------------------------------------------------------------
!    register the diagnostic fields for output.
!---------------------------------------------------------------------
          nfields = size(names(:))
          call register_fields (Time, axes, nfields, names,  &
                                family_names)
        endif

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine rad_output_file_init



!################################################################
! <SUBROUTINE NAME="write_rad_output_file">
!  <OVERVIEW>
!   write_rad_output_file produces a netcdf output file containing
!   the user-specified radiation-related variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   write_rad_output_file produces a netcdf output file containing
!   the user-specified radiation-related variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call write_rad_output_file (is, ie, js, je, Atmos_input, Surface, &
!                                  Rad_output, &
!                                  Sw_output, Lw_output, Rad_gases, &
!                                  Cldrad_props, Cld_spec,  &
!                                  Time_diag, aerosol_in)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data 
!   in the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the model grid
!  </IN>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </IN>
!  <IN NAME="Rad_output" TYPE="rad_output_type">
!   rad_output_type variable containing radiation
!                        output data needed by other modules
!  </IN>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   sw_output_type variable containing shortwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid
!  </IN>
!  <IN NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid 
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing rad-
!                        iative gas input data for the radiation package
!                        on the model grid
!  </IN>
!  <IN NAME="Cldrad_pros" TYPE="cldrad_properties_type">
!   cldrad_properties_type variable containing 
!                        cloud radiative property input data for the 
!                        radiation package on the model grid
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_diagnostics_type">
!   cld_specification_type variable containing cloud
!                        microphysical data
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!   time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ]
!  </IN>
!  <IN NAME="aerosol_in" TYPE="real">
!   optional aerosol data
!  </IN>
! </SUBROUTINE>
!
subroutine write_rad_output_file (is, ie, js, je, Atmos_input, Surface,&
                                  Rad_output, Sw_output, Lw_output,    &
                                  Rad_gases, Cldrad_props, Cld_spec,   &
                                  Time_diag, Aerosol, Aerosol_props, &
                                  Aerosol_diags)

!----------------------------------------------------------------
!    write_rad_output_file produces a netcdf output file containing
!    the user-specified radiation-related variables.
!----------------------------------------------------------------

integer,                      intent(in)            ::  is, ie, js, je
type(atmos_input_type),       intent(in)            ::  Atmos_input
type(surface_type),           intent(in)            ::  Surface
type(rad_output_type),        intent(in)            ::  Rad_output
type(sw_output_type),         intent(in)            ::  Sw_output
type(lw_output_type),         intent(in)            ::  Lw_output
type(radiative_gases_type),   intent(in)            ::  Rad_gases
type(cldrad_properties_type), intent(in)            ::  Cldrad_props
type(cld_specification_type), intent(in)            ::  Cld_spec      
type(time_type),              intent(in)            ::  Time_diag
type(aerosol_type),           intent(in), optional  ::  Aerosol
type(aerosol_properties_type), intent(in), optional ::  Aerosol_props
type(aerosol_diagnostics_type), intent(in), optional :: Aerosol_diags

!------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      Atmos_input       atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the model grid
!      Surface           surface input fields to radiation package
!                        [ surface_type ]
!      Rad_output        rad_output_type variable containing radiation
!                        output data needed by other modules
!      Sw_output         sw_output_type variable containing shortwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid
!      Lw_output         lw_output_type variable containing longwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid 
!      Rad_gases         radiative_gases_type variable containing rad-
!                        iative gas input data for the radiation package
!                        on the model grid
!      Cldrad_props      cldrad_properties_type variable containing 
!                        cloud radiative property input data for the 
!                        radiation package on the model grid
!      Cld_spec          cld_specification_type variable containing 
!                        cloud specification input data for the 
!                        radiation package on the model grid
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ]  
!
!  intent(in), optional variables:
!
!      aerosol_in        active aerosol distributions
!                        [ kg / m**2 ]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real, dimension (size(Atmos_input%press,1),   &
                       size(Atmos_input%press,2) )  ::  &
                           tmpsfc, psj, cvisrfgd_dir, cirrfgd_dir,  &
                           cvisrfgd_dif, cirrfgd_dif, qo3_col, &
                           ptop

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3)  ) ::   &
                          fsw, ufsw, fswcf, ufswcf, flxnet, flxnetcf, &
                          phalfm, pfluxm

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3)-1 ) ::  &
                       temp, rh2o, qo3, cmxolw, crndlw, radp, radswp, &
                       radpcf, radswpcf, pressm, dphalf, dpflux, deltaz

      real, dimension(:,:,:),   allocatable :: aerosol_col
      real, dimension(:,:,:,:),   allocatable :: extopdep_col
      real, dimension(:,:,:,:),   allocatable :: absopdep_col
      real, dimension(:,:,:),   allocatable :: lw_extopdep_vlcno_col
      real, dimension(:,:,:),   allocatable :: lw_absopdep_vlcno_col
      real, dimension(:,:,:),   allocatable :: extopdep_vlcno_col
      real, dimension(:,:,:),   allocatable :: absopdep_vlcno_col
      real, dimension(:,:,:,:),   allocatable :: absopdep_fam_col
      real, dimension(:,:,:,:),   allocatable :: extopdep_fam_col
      real, dimension(:,:,:),   allocatable :: aerosol_fam_col
      real, dimension(:,:,:,:,:), allocatable :: absopdep_fam
      real, dimension(:,:,:,:,:), allocatable :: extopdep_fam
      real, dimension(:,:,:,:), allocatable :: aerosol_fam

      logical   :: used  
      integer   :: kerad ! number of model layers
      integer   :: n, k, na, nfamilies, nl
      integer   :: nv, vis_indx, nir_indx
      integer   :: co_indx, bnd_indx

!----------------------------------------------------------------------
!  local variables:
!
!      tmpsfc         surface temperature [ deg K ]
!      psj            surface pressure [ hPa ]
!      cvisrfgd       surface visible light albedo [ dimensionless ]
!      cirrfgd        surface ir albedo [ dimensionless ]
!      tot_clds       total column isccp clouds [ percent ]
!      cld_isccp_hi   number of isccp high clouds [ percent ]
!      cld_isccp_mid  number of isccp middle clouds [ percent ]
!      cld_isccp_low  number of isccp low clouds [ percent ]
!      qo3_col        ozone column [ DU ]
!      fsw            net shortwave flux [ W / m**2 ]
!      ufsw           upward shortwave flux [ W / m**2 ]
!      fswcf          net sw flux in the absence of clouds [ W / m**2 ]
!      ufswcf         upward sw flux in absence of clouds [ W / m**2]
!      flxnet         net longwave flux [ W / m**2 ]
!      flxnetcf       net lw flux in the absence of clouds [ W / m**2 ]
!      phalfm         model interface level pressure [ Pa ]
!      pfluxm         avg of adjacent model level pressures [ Pa ]
!      temp           temperature [ deg K ]
!      rh2o           water vapor specific humidity [ g / g ]
!      qo3            ozone mixing ratio [ g / g ]
!      heatra         lw heating rate [ deg K / day ]
!      heatracf       lw heating rate without cloud [ deg K / day ]
!      cmxolw         amount of maximal overlap clouds [ percent]
!      crndlw         amount of ramndom overlap clouds [ percent]
!      radp           lw + sw heating rate [ deg K / sec ]
!      radswp         sw heating rate [ deg K / sec ]
!      radpcf         lw + sw heating rate w/o clouds [ deg K / sec ]
!      radswpcf       sw heating rate w/o clouds [ deg K / sec ]
!      pressm         pressure at model levels [ Pa ]
!      aerosol_col
!      used
!      kerad
!      n,k
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_output_file_mod', &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    if the file is not to be written, do nothing.
!--------------------------------------------------------------------
      if (write_data_file) then

!--------------------------------------------------------------------
!    if the file is to be written, define the number of model layers.
!--------------------------------------------------------------------
        kerad = ubound(Atmos_input%temp,3) - 1

!--------------------------------------------------------------------
!    retrieve the desired fields from the input derived data types.
!--------------------------------------------------------------------
        tmpsfc(:,:)     = Atmos_input%temp(:,:,kerad+1)
        psj   (:,:)     = Atmos_input%press(:,:,kerad+1)
        pressm(:,:,:)   = Atmos_input%press(:,:,1:kerad)
        phalfm(:,:,:)   = Atmos_input%phalf(:,:,:)
        pfluxm(:,:,:)   = Atmos_input%pflux(:,:,:)
        temp(:,:,:)     = Atmos_input%temp(:,:,1:kerad)
        rh2o(:,:,:)     = Atmos_input%rh2o(:,:,:)
        radp(:,:,:)     = Rad_output%tdt_rad(is:ie,js:je,:)
        radswp(:,:,:)   = Rad_output%tdtsw  (is:ie,js:je,:)
!       cirrfgd(:,:)    = Surface%asfc(:,:)
!       cvisrfgd(:,:)   = Surface%asfc(:,:)
        cirrfgd_dir(:,:)    = Surface%asfc_nir_dir(:,:)
        cvisrfgd_dir(:,:)   = Surface%asfc_vis_dir(:,:)
        cirrfgd_dif(:,:)    = Surface%asfc_nir_dif(:,:)
        cvisrfgd_dif(:,:)   = Surface%asfc_vis_dif(:,:)
        fsw(:,:,:)      = Sw_output% fsw(:,:,:)
        ufsw(:,:,:)     = Sw_output%ufsw(:,:,:)
        flxnet(:,:,:)   = Lw_output%flxnet(:,:,:)
        qo3(:,:,:)      = Rad_gases%qo3(:,:,:)
        cmxolw(:,:,:)   = 100.0*Cld_spec%cmxolw(:,:,:)
        crndlw(:,:,:)   = 100.0*Cld_spec%crndlw(:,:,:)
        deltaz(:,:,:)   = Atmos_input%deltaz(:,:,:)

          do k = 1,kerad
            dphalf(:,:,k)   = phalfm(:,:,k+1) - phalfm(:,:,k)
            dpflux(:,:,k)   = pfluxm(:,:,k+1) - pfluxm(:,:,k)
          enddo
          ptop(:,:) = 0.01*Atmos_input%phalf(:,:,1)

        if (Rad_control%do_totcld_forcing) then 
          fswcf(:,:,:)    = Sw_output%fswcf(:,:,:)
          ufswcf(:,:,:)   = Sw_output%ufswcf(:,:,:)
          flxnetcf(:,:,:) = Lw_output%flxnetcf(:,:,:)
          radpcf(:,:,:)   = Rad_output%tdt_rad_clr(is:ie,js:je,:)
          radswpcf(:,:,:) = Rad_output%tdtsw_clr  (is:ie,js:je,:)
        endif

!---------------------------------------------------------------------
!    calculate the column ozone in DU (Dobson units). convert from 
!    (kg/kg) * (N/m2) to DU (1DU = 2.687E16 molec cm^-2).
!---------------------------------------------------------------------
        qo3_col(:,:) = 0.
        do k = 1,size(qo3,3)
          qo3_col(:,:) = qo3_col(:,:) + &
                         qo3(:,:,k)*(Atmos_input%pflux(:,:,k+1) -  &
                                     Atmos_input%pflux(:,:,k))
        end do
        qo3_col(:,:) = qo3_col(:,:)*DU_factor2

!---------------------------------------------------------------------
!    define the aerosol fields and calculate the column aerosol. 
!---------------------------------------------------------------------
        if (Rad_control%do_aerosol) then
          allocate ( aerosol_col(size(Aerosol%aerosol, 1), &
                                 size(Aerosol%aerosol, 2), &
                                 size(Aerosol%aerosol, 4)) )       
          aerosol_col(:,:,:) = SUM (Aerosol%aerosol(:,:,:,:), 3)
!         if (Sw_control%do_swaerosol) then
            allocate ( extopdep_col(size(Aerosol_diags%extopdep  , 1), &
                                    size(Aerosol_diags%extopdep  , 2), &
                               size(Aerosol_diags%extopdep  , 4), &
                                    N_DIAG_BANDS) )
            extopdep_col(:,:,:,:) =    &
                          SUM (Aerosol_diags%extopdep  (:,:,:,:,:), 3)
            allocate ( absopdep_col(size(Aerosol_diags%absopdep  , 1), &
                                    size(Aerosol_diags%absopdep  , 2), &
                              size(Aerosol_diags%absopdep  , 4), &
                                   N_DIAG_BANDS) )
            absopdep_col(:,:,:,:) =    &
                           SUM (Aerosol_diags%absopdep  (:,:,:,:,:), 3)
          if (Sw_control%do_swaerosol) then
            if (Rad_control%volcanic_sw_aerosols) then
              allocate ( extopdep_vlcno_col(   &
                           size(Aerosol_diags%extopdep_vlcno  , 1), &
                           size(Aerosol_diags%extopdep_vlcno  , 2),2))
              extopdep_vlcno_col(:,:,:) =    &
                        SUM (Aerosol_diags%extopdep_vlcno  (:,:,:,:), 3)
              allocate ( absopdep_vlcno_col(    &
                            size(Aerosol_diags%absopdep_vlcno  , 1), &
                            size(Aerosol_diags%absopdep_vlcno  , 2),2))
              absopdep_vlcno_col(:,:,:) =    &
                        SUM (Aerosol_diags%absopdep_vlcno  (:,:,:,:), 3)
            endif
              
          endif
          if (Lw_control%do_lwaerosol) then
            if (Rad_control%volcanic_lw_aerosols) then
              allocate ( lw_extopdep_vlcno_col(   &
                         size(Aerosol_diags%lw_extopdep_vlcno  , 1), &
                         size(Aerosol_diags%lw_extopdep_vlcno  , 2),2))
              lw_extopdep_vlcno_col(:,:,:) =    &
                    SUM (Aerosol_diags%lw_extopdep_vlcno  (:,:,:,:), 3)
              allocate ( lw_absopdep_vlcno_col(    &
                          size(Aerosol_diags%lw_absopdep_vlcno  , 1), &
                          size(Aerosol_diags%lw_absopdep_vlcno  , 2),2))
              lw_absopdep_vlcno_col(:,:,:) =    &
                    SUM (Aerosol_diags%lw_absopdep_vlcno  (:,:,:,:), 3)
            endif
          endif
        endif
        
!---------------------------------------------------------------------
!    define the aerosol family output fields.
!---------------------------------------------------------------------
        if (Rad_control%do_aerosol) then
          nfamilies = size(Aerosol%family_members,2)
          if (nfamilies > 0) then
            allocate (aerosol_fam (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    size(Aerosol%aerosol,3), &
                                    nfamilies))
            allocate (aerosol_fam_col (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    nfamilies))
            allocate (extopdep_fam (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    size(Aerosol%aerosol,3), &
                                    nfamilies, N_DIAG_BANDS))
            allocate (absopdep_fam (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    size(Aerosol%aerosol,3), &
                                    nfamilies, N_DIAG_BANDS))
            allocate (extopdep_fam_col (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    nfamilies, N_DIAG_BANDS))
            allocate (absopdep_fam_col (     &
                                    size(Aerosol%aerosol,1), &
                                    size(Aerosol%aerosol,2), &
                                    nfamilies, N_DIAG_BANDS))
            aerosol_fam = 0.
            aerosol_fam_col = 0.
            extopdep_fam = 0.
            absopdep_fam = 0.
            extopdep_fam_col = 0.
            absopdep_fam_col = 0.
            do n = 1, nfamilies                      
              do na = 1, naerosol                
                if (Aerosol%family_members(na,n)) then
                  aerosol_fam(:,:,:,n) = aerosol_fam(:,:,:,n) +  &
                                         Aerosol%aerosol(:,:,:,na)
                  aerosol_fam_col(:,:,n) = aerosol_fam_col(:,:,n) +  &
                                         aerosol_col(:,:,na)
                  do nl = 1,N_DIAG_BANDS
                extopdep_fam(:,:,:,n,nl) = extopdep_fam(:,:,:,n,nl) +  &
                                    Aerosol_diags%extopdep(:,:,:,na,nl)
            extopdep_fam_col(:,:,n,nl) = extopdep_fam_col(:,:,n,nl) +  &
                                      extopdep_col(:,:,na,nl)
            absopdep_fam(:,:,:,n,nl) = absopdep_fam(:,:,:,n,nl) +  &
                                    Aerosol_diags%absopdep(:,:,:,na,nl)
            absopdep_fam_col(:,:,n,nl) = absopdep_fam_col(:,:,n,nl) +  &
                                      absopdep_col(:,:,na,nl)
                  end do ! (nl)
                endif
              end do

              if (Aerosol%family_members(naerosol+1,n)) then
           if (Sw_control%do_swaerosol) then
             if (Rad_control%volcanic_sw_aerosols) then
             extopdep_fam_col(:,:,n,1) = extopdep_fam_col(:,:,n,1) +  &
                                       extopdep_vlcno_col(:,:,1)
             absopdep_fam_col(:,:,n,1) = absopdep_fam_col(:,:,n,1) +  &
                                       absopdep_vlcno_col(:,:,1)
             extopdep_fam_col(:,:,n,2) = extopdep_fam_col(:,:,n,2) +  &
                                       extopdep_vlcno_col(:,:,2)
             absopdep_fam_col(:,:,n,2) = absopdep_fam_col(:,:,n,2) +  &
                                       absopdep_vlcno_col(:,:,2)
           endif
          endif
          if (Lw_control%do_lwaerosol) then
            if (Rad_control%volcanic_lw_aerosols) then
             extopdep_fam_col(:,:,n,4) = extopdep_fam_col(:,:,n,4) +  &
                                    lw_extopdep_vlcno_col(:,:,1)
             absopdep_fam_col(:,:,n,4) = absopdep_fam_col(:,:,n,4) +  &
                                    lw_absopdep_vlcno_col(:,:,1)
             extopdep_fam_col(:,:,n,5) = extopdep_fam_col(:,:,n,5) +  &
                                    lw_extopdep_vlcno_col(:,:,2)
             absopdep_fam_col(:,:,n,5) = absopdep_fam_col(:,:,n,5) +  &
                                    lw_absopdep_vlcno_col(:,:,2)
            endif
           endif
          endif
           end do

          do n = 1,nfamilies
            if (id_aerosol_fam(n)  > 0 ) then
              used = send_data (id_aerosol_fam(n),  &
                                aerosol_fam(:,:,:,n)/deltaz(:,:,:),   &
                                Time_diag, is, js, 1)
            endif
            if (id_aerosol_fam_column(n)  > 0 ) then
              used = send_data (id_aerosol_fam_column(n),     &
                                aerosol_fam_col(:,:,n), Time_diag, is, js)
            endif
!           if (Sw_control%do_swaerosol) then
            do nl=1,N_DIAG_BANDS
              if (id_extopdep_fam(n,nl)  > 0 ) then
                used = send_data (id_extopdep_fam(n,nl),    &
                                  extopdep_fam  (:,:,:,n,nl), &
                                  Time_diag, is, js, 1)
              endif
              if (id_extopdep_fam_column(n,nl)  > 0 ) then
                used = send_data (id_extopdep_fam_column(n,nl),     &
                                 extopdep_fam_col(:,:,n,nl), Time_diag, is, js)
              endif
              if (id_absopdep_fam(n,nl)  > 0 ) then
                used = send_data (id_absopdep_fam(n,nl),    &
                                  absopdep_fam  (:,:,:,n,nl), &
                                  Time_diag, is, js, 1)
              endif
              if (id_absopdep_fam_column(n,nl)  > 0 ) then
                used = send_data (id_absopdep_fam_column(n,nl),     &
                                 absopdep_fam_col(:,:,n,nl), Time_diag, is, js)
              endif
            end do  
!           endif
        end do
          deallocate (aerosol_fam)
          deallocate (aerosol_fam_col)
!         if (Sw_control%do_swaerosol) then
            deallocate (extopdep_fam)
            deallocate (absopdep_fam)
            deallocate (extopdep_fam_col)
            deallocate (absopdep_fam_col)
!         endif
      endif
    endif
        
!---------------------------------------------------------------------
!    send the user-designated data to diag_manager_mod for processing.
!---------------------------------------------------------------------
        if (id_radswp > 0 ) then
          used = send_data (id_radswp, radswp, Time_diag, is, js, 1)
        endif

        if (id_radp > 0 ) then
          used = send_data (id_radp, radp, Time_diag, is, js, 1)
        endif

        if (id_temp > 0 ) then
          used = send_data (id_temp, temp, Time_diag, is, js, 1)
        endif

        if (id_pressm > 0 ) then
          used = send_data (id_pressm, pressm, Time_diag, is, js, 1)
        endif

        if (id_phalfm > 0 ) then
          used = send_data (id_phalfm, phalfm, Time_diag, is, js, 1)
        endif
 
        if (id_pfluxm > 0 ) then
          used = send_data (id_pfluxm, pfluxm, Time_diag, is, js, 1)
        endif

        if (id_rh2o > 0 ) then
          used = send_data (id_rh2o, rh2o, Time_diag, is, js, 1)
        endif

        if (id_cldwater > 0 ) then
          used = send_data (id_cldwater, Cld_spec%cloud_water,  &
                            Time_diag, is, js, 1)
        endif

        if (id_cldice > 0 ) then
          used = send_data (id_cldice, Cld_spec%cloud_ice,  &
                            Time_diag, is, js, 1)
        endif

        if (id_cldarea > 0 ) then
          used = send_data (id_cldarea, Cld_spec%cloud_area,  &
                            Time_diag, is, js, 1)
        endif

        if (id_qo3  > 0 ) then
          used = send_data (id_qo3, qo3, Time_diag, is, js, 1)
        endif

        if (id_qo3_col  > 0 ) then
          used = send_data (id_qo3_col, qo3_col, Time_diag, is, js)
        endif

          if (id_dphalf > 0 ) then
            used = send_data (id_dphalf, dphalf, Time_diag, is, js, 1)
          endif

          if (id_dpflux > 0 ) then
            used = send_data (id_dpflux, dpflux, Time_diag, is, js, 1)
          endif

          if (id_ptop  > 0 ) then
            used = send_data (id_ptop, ptop, Time_diag, is, js)
          endif
        if (Rad_control%do_aerosol) then
          do n = 1,naerosol
            if (id_aerosol(n)  > 0 ) then
              used = send_data (id_aerosol(n),  &
                                Aerosol%aerosol(:,:,:,n)/deltaz(:,:,:),&
                                Time_diag, is, js, 1)
            endif
            if (id_aerosol_column(n)  > 0 ) then
              used = send_data (id_aerosol_column(n),     &
                                aerosol_col(:,:,n), Time_diag, is, js)
            endif
!           if (Sw_control%do_swaerosol) then
            do nl=1,N_DIAG_BANDS
              if (id_extopdep(n,nl)  > 0 ) then
                used = send_data (id_extopdep(n,nl),    &
                                  Aerosol_diags%extopdep  (:,:,:,n,nl), &
                                  Time_diag, is, js, 1)
              endif
              if (id_extopdep_column(n,nl)  > 0 ) then
                used = send_data (id_extopdep_column(n,nl),     &
                                 extopdep_col(:,:,n,nl), Time_diag, is, js)
              endif
              if (id_absopdep(n,nl)  > 0 ) then
                used = send_data (id_absopdep(n,nl),    &
                                  Aerosol_diags%absopdep  (:,:,:,n,nl), &
                                  Time_diag, is, js, 1)
              endif
              if (id_absopdep_column(n,nl)  > 0 ) then
                used = send_data (id_absopdep_column(n,nl),     &
                                 absopdep_col(:,:,n,nl), Time_diag, is, js)
              endif
!           endif
            end do
          end do
          if (Rad_control%volcanic_lw_aerosols) then
!           co_indx = size(Aerosol_props%lw_ext,4)
            co_indx = 5
            if (id_lwext_vlcno(1)  > 0 ) then
              used = send_data (id_lwext_vlcno(1),     &
                         Aerosol_props%lw_ext(:,:,:,co_indx)*  &
                         Atmos_input%deltaz(:,:,:),  &
                         Time_diag, is, js,1)
            endif
            if (id_lw_xcoeff_vlcno(1)  > 0 ) then
              used = send_data (id_lw_xcoeff_vlcno(1),     &
                            Aerosol_props%lw_ext(:,:,:,co_indx),  &
                            Time_diag, is, js,1)
            endif
            if (id_lwssa_vlcno(1)  > 0 ) then
              used = send_data (id_lwssa_vlcno(1),     &
                            Aerosol_props%lw_ssa(:,:,:,co_indx),  &
                            Time_diag, is, js,1)
            endif
            if (id_lwasy_vlcno(1)  > 0 ) then
              used = send_data (id_lwasy_vlcno(1),     &
                           Aerosol_props%lw_asy(:,:,:,co_indx),  &
                           Time_diag, is, js,1)
            endif
!           bnd_indx = 4
            bnd_indx = 6
            if (id_lwext_vlcno(2)  > 0 ) then
              used = send_data (id_lwext_vlcno(2),     &
                               Aerosol_props%lw_ext(:,:,:,bnd_indx)* &
                               Atmos_input%deltaz(:,:,:),  &
                               Time_diag, is, js,1)
            endif
            if (id_lw_xcoeff_vlcno(2)  > 0 ) then
              used = send_data (id_lw_xcoeff_vlcno(2),     &
                               Aerosol_props%lw_ext(:,:,:,bnd_indx),  &
                               Time_diag, is, js,1)
            endif
            if (id_lwssa_vlcno(2)  > 0 ) then
              used = send_data (id_lwssa_vlcno(2),     &
                               Aerosol_props%lw_ssa(:,:,:,bnd_indx), &
                               Time_diag, is, js,1)
            endif
            if (id_lwasy_vlcno(2)  > 0 ) then
              used = send_data (id_lwasy_vlcno(2),     &
                              Aerosol_props%lw_asy(:,:,:,bnd_indx), &
                              Time_diag, is, js,1)
            endif
            do nv=1,2
              if (id_lw_extopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_lw_extopdep_vlcno_column(nv),     &
                                  lw_extopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
              if (id_lw_absopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_lw_absopdep_vlcno_column(nv),     &
                                  lw_absopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
            end do
            deallocate (lw_absopdep_vlcno_col, lw_extopdep_vlcno_col)
          endif
          if (Rad_control%volcanic_sw_aerosols) then
            vis_indx = Solar_spect%visible_band_indx
            if (id_swext_vlcno(1)  > 0 ) then
              used = send_data (id_swext_vlcno(1),     &
                                Aerosol_props%sw_ext(:,:,:,vis_indx)* &
                                Atmos_input%deltaz(:,:,:),  &
                                Time_diag, is, js,1)
            endif
            if (id_sw_xcoeff_vlcno(1)  > 0 ) then
              used = send_data (id_sw_xcoeff_vlcno(1),     &
                               Aerosol_props%sw_ext(:,:,:,vis_indx),  &
                               Time_diag, is, js,1)
            endif
            if (id_swssa_vlcno(1)  > 0 ) then
              used = send_data (id_swssa_vlcno(1),     &
                                 Aerosol_props%sw_ssa(:,:,:,vis_indx), &
                                 Time_diag, is, js,1)
            endif
            if (id_swasy_vlcno(1)  > 0 ) then
              used = send_data (id_swasy_vlcno(1),     &
                                Aerosol_props%sw_asy(:,:,:,vis_indx), &
                                Time_diag, is, js,1)
            endif
            if (id_swheat_vlcno  > 0 ) then
              used = send_data (id_swheat_vlcno ,    &
                                Aerosol_diags%sw_heating_vlcno(:,:,:), &
                                Time_diag, is, js, 1)
            endif
            nir_indx = Solar_spect%one_micron_indx
            if (id_swext_vlcno(2)  > 0 ) then
              used = send_data (id_swext_vlcno(2),     &
                               Aerosol_props%sw_ext(:,:,:,nir_indx)*  &
                               Atmos_input%deltaz(:,:,:),  &
                               Time_diag, is, js,1)
            endif
            if (id_sw_xcoeff_vlcno(2)  > 0 ) then
              used = send_data (id_sw_xcoeff_vlcno(2),     &
                               Aerosol_props%sw_ext(:,:,:,nir_indx),  &
                               Time_diag, is, js,1)
            endif
            if (id_swssa_vlcno(2)  > 0 ) then
              used = send_data (id_swssa_vlcno(2),     &
                                Aerosol_props%sw_ssa(:,:,:,nir_indx),  &
                                Time_diag, is, js,1)
            endif
            if (id_swasy_vlcno(2)  > 0 ) then
              used = send_data (id_swasy_vlcno(2),     &
                                Aerosol_props%sw_asy(:,:,:,vis_indx), &
                                Time_diag, is, js,1)
            endif
            do nv=1,2
              if (id_extopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_extopdep_vlcno_column(nv),     &
                                  extopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
              if (id_absopdep_vlcno_column(nv)  > 0 ) then
                used = send_data (id_absopdep_vlcno_column(nv),     &
                                  absopdep_vlcno_col(:,:,nv),  &
                                  Time_diag, is, js)
              endif
            end do
            deallocate (absopdep_vlcno_col, extopdep_vlcno_col)
          endif
          deallocate (aerosol_col)
          if (Sw_control%do_swaerosol) then
            deallocate (extopdep_col)
            deallocate (absopdep_col)
          endif
        endif

        if (id_cmxolw > 0 ) then
          used = send_data (id_cmxolw, cmxolw, Time_diag, is, js, 1)
        endif

        if (id_crndlw > 0 ) then
          used = send_data (id_crndlw, crndlw, Time_diag, is, js, 1)
        endif

        if (id_flxnet > 0 ) then
          used = send_data (id_flxnet, flxnet, Time_diag, is, js, 1)
        endif

        if (id_fsw    > 0 ) then
          used = send_data (id_fsw   , fsw   , Time_diag, is, js, 1)
        endif

        if (id_ufsw   > 0 ) then
          used = send_data (id_ufsw  , ufsw  , Time_diag, is, js, 1)
        endif

        if (id_psj    > 0 ) then
          used = send_data (id_psj   , psj   , Time_diag, is, js)
        endif

        if (id_tmpsfc > 0 ) then
          used = send_data (id_tmpsfc, tmpsfc, Time_diag, is, js)
        endif

        if (id_cvisrfgd_dir > 0 ) then
          used = send_data (id_cvisrfgd_dir, cvisrfgd_dir, Time_diag, is, js)
        endif
 
        if (id_cvisrfgd_dif > 0 ) then
          used = send_data (id_cvisrfgd_dif, cvisrfgd_dif, Time_diag, is, js)
        endif

        if (id_cirrfgd_dir > 0 ) then
          used = send_data (id_cirrfgd_dir , cirrfgd_dir , Time_diag, is,js)
        endif
 
        if (id_cirrfgd_dif > 0 ) then
          used = send_data (id_cirrfgd_dif , cirrfgd_dif , Time_diag, is,js)
        endif

     do n=1, 4
        if (id_lw_bdyflx(n) > 0 ) then
          used = send_data (id_lw_bdyflx(n) , Lw_output%bdy_flx(:,:,n),&
                            Time_diag, is,js)
        endif
        if (id_sw_bdyflx(n) > 0 ) then
          used = send_data (id_sw_bdyflx(n) , Sw_output%bdy_flx(:,:,n),&
                            Time_diag, is,js)
        endif
     end do

        if (Rad_control%do_totcld_forcing) then
     do n=1, 4
        if (id_lw_bdyflx_clr(n) > 0 ) then
          used = send_data (id_lw_bdyflx_clr(n) ,   &
                            Lw_output%bdy_flx_clr(:,:,n),&
                            Time_diag, is,js)
        endif
        if (id_sw_bdyflx_clr(n) > 0 ) then
          used = send_data (id_sw_bdyflx_clr(n) ,   &
                            Sw_output%bdy_flx_clr(:,:,n),&
                            Time_diag, is,js)
        endif
     end do

          if (id_radswpcf > 0 ) then
            used = send_data (id_radswpcf, radswpcf, Time_diag,   &
                              is, js, 1)
          endif

          if (id_radpcf > 0 ) then
            used = send_data (id_radpcf, radpcf, Time_diag, is, js, 1)
          endif

          if (id_flxnetcf > 0 ) then
            used = send_data (id_flxnetcf, flxnetcf, Time_diag,   &
                              is, js, 1)
          endif

          if (id_fswcf  > 0 ) then
            used = send_data (id_fswcf , fswcf , Time_diag, is, js, 1)
          endif

          if (id_ufswcf  > 0 ) then
            used = send_data (id_ufswcf , ufswcf , Time_diag, is, js, 1)
          endif
        endif
      endif

!------------------------------------------------------------------


end subroutine write_rad_output_file



!#####################################################################
! <SUBROUTINE NAME="rad_output_file_end">
!  <OVERVIEW>
!   rad_output_file_end is the destructor for rad_output_file_mod
!  </OVERVIEW>
!  <DESCRIPTION>
!   rad_output_file_end is the destructor for rad_output_file_mod
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rad_output_file_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rad_output_file_end

!-------------------------------------------------------------------
!    rad_output_file_end is the destructor for rad_output_file_mod.
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_output_file_mod', &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized= .false. 

!----------------------------------------------------------------------

end subroutine rad_output_file_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#################################################################
! <SUBROUTINE NAME="register_fields">
!  <OVERVIEW>
!   register_fields send the relevant information concerning the 
!    user-desired output fields to diag_manager_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   register_fields send the relevant information concerning the 
!    user-desired output fields to diag_manager_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call register_fields (Time, axes, nfilds, names)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="nfields" TYPE="integer">
!   number of aerosol fields
!  </IN>
!  <IN NAME="names" TYPE="character">
!   names of aerosol fields
!  </IN>
! </SUBROUTINE>
!
subroutine register_fields (Time, axes, nfields, names, family_names)

!--------------------------------------------------------------------
!    register_fields send the relevant information concerning the 
!    user-desired output fields to diag_manager_mod.
!--------------------------------------------------------------------

type(time_type),                 intent(in) :: Time
integer, dimension(4),           intent(in) :: axes
integer,                         intent(in) :: nfields
!character(len=64), dimension(:), intent(in) :: names
character(len=*), dimension(:), intent(in) :: names, family_names

!--------------------------------------------------------------------
!  intent(in) variables:
!
!       Time      current time [ time_type(days, seconds) ]
!       axes      diagnostic variable axes for netcdf files
!       nfields   number of active aerosol species
!       names     names of active aerosol species
!       family_names  
!                 names of active aerosol families
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      character(len=64), dimension(:), allocatable ::   & 
                                               aerosol_column_names, &
                                              extopdep_column_names, &
                                              absopdep_column_names, &
                                                 extopdep_names, &
                                                 absopdep_names, &
                                          aerosol_fam_column_names, &
                                         extopdep_fam_column_names, &
                                        absopdep_fam_column_names, &
                                        extopdep_fam_names, &
                                        absopdep_fam_names
      integer, dimension(4)    :: bxes
      integer                  :: n, nl
      integer                  :: nfamilies
      real                     :: trange(2)

!---------------------------------------------------------------------
!   local variables:
!
!       aerosol_column_names
!       bxes                   diagnostic variable axes with elements 
!                              (1:3) valid for variables defined at
!                              flux levels
!       n
!
!--------------------------------------------------------------------
 
!-------------------------------------------------------------------
!    define variable axis array with elements (1:3) valid for variables
!    defined at flux levels.
!-------------------------------------------------------------------
      bxes(1:2) = axes(1:2)
      bxes(3) = axes(4)
      bxes(4) = axes(4)
      trange =(/ 100., 400. /)

!---------------------------------------------------------------------
!    register the potential diagnostic variables from this module.
!    the variables will actually be saved and then output only if
!    they are activated through the input diag_file.
!---------------------------------------------------------------------
      id_radswp = &
         register_diag_field (mod_name, 'radswp', axes(1:3), Time, &
                          'temperature tendency for SW radiation', &
                          'deg_K/sec', missing_value=missing_value)

      id_radp = &
         register_diag_field (mod_name, 'radp', axes(1:3), Time, &
                          'temperature tendency for radiation', &
                          'deg_K/sec', missing_value=missing_value)

      id_temp   = &
         register_diag_field (mod_name, 'temp', axes(1:3), Time, &
                          'temperature field', &
                          'deg_K', missing_value=trange(1), &
                          range=trange)

      id_pressm  = &
         register_diag_field (mod_name, 'pressm', axes(1:3), Time, &
                           'model level pressure', &
                           'Pa', missing_value=missing_value)
 
      id_phalfm  = &
         register_diag_field (mod_name, 'phalfm', bxes(1:3), Time, &
                           'model interface level pressure', &
                           'Pa', missing_value=missing_value)

      id_pfluxm  = &
          register_diag_field (mod_name, 'pfluxm', bxes(1:3), Time, &
                           'radiation flux level pressures', &
                           'Pa', missing_value=missing_value)

        id_dpflux  = &
          register_diag_field (mod_name, 'dpflux', axes(1:3), Time,  &
                           'radiation flux layer thickness', &
                           'hPa', missing_value=missing_value)
   
        id_dphalf  = &
          register_diag_field (mod_name, 'dphalf', axes(1:3), Time,  &
                           'radiation model layer thickness', &
                           'hPa', missing_value=missing_value)
   
        id_ptop  = &
          register_diag_field (mod_name, 'ptop', axes(1:2), Time, &
                           'pressure at model top', &
                           'hPa', missing_value=missing_value)
 

      id_rh2o   = &
         register_diag_field (mod_name, 'rh2o', axes(1:3), Time, &
                          'water vapor mixing ratio', &
                          'kg/kg', missing_value=missing_value)

      id_cldwater = &
         register_diag_field (mod_name, 'cloud_water', axes(1:3), Time,&
                          'cloud water specific humidity', &
                          'kg/kg', missing_value=missing_value)

      id_cldice = &
         register_diag_field (mod_name, 'cloud_ice', axes(1:3), Time,&
                          'cloud ice specific humidity', &
                          'kg/kg', missing_value=missing_value)

      id_cldarea = &
         register_diag_field (mod_name, 'cloud_area', axes(1:3), Time,&
                          'cloud fractional area', &
                          'fraction', missing_value=missing_value)

      id_qo3    = &
         register_diag_field (mod_name, 'qo3', axes(1:3), Time, &
                          'ozone mixing ratio', &
                          'kg/kg', missing_value=missing_value)

      id_qo3_col = &
         register_diag_field (mod_name, 'qo3_col', axes(1:2), Time, &
                          'ozone column', &
                          'DU', missing_value=missing_value)

!--------------------------------------------------------------------
!    allocate space for and save aerosol name information.
!--------------------------------------------------------------------
      if (nfields /= 0) then
        naerosol = nfields
        allocate (id_aerosol(naerosol))
        allocate (id_aerosol_column(naerosol)) 
        allocate (aerosol_column_names(naerosol))
        do n = 1,naerosol                           
          aerosol_column_names(n) = TRIM(names(n) ) // "_col"
        end do
        do n = 1,naerosol
          id_aerosol(n)    = &
             register_diag_field (mod_name, TRIM(names(n)), axes(1:3), &
                                  Time, TRIM(names(n)),&
                                  'kg/m3', missing_value=missing_value)
          id_aerosol_column(n)    = &
             register_diag_field (mod_name,   &
                      TRIM(aerosol_column_names(n)), axes(1:2), Time, &
                      TRIM(aerosol_column_names(n)), &
                      'kg/m2', missing_value=missing_value)
        end do
        deallocate (aerosol_column_names)

        allocate (extopdep_names(naerosol))
        allocate (extopdep_column_names(naerosol))
        allocate (absopdep_names(naerosol))
        allocate (absopdep_column_names(naerosol))
        allocate (id_extopdep(naerosol, N_DIAG_BANDS))
        allocate (id_extopdep_column(naerosol, N_DIAG_BANDS))
        allocate (id_absopdep(naerosol, N_DIAG_BANDS))
        allocate (id_absopdep_column(naerosol, N_DIAG_BANDS))
     do nl=1,N_DIAG_BANDS
        do n = 1,naerosol                           
          extopdep_names(n) =   &
                TRIM(names(n) ) // "_exopdep" // TRIM(band_suffix(nl))
          extopdep_column_names(n) =   &
             TRIM(names(n) ) // "_exopdep_col" // TRIM(band_suffix(nl))
          absopdep_names(n) =   &
             TRIM(names(n) ) // "_abopdep" // TRIM(band_suffix(nl))
          absopdep_column_names(n) =   &
             TRIM(names(n) ) // "_abopdep_col" // TRIM(band_suffix(nl))
        end do
        do n = 1,naerosol
          id_extopdep(n,nl)    = &
             register_diag_field (mod_name, TRIM(extopdep_names(n)), axes(1:3), &
                                  Time, TRIM(extopdep_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_extopdep_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(extopdep_column_names(n)), axes(1:2), Time, &
                      TRIM(extopdep_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
          id_absopdep(n,nl)    = &
             register_diag_field (mod_name, TRIM(absopdep_names(n)), axes(1:3), &
                                  Time, TRIM(absopdep_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_absopdep_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(absopdep_column_names(n)), axes(1:2), Time, &
                      TRIM(absopdep_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
        end do
      end do
        deallocate (extopdep_names)
        deallocate (extopdep_column_names)
        deallocate (absopdep_names)
        deallocate (absopdep_column_names)
      endif

      if (size(family_names(:)) /= 0) then
        nfamilies = size(family_names(:))
        allocate (id_aerosol_fam(nfamilies))
        allocate (id_aerosol_fam_column(nfamilies)) 
        allocate (aerosol_fam_column_names(naerosol))
        do n=1,nfamilies      
          aerosol_fam_column_names(n) = TRIM(family_names(n) ) // "_col"
        end do
        do n = 1,nfamilies
          id_aerosol_fam(n)    = &
             register_diag_field (mod_name, TRIM(family_names(n)), axes(1:3), &
                                  Time, TRIM(family_names(n)),&
                                  'kg/m3', missing_value=missing_value)
          id_aerosol_fam_column(n)    = &
             register_diag_field (mod_name,   &
                      TRIM(aerosol_fam_column_names(n)), axes(1:2), Time, &
                      TRIM(aerosol_fam_column_names(n)), &
                      'kg/m2', missing_value=missing_value)
        end do
        deallocate (aerosol_fam_column_names)


        allocate (id_extopdep_fam(nfamilies, N_DIAG_BANDS))
        allocate (id_extopdep_fam_column(nfamilies, N_DIAG_BANDS))
        allocate (id_absopdep_fam(nfamilies, N_DIAG_BANDS))
        allocate (id_absopdep_fam_column(nfamilies, N_DIAG_BANDS))
        allocate (extopdep_fam_names(naerosol))
        allocate (extopdep_fam_column_names(naerosol))
        allocate (absopdep_fam_names(naerosol))
        allocate (absopdep_fam_column_names(naerosol))
   do nl=1,N_DIAG_BANDS
        do n=1,nfamilies      
          extopdep_fam_names(n) =   &
           TRIM(family_names(n) ) // "_exopdep" // TRIM(band_suffix(nl))
          extopdep_fam_column_names(n) =   &
       TRIM(family_names(n) ) // "_exopdep_col" // TRIM(band_suffix(nl))
          absopdep_fam_names(n) =   &
          TRIM(family_names(n) ) // "_abopdep" // TRIM(band_suffix(nl))
          absopdep_fam_column_names(n) =  &
       TRIM(family_names(n) ) // "_abopdep_col" // TRIM(band_suffix(nl))
        end do
        do n = 1,nfamilies
          id_extopdep_fam(n,nl)    = &
             register_diag_field (mod_name, TRIM(extopdep_fam_names(n)), axes(1:3), &
                                  Time, TRIM(extopdep_fam_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_extopdep_fam_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(extopdep_fam_column_names(n)), axes(1:2), Time, &
                      TRIM(extopdep_fam_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
          id_absopdep_fam(n,nl)    = &
             register_diag_field (mod_name, TRIM(absopdep_fam_names(n)), axes(1:3), &
                                  Time, TRIM(absopdep_fam_names(n)),&
                                  'dimensionless', missing_value=missing_value)
          id_absopdep_fam_column(n,nl)    = &
             register_diag_field (mod_name,   &
                      TRIM(absopdep_fam_column_names(n)), axes(1:2), Time, &
                      TRIM(absopdep_fam_column_names(n)), &
                      'dimensionless', missing_value=missing_value)
        end do
   end do
        deallocate (extopdep_fam_names)
        deallocate (extopdep_fam_column_names)
        deallocate (absopdep_fam_names)
        deallocate (absopdep_fam_column_names)
      endif
      
      if (Rad_control%volcanic_lw_aerosols_iz) then
        if (Rad_control%volcanic_lw_aerosols) then
          id_lw_extopdep_vlcno_column(1) = &
             register_diag_field (mod_name,   &
                    'lw_b5_extopdep_vlcno_c', axes(1:2), Time, &
                    'lw 900-990 band column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_absopdep_vlcno_column(1)    = &
             register_diag_field (mod_name,   &
                    'lw_b5_absopdep_vlcno_c', axes(1:2), Time, &
                    'lw 900-990 band column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lwext_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_extopdep_vlcno', axes(1:3), Time, &
                    '900-990 band volcanic lw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_xcoeff_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_lwext_vlcno', axes(1:3), Time, &
                    '900-990 band volcanic lw extinction',  &
                      'meter**(-1)  ', missing_value=missing_value)
          id_lwssa_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_lwssa_vlcno', axes(1:3), Time, &
                    '900-990   band volcanic lw scattering albedo', &
                      'dimensionless', missing_value=missing_value)
          id_lwasy_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'bnd5_lwasy_vlcno', axes(1:3), Time, &
                    '900-990 band volcanic lw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_extopdep_vlcno_column(2) = &
             register_diag_field (mod_name,   &
                    'bnd6_extopdep_vlcno_c', axes(1:2), Time, &
                    '990-1070 column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_absopdep_vlcno_column(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_absopdep_vlcno_c', axes(1:2), Time, &
                    '990-1070 column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_lwext_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_extopdep_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_lw_xcoeff_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_lwext_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw extinction',  &
                      'meter**(-1)  ', missing_value=missing_value)
          id_lwssa_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_lwssa_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw scattering albedo',  &
                      'dimensionless', missing_value=missing_value)
          id_lwasy_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'bnd6_lwasy_vlcno', axes(1:3), Time, &
                    '990-1070 volcanic lw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
        endif
      else 
        call error_mesg ('rad_output_file_mod', &
            'Rad_control%volcanic_lw_aerosols not yet defined', FATAL)
      endif
      if (Rad_control%volcanic_sw_aerosols_iz) then
        if (Rad_control%volcanic_sw_aerosols) then
          id_extopdep_vlcno_column(1) = &
             register_diag_field (mod_name,   &
                    'vis_extopdep_vlcno_c', axes(1:2), Time, &
                    'visband column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_absopdep_vlcno_column(1)    = &
             register_diag_field (mod_name,   &
                    'vis_absopdep_vlcno_c', axes(1:2), Time, &
                    'visband column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_swext_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swextopdep_vlcno', axes(1:3), Time, &
                    'visband volcanic sw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_sw_xcoeff_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swext_vlcno', axes(1:3), Time, &
                    'visband volcanic sw extinction',  &
                      'meter**(-1)', missing_value=missing_value)
          id_swssa_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swssa_vlcno', axes(1:3), Time, &
                    'visband volcanic sw scattering albedo',  &
                      'dimensionless', missing_value=missing_value)
          id_swasy_vlcno(1)    = &
             register_diag_field (mod_name,   &
                    'visband_swasy_vlcno', axes(1:3), Time, &
                    'visband volcanic sw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
          id_extopdep_vlcno_column(2) = &
             register_diag_field (mod_name,   &
                    'nir_extopdep_vlcno_c', axes(1:2), Time, &
                    'nirband column volcanic extopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_absopdep_vlcno_column(2)    = &
             register_diag_field (mod_name,   &
                    'nir_absopdep_vlcno_c', axes(1:2), Time, &
                    'nirband column volcanic absopdep',  &
                      'dimensionless', missing_value=missing_value)
          id_swext_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swextopdep_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw extopdep  ',  &
                      'dimensionless', missing_value=missing_value)
          id_sw_xcoeff_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swext_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw extinction',  &
                      'meter**(-1)', missing_value=missing_value)
          id_swssa_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swssa_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw scattering albedo',  &
                      'dimensionless', missing_value=missing_value)
          id_swasy_vlcno(2)    = &
             register_diag_field (mod_name,   &
                    'nirband_swasy_vlcno', axes(1:3), Time, &
                    'nirband volcanic sw asymmetry',  &
                      'dimensionless', missing_value=missing_value)
          id_swheat_vlcno    = &
             register_diag_field (mod_name,   &
                    'sw_heating_vlcno', axes(1:3), Time, &
                    'sw heating due to vlcnic aero',  &
                      'deg K per day', missing_value=missing_value)
        endif
      else 
        call error_mesg ('rad_output_file_mod', &
           'Rad_control%volcanic_sw_aerosols not yet defined', FATAL)
      endif


      id_cmxolw = &
         register_diag_field (mod_name, 'cmxolw', axes(1:3), Time, &
                          'maximum overlap cloud amount', &
                          'percent', missing_value=missing_value)

      id_crndlw = &
         register_diag_field (mod_name, 'crndlw', axes(1:3), Time, &
                          'random overlap cloud amount', &
                          'percent', missing_value=missing_value)

      id_flxnet = &
         register_diag_field (mod_name, 'flxnet', bxes(1:3), Time, &
                          'net longwave radiative flux', &
                          'W/m**2', missing_value=missing_value)

      id_fsw    = &
         register_diag_field (mod_name, 'fsw', bxes(1:3), Time, &
                          'net shortwave radiative flux', &
                          'W/m**2', missing_value=missing_value)

      id_ufsw   = &
         register_diag_field (mod_name, 'ufsw', bxes(1:3), Time, &
                          'upward shortwave radiative flux ', &
                          'W/m**2', missing_value=missing_value)

      id_psj    = &
         register_diag_field (mod_name, 'psj', axes(1:2), Time, &
                          'surface pressure', &
                          'Pa', missing_value=missing_value)

      id_tmpsfc = &
         register_diag_field (mod_name, 'tmpsfc', axes(1:2), Time, &
                          'surface temperature', &
                          'deg_K', missing_value=missing_value)

      id_cvisrfgd_dir = &
         register_diag_field (mod_name, 'cvisrfgd_dir', axes(1:2), Time , &
                         'direct visible surface albedo', &
                        'dimensionless', missing_value=missing_value)

       id_cvisrfgd_dif = &
       register_diag_field (mod_name, 'cvisrfgd_dif', axes(1:2), Time, &
                          'diffuse visible surface albedo', &
                          'dimensionless', missing_value=missing_value)
 
       id_cirrfgd_dir = &
        register_diag_field (mod_name, 'cirrfgd_dir', axes(1:2), Time, &
                       'direct infra-red surface albedo', &
                       'dimensionless', missing_value=missing_value)

       id_cirrfgd_dif = &
         register_diag_field (mod_name, 'cirrfgd_dif', axes(1:2), Time, &
                       'diffuse infra-red surface albedo', &
                      'dimensionless', missing_value=missing_value)

       id_lw_bdyflx(1) = &
         register_diag_field (mod_name, 'olr_800_1200', axes(1:2), Time, &
                       'olr in 800_1200  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(2) = &
         register_diag_field (mod_name, 'olr_900_990', axes(1:2), Time, &
                       'olr in 800_900  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(3) = &
         register_diag_field (mod_name, 'sfc_800_1200', axes(1:2),&
                             Time, 'lw sfc flx in 800_1200  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx(4) = &
         register_diag_field (mod_name, 'sfc_900_990', axes(1:2), &
                              Time, 'lw sfc flx in 900_990 band', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx(1) = &
         register_diag_field (mod_name, 'swup_toa_vis', axes(1:2),  &
                       Time, &
                       'sw up flx in vis band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx(2) = &
         register_diag_field (mod_name, 'swup_toa_1p6', axes(1:2),  &
                       Time, &
                       'sw up flx in 1.6 micron band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx(3) = &
         register_diag_field (mod_name, 'swnt_sfc_vis', axes(1:2),  &
                       Time, &
                       'net sw flx in vis band at sfc', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx(4) = &
         register_diag_field (mod_name, 'swnt_sfc_1p6', axes(1:2),  &
                       Time, &
                       'net sw flx in 1.6 micron band at sfc', &
                      'W/m**2', missing_value=missing_value)

!-------------------------------------------------------------------
!    verify that Rad_control%do_totcld_forcing has been initialized.
!-------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing_iz) then
      else
        call error_mesg ('rad_output_file_mod', &
         ' attempting to use Rad_control%do_totcld_forcing before'//&
                                                ' it is set', FATAL)
      endif

      if (Rad_control%do_totcld_forcing) then
        id_radswpcf = &
           register_diag_field (mod_name, 'radswpcf', axes(1:3), Time, &
                            'temperature forcing from sw w/o clouds', &
                            'deg_K/sec', missing_value=missing_value)

        id_radpcf = &
           register_diag_field (mod_name, 'radpcf', axes(1:3), Time, &
                            'temperature forcing w/o clouds', &
                            'deg_K/sec', missing_value=missing_value)

        id_flxnetcf = &
           register_diag_field (mod_name, 'flxnetcf', bxes(1:3), Time, &
                            'net longwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

        id_fswcf = &
           register_diag_field (mod_name, 'fswcf', bxes(1:3), Time, &
                            'net shortwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

        id_ufswcf   = &
           register_diag_field (mod_name, 'ufswcf', bxes(1:3), Time, &
                            'upward shortwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(1) = &
         register_diag_field (mod_name, 'olr_800_1200_cf', axes(1:2), Time, &
                       'clr sky olr in 800_1200  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(2) = &
         register_diag_field (mod_name, 'olr_900_990_cf', axes(1:2), Time, &
                       'clr sky olr in 800_900  band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(3) = &
         register_diag_field (mod_name, 'sfc_800_1200_cf', axes(1:2),&
                             Time, 'clr sky lw sfc flx in 800_1200 band', &
                      'W/m**2', missing_value=missing_value)

       id_lw_bdyflx_clr(4) = &
         register_diag_field (mod_name, 'sfc_900_990_cf', axes(1:2), &
                              Time, 'clr sky lw sfc flx in 900_990 band', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx_clr(1) = &
         register_diag_field (mod_name, 'swup_toa_vis_cf', axes(1:2),  &
                       Time, &
                       'clr sky sw up flx in vis band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx_clr(2) = &
         register_diag_field (mod_name, 'swup_toa_1p6_cf', axes(1:2),  &
                       Time, &
                       'clr sky sw up flx in 1.6 micron band at toa', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx_clr(3) = &
         register_diag_field (mod_name, 'swnt_sfc_vis_cf', axes(1:2),  &
                       Time, &
                       'clr sky net sw flx in vis band at sfc', &
                      'W/m**2', missing_value=missing_value)

       id_sw_bdyflx_clr(4) = &
         register_diag_field (mod_name, 'swnt_sfc_1p6_cf', axes(1:2),  &
                       Time, &
                       'clr sky net sw flx in 1.6 micron band at sfc', &
                      'W/m**2', missing_value=missing_value)

      endif

!---------------------------------------------------------------------


end subroutine register_fields


!####################################################################





                  end module rad_output_file_mod
