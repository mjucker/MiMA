module flux_exchange_mod
!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman </CONTACT>
! <CONTACT EMAIL="V.Balaji@noaa.gov"> V. Balaji </CONTACT>
! <CONTACT EMAIL="Sergey.Malyshev@noaa.gov"> Sergey Malyshev </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   The flux_exchange module provides interfaces to couple the following component 
!   models: atmosphere, ocean, land, and ice. All interpolation between physically 
!   distinct model grids is handled by the exchange grid (xgrid_mod) with the 
!   interpolated quantities being conserved.
! </OVERVIEW>

! <DESCRIPTION>
!  <PRE>
!  1.This version of flux_exchange_mod allows the definition of physically independent
!    grids for atmosphere, land and sea ice. Ice and ocean must share the same physical
!    grid (though the domain decomposition on parallel systems may be different). 
!    Grid information is input through the grid_spec file (URL). The masked region of the
!    land grid and ice/ocean grid must "tile" each other. The masked region of the ice grid
!    and ocean grid must be identical. 
!
!         ATMOSPHERE  |----|----|----|----|----|----|----|----|
!
!               LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
!
!                ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!
!               OCEAN |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!
!              where  |xxx| = masked grid point
!         
!
!    The atmosphere, land, and ice grids exchange information using the exchange grid xmap_sfc.
!
!    The land and ice grids exchange runoff data using the exchange grid xmap_runoff.
!
!    Transfer of data between the ice bottom and ocean does not require an exchange 
!    grid as the grids are physically identical. The flux routines will automatically
!    detect and redistribute data if their domain decompositions are different.
!
!    To get information from the atmosphere to the ocean it must pass through the 
!    ice model, first by interpolating from the atmospheric grid to the ice grid, 
!    and then transferring from the ice grid to the ocean grid.

!  2.Each component model must have a public defined data type containing specific 
!    boundary fields. A list of these quantities is located in the NOTES of this document. 
!
!  3.The surface flux of sensible heat and surface evaporation can be implicit functions
!    of surface temperature. As a consequence, the parts of the land and sea-ice models 
!    that update the surface temperature must be called on the atmospheric time step 
!
!  4.The surface fluxes of all other tracers and of momentum are assumed to be explicit
!    functions of all surface parameters 
!
!  5.While no explicit reference in made within this module to the implicit treatment 
!    of vertical diffusion in the atmosphere and in the land or sea-ice models, the 
!    module is designed to allow for simultaneous implicit time integration on both 
!    sides of the surface interface. 
!
!  6.Due to #5, the diffusion part of the land and ice models must be called on the 
!    atmospheric time step.
  
!7. Any field passed from one component to another may be "faked" to a
!   constant value, or to data acquired from a file, using the
!   data_override feature of FMS. The fields to override are runtime
!   configurable, using the text file <tt>data_table</tt> for input.
!   See the data_override_mod documentation for more details.
!
!   We DO NOT RECOMMEND exercising the data override capabilities of
!   the FMS coupler until the user has acquired considerable
!   sophistication in running FMS.
!
!   Here is a listing of the override capabilities of the flux_exchange
!   module:
!
!   FROM the atmosphere boundary TO the exchange grid (in sfc_boundary_layer):
!  
!        t_bot, q_bot, z_bot, p_bot, u_bot, v_bot, p_surf, gust
!
!   FROM the ice boundary TO the exchange grid (in sfc_boundary_layer):
!
!        t_surf, rough_mom, rough_heat, rough_moist, albedo, u_surf, v_surf
!     
!   FROM the land boundary TO the exchange grid (in sfc_boundary_layer):
!
!        t_surf, t_ca, q_ca, rough_mom, rough_heat, albedo
!
!   FROM the exchange grid TO land_ice_atmos_boundary (in
!   sfc_boundary_layer):
!
!        t, albedo, land_frac, dt_t, dt_q, u_flux, v_flux, dtaudu, dtaudv,
!        u_star, b_star, rough_mom
!   
!   FROM the atmosphere boundary TO the exchange grid (in
!    flux_down_from_atmos):
!
!        flux_sw, flux_lw, lprec, fprec, coszen, dtmass, delta_t,
!        delta_q, dflux_t, dflux_q
!        
!   FROM the exchange grid TO the land boundary (in
!    flux_down_from_atmos):
!
!    t_flux, q_flux, lw_flux, sw_flux, lprec, fprec, dhdt, dedt, dedq,
!    drdt, drag_q, p_surf
!    
!   FROM the exchange grid TO the ice boundary (in flux_down_from_atmos):
!
!        u_flux, v_flux, t_flux, q_flux, lw_flux, lw_flux_dn, sw_flux,
!        sw_flux_dn, lprec, fprec, dhdt, dedt, drdt, coszen, p 
!
!   FROM the land boundary TO the ice boundary (in flux_land_to_ice):
!
!        runoff, calving
!
!   FROM the ice boundary TO the ocean boundary (in flux_ice_to_ocean):
! 
!        u_flux, v_flux, t_flux, q_flux, salt_flux, lw_flux, sw_flux,
!        lprec, fprec, runoff, calving, p
!        
!   FROM the ocean boundary TO the ice boundary (in flux_ocean_to_ice):
!
!        u, v, t, s, frazil, sea_level
!
!   FROM the ice boundary TO the atmosphere boundary (in flux_up_to_atmos):
!
!        t_surf
!
!   FROM the land boundary TO the atmosphere boundary (in
!    flux_up_to_atmos):
!  
!        t_ca, t_surf, q_ca
!
!  See NOTES below for an explanation of the field names.
!  </PRE>
! </DESCRIPTION>

  use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, &
       mpp_error, stderr, stdlog, FATAL, mpp_set_current_pelist, &
       mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
       CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_ROUTINE
                    
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                             mpp_global_sum, mpp_redistribute, operator(.EQ.)
  use mpp_io_mod,      only: mpp_close

!model_boundary_data_type contains all model fields at the boundary.
!model1_model2_boundary_type contains fields that model2 gets
!from model1, may also include fluxes. These are declared by
!flux_exchange_mod and have private components. All model fields in
!model_boundary_data_type may not be exchanged.
!will support 3 types of flux_exchange:
!REGRID: physically distinct grids, via xgrid
!REDIST: same grid, transfer in index space only
!DIRECT: same grid, same decomp, direct copy
  use atmos_model_mod, only: atmos_data_type, land_ice_atmos_boundary_type
  use ocean_model_mod, only: ocean_data_type, ice_ocean_boundary_type
  use ice_model_mod,   only: ice_data_type, land_ice_boundary_type, &
       ocean_ice_boundary_type, atmos_ice_boundary_type
  use    land_model_mod, only:  land_data_type, atmos_land_boundary_type

  use  surface_flux_mod, only: surface_flux
  use monin_obukhov_mod, only: mo_profile     

  use xgrid_mod, only: xmap_type, setup_xmap, set_frac_area, &
       put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init

  use diag_integral_mod, only:     diag_integral_field_init, &
       sum_diag_integral_field

  use  diag_manager_mod, only: register_diag_field,  &
       register_static_field, send_data, send_tile_averaged_data

  use  time_manager_mod, only: time_type

  use sat_vapor_pres_mod, only: escomp

  use      constants_mod, only: rdgas, rvgas, cp_air, stefan

!Balaji
!utilities stuff into use fms_mod
  use fms_mod, only: clock_flag_default, check_nml_error, error_mesg, &
       open_namelist_file, write_version_number

  use data_override_mod, only: data_override

  implicit none
  include 'netcdf.inc'
private

  public :: flux_exchange_init,   &
     sfc_boundary_layer,   &
     generate_sfc_xgrid,   &
     flux_down_from_atmos, &
     flux_up_to_atmos,     &
     flux_land_to_ice,     &
     flux_ice_to_ocean,    &
     flux_ocean_to_ice

!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: flux_exchange.f90,v 12.0 2005/04/14 15:59:14 fms Exp $'
  character(len=128) :: tag = '$Name: lima $'
!-----------------------------------------------------------------------
!---- exchange grid maps -----

type(xmap_type), save :: xmap_sfc, xmap_runoff

integer         :: n_xgrid_sfc,  n_xgrid_runoff

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=4), parameter :: mod_name = 'flux'

  integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,     &
     id_rough_moist, id_rough_heat, id_rough_mom,    &
     id_land_mask,   id_ice_mask,     &
     id_u_star, id_b_star, id_q_star, id_u_flux, id_v_flux,   &
     id_t_surf, id_t_flux, id_q_flux, id_r_flux,              &
     id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
     id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref,               &
     id_del_h,  id_del_m,  id_del_q,  id_rough_scale,         &
     id_t_ca,   id_q_surf, id_q_atm, id_z_atm, id_p_atm, id_gust, &
     id_t_ref_land, id_rh_ref_land, id_u_ref_land, id_v_ref_land, &
     id_q_ref,  id_q_ref_land

logical :: first_static = .true.
logical :: do_init = .true.
integer :: remap_method = 1

real, parameter :: bound_tol = 1e-7

real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622

!--- namelist interface ------------------------------------------------------
! <NAMELIST NAME="flux_exchange_nml">
!   <DATA NAME="z_ref_heat"  TYPE="real"  DEFAULT="2.0">
!    eference height (meters) for temperature and relative humidity 
!    diagnostics (t_ref,rh_ref,del_h,del_q)
!   </DATA>
!   <DATA NAME="z_ref_mom"  TYPE="real"  DEFAULT="10.0">
!    reference height (meters) for momentum diagnostics (u_ref,v_ref,del_m)
!   </DATA>
!   <DATA NAME="ex_u_star_smooth_bug"  TYPE="logical"  DEFAULT="false">
!    By default, the global exchange grid u_star will not be interpolated from 
!    atmospheric grid, this is different from Jakarta behavior and will
!    change answers. So to perserve Jakarta behavior and reproduce answers
!    explicitly set this namelist variable to .true. in input.nml.
!    Talk to mw, ens for details.
!   </DATA>

  real ::  z_ref_heat =  2.,  &
           z_ref_mom  = 10.
  logical :: ex_u_star_smooth_bug = .false.

namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom, ex_u_star_smooth_bug
! </NAMELIST>

! ---- allocatable module storage --------------------------------------------
real, allocatable, dimension(:) :: &
     ! NOTE: T canopy is only differet from t_surf over vegetated land
     ex_t_surf,    &   ! surface temperature for radiation calc, degK
     ex_t_ca,      &   ! near-surface (canopy) air temperature, degK
     ex_q_surf,    &   ! near-surface (canopy) air specific humidity, kg/kg
     ex_p_surf,    &   ! surface pressure

     ex_flux_t,    &   ! sens heat flux
     ex_flux_q,    &   ! water vapor flux
     ex_flux_lw,   &   ! longwave radiation flux

     ex_dhdt_surf, &   ! d(sens.heat.flux)/d(T canopy)
     ex_dedt_surf, &   ! d(water.vap.flux)/d(T canopy)
     ex_dedq_surf, &   ! d(water.vap.flux)/d(q canopy)
     ex_drdt_surf, &   ! d(LW flux)/d(T surf)
     ex_dhdt_atm,  &   ! d(sens.heat.flux)/d(T atm)
     ex_dedq_atm,  &   ! d(water.vap.flux)/d(q atm)
     ex_flux_u,    &   ! u stress on atmosphere
     ex_flux_v,    &   ! v stress on atmosphere
     ex_dtaudu_atm,&   ! d(stress)/d(u)
     ex_dtaudv_atm,&   ! d(stress)/d(v)
     ex_albedo_fix,&
     ex_albedo_vis_dir_fix,&
     ex_albedo_nir_dir_fix,&
     ex_albedo_vis_dif_fix,&
     ex_albedo_nir_dif_fix,&
     ex_old_albedo,&   ! old value of albedo for downward flux calculations
     ex_drag_q,    &   ! q drag.coeff.
     ex_cd_t,      &
     ex_cd_m,      &
     ex_b_star,    &
     ex_u_star,    &
     ex_wind,      &
     ex_z_atm

logical, allocatable, dimension(:) :: &
     ex_avail,     &   ! true where data on exchange grid are available
     ex_land           ! true if exchange grid cell is over land
real, allocatable, dimension(:) :: &
     ex_e_t_n,      &
     ex_f_t_delt_n, &
     ex_e_q_n,      &
     ex_f_q_delt_n

integer :: ni_atm, nj_atm ! to do atmos diagnostic from flux_ocean_to_ice
real, dimension(3) :: ccc ! for conservation checks
!Balaji, sets boundary_type%xtype
!  REGRID: grids are physically different, pass via exchange grid
!  REDIST: same physical grid, different decomposition, must move data around
!  DIRECT: same physical grid, same domain decomposition, can directly copy data
integer, parameter :: REGRID=1, REDIST=2, DIRECT=3
!Balaji: clocks moved into flux_exchange
  integer :: cplClock, sfcClock, fluxAtmDnClock, fluxLandIceClock, &
             fluxIceOceanClock, fluxOceanIceClock, regenClock, fluxAtmUpClock, &
             cplOcnClock
contains

!#######################################################################
! <SUBROUTINE NAME="flux_exchange_init">
!  <OVERVIEW>
!   Initialization routine.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initializes the interpolation routines,diagnostics and boundary data
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
!		atmos_ice_boundary, land_ice_atmos_boundary, &
!		land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </IN>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Ocean" TYPE="ocean_data_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <INOUT NAME="atmos_ice_boundary" TYPE="atmos_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to ice.
!  </INOUT>
!  <INOUT NAME="land_ice_atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.
!  </INOUT>
!  <INOUT NAME="land_ice_boundary" TYPE="land_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from land to ice.
!  </INOUT>
!  <INOUT NAME="ice_ocean_boundary" TYPE="ice_ocean_boundary_type">
!  A derived data type to specify properties and fluxes passed from ice to ocean.
!  </INOUT>
!  <INOUT NAME="ocean_ice_boundary" TYPE="ocean_ice_boundary_type">
!  A derived data type to specify properties and fluxes passed from ocean to ice.
!  </INOUT>
!
subroutine flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
       atmos_ice_boundary, land_ice_atmos_boundary, &
       land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )

  type(time_type),                   intent(in)  :: Time
  type(atmos_data_type),             intent(in)  :: Atm
  type(land_data_type),              intent(in)  :: Land
  type(ice_data_type),               intent(in)  :: Ice
  type(ocean_data_type),             intent(in)  :: Ocean
! All intent(OUT) derived types with pointer components must be 
! COMPLETELY allocated here and in subroutines called from here;
! NO pointer components should have been allocated before entry if the
! derived type has intent(OUT) otherwise they may be lost.
  type(atmos_ice_boundary_type),     intent(inout) :: atmos_ice_boundary
  type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary
  type(land_ice_boundary_type),      intent(inout) :: land_ice_boundary
  type(ice_ocean_boundary_type),     intent(inout) :: ice_ocean_boundary
  type(ocean_ice_boundary_type),     intent(inout) :: ocean_ice_boundary
  
  integer :: unit, ierr, io, iref, i
  integer :: rcode, ncid, varid, dims(4), start(4), nread(4)
  integer :: nlon, nlat
  integer :: ndims, n_ocean_areas, n_land_areas
  real, dimension(:), allocatable :: atmlonb, atmlatb
  integer :: is, ie, js, je, kd
  integer, dimension(:), allocatable :: pelist  

!-----------------------------------------------------------------------
!----- read namelist -------

    unit = open_namelist_file()
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=flux_exchange_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'flux_exchange_nml')
    enddo
10  call mpp_close(unit)

!----- write namelist to logfile -----
    call write_version_number (version, tag)
    if( mpp_pe() == mpp_root_pe() )write( stdlog(), nml=flux_exchange_nml )


!--------- read gridspec file ------------------
!only atmos pelists needs to do it here, ocean model will do it elsewhere

    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        rcode = nf_open('INPUT/grid_spec.nc',0,ncid)
!   <ERROR MSG="cannot open INPUT/grid_spec.nc" STATUS="FATAL">
!      Can not open the file INPUT/grid_spec.nc. The possible reason is 
!      that file grid_spec.nc is not in the directory INPUT.
!   </ERROR>
        if (rcode/=0) call error_mesg ('flux_exchange_mod', &
             'cannot open INPUT/grid_spec.nc', FATAL)

!
! check atmosphere and grid_spec.nc have same atmosphere lat/lon boundaries
!
        rcode = nf_inq_varid(ncid, 'AREA_ATM', varid)
!   <ERROR MSG="cannot find AREA_ATM on INPUT/grid_spec.nc" STATUS="FATAL">
!      File INPUT/grid_spec.nc does not contain the field AREA_ATM.
!   </ERROR>
        if (rcode/=0) call error_mesg ('flux_exchange_mod', &
             'cannot find AREA_ATM on INPUT/grid_spec.nc', &
             FATAL)
        rcode = nf_inq_vardimid(ncid, varid, dims)
        rcode = nf_inq_dimlen(ncid, dims(1), nlon)
        rcode = nf_inq_dimlen(ncid, dims(2), nlat)
        if (nlon+1/=size(Atm%glon_bnd(:)).or.nlat+1/=size(Atm%glat_bnd(:))) then
            if (mpp_pe()==mpp_root_pe()) then
                print *, 'grid_spec.nc has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                     'atmosphere has', size(Atm%glon_bnd(:))-1, 'longitudes,', &
                     size(Atm%glat_bnd(:))-1, 'latitudes (see xba.dat and yba.dat)'
            end if
!   <ERROR MSG="grid_spec.nc incompatible with atmosphere resolution" STATUS="FATAL">
!      The atmosphere grid size from file grid_spec.nc is not compatible with the atmosphere 
!      resolution from atmosphere model.
!   </ERROR>
            call error_mesg ('flux_exchange_mod',  &
                 'grid_spec.nc incompatible with atmosphere resolution', FATAL)
        end if
        allocate( atmlonb(size(Atm%glon_bnd(:))) )
        allocate( atmlatb(size(Atm%glat_bnd(:))) )
        rcode = nf_inq_varid(ncid, 'xba', varid)
        start = 1; nread = 1; nread(1) = nlon+1;
        rcode = nf_get_vara_double(ncid, varid, start, nread, atmlonb)
        rcode = nf_inq_varid(ncid, 'yba', varid)
        start = 1; nread = 1; nread(1) = nlat+1;
        rcode = nf_get_vara_double(ncid, varid, start, nread, atmlatb)
        if (maxval(abs(atmlonb-Atm%glon_bnd*45/atan(1.0)))>bound_tol) then
            if (mpp_pe() == mpp_root_pe()) then
                print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY'
                do i=1,size(atmlonb(:))
                   print *,atmlonb(i),Atm%glon_bnd(i)*45/atan(1.0)
                end do
            end if
!   <ERROR MSG="grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)" STATUS="FATAL">
!      longitude from file grid_spec.nc ( from field xba ) is different from the longitude from atmosphere model.
!   </ERROR>
            call error_mesg ('flux_exchange_mod', &
                 'grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)'&
                 ,FATAL)
        end if
        if (maxval(abs(atmlatb-Atm%glat_bnd*45/atan(1.0)))>bound_tol) then
            if (mpp_pe() == mpp_root_pe()) then
                print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY'
                do i=1,size(atmlatb(:))
                   print *,atmlatb(i),Atm%glat_bnd(i)*45/atan(1.0)
                end do
            end if
!   <ERROR MSG="grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)" STATUS="FATAL">
!      longitude from file grid_spec.nc ( from field yba ) is different from the longitude from atmosphere model.
!   </ERROR>
            call error_mesg ('flux_exchange_mod', &
                 'grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)'&
                 , FATAL)
        end if

        call xgrid_init(remap_method)

        call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /),   &
             (/ Atm%Domain, Ice%Domain, Land%Domain /),        &
             "INPUT/grid_spec.nc")
        call generate_sfc_xgrid( Land, Ice )

        call setup_xmap(xmap_runoff, (/ 'LND', 'OCN' /),       &
             (/ Land%Domain, Ice%Domain /),                    &
             "INPUT/grid_spec.nc"             )
        n_xgrid_runoff = max(xgrid_count(xmap_runoff),1)

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!----- initialize quantities for global integral package -----

!! call diag_integral_field_init ('prec', 'f6.3')
        call diag_integral_field_init ('evap', 'f6.3')

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----
!----- all fields will be output on the atmospheric grid -----

        call diag_field_init ( Time, Atm%axes(1:2), Land%axes )
        ni_atm = size(Atm%lon_bnd(:))-1 ! to dimension "diag_atm"
        nj_atm = size(Atm%lat_bnd(:))-1 ! in flux_ocean_to_ice

!Balaji
        
!allocate atmos_ice_boundary
        call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
        kd = size(Ice%ice_mask,3)
        allocate( atmos_ice_boundary%u_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%v_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%u_star(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%t_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%q_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lw_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_dir(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis_dir(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_dif(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis_dif(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%fprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dhdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dedt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%drdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%coszen(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%p(is:ie,js:je,kd) )
! initialize boundary values for override experiments (mjh)
        atmos_ice_boundary%u_flux=0.0
        atmos_ice_boundary%v_flux=0.0
        atmos_ice_boundary%u_star=0.0
        atmos_ice_boundary%t_flux=0.0
        atmos_ice_boundary%q_flux=0.0
        atmos_ice_boundary%lw_flux=0.0
        atmos_ice_boundary%sw_flux=0.0
        atmos_ice_boundary%sw_flux_vis=0.0
        atmos_ice_boundary%sw_flux_dir=0.0
        atmos_ice_boundary%sw_flux_vis_dir=0.0
        atmos_ice_boundary%sw_flux_dif=0.0
        atmos_ice_boundary%sw_flux_vis_dif=0.0
        atmos_ice_boundary%lprec=0.0
        atmos_ice_boundary%fprec=0.0
        atmos_ice_boundary%dhdt=0.0
        atmos_ice_boundary%dedt=0.0
        atmos_ice_boundary%drdt=0.0
        atmos_ice_boundary%coszen=0.0
        atmos_ice_boundary%p=0.0
        
!allocate land_ice_boundary
        allocate( land_ice_boundary%runoff(is:ie,js:je) )
        allocate( land_ice_boundary%calving(is:ie,js:je) )
! initialize values for override experiments (mjh)
        land_ice_boundary%runoff=0.0
        land_ice_boundary%calving=0.0
!allocate land_ice_atmos_boundary
        call mpp_get_compute_domain( Atm%domain, is, ie, js, je )
        allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%land_frac(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_q(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%u_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%v_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudu(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudv(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%u_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%b_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%q_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%rough_mom(is:ie,js:je) )
! initialize boundary values for override experiments (mjh)
        land_ice_atmos_boundary%t=273.0
        land_ice_atmos_boundary%albedo=0.0
        land_ice_atmos_boundary%albedo_vis_dir=0.0
        land_ice_atmos_boundary%albedo_nir_dir=0.0
        land_ice_atmos_boundary%albedo_vis_dif=0.0
        land_ice_atmos_boundary%albedo_nir_dif=0.0
        land_ice_atmos_boundary%land_frac=0.0
        land_ice_atmos_boundary%dt_t=0.0
        land_ice_atmos_boundary%dt_q=0.0
        land_ice_atmos_boundary%u_flux=0.0
        land_ice_atmos_boundary%v_flux=0.0
        land_ice_atmos_boundary%dtaudu=0.0
        land_ice_atmos_boundary%dtaudv=0.0
        land_ice_atmos_boundary%u_star=0.0
        land_ice_atmos_boundary%b_star=0.0
        land_ice_atmos_boundary%q_star=0.0
        land_ice_atmos_boundary%rough_mom=0.01
!Balaji: clocks on atm%pe only        
    cplClock = mpp_clock_id( 'Land-ice-atm coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    sfcClock = mpp_clock_id( 'SFC boundary layer', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    fluxAtmDnClock = mpp_clock_id( 'Flux DN from atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxLandIceClock = mpp_clock_id( 'Flux land to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    regenClock = mpp_clock_id( 'XGrid generation', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxAtmUpClock = mpp_clock_id( 'Flux UP to atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    end if

    allocate(pelist(size(Atm%pelist(:))+size(Ocean%pelist(:))))
    pelist(1:size(Atm%pelist(:))) = Atm%pelist(:)
    pelist(size(Atm%pelist(:))+1:size(pelist(:))) = Ocean%pelist(:)
    call mpp_set_current_pelist(pelist)
    deallocate(pelist)

!ocean_ice_boundary and ice_ocean_boundary must be done on all PES
!domain boundaries will assure no space is allocated on non-relevant PEs.
    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
!allocate ocean_ice_boundary
    allocate( ocean_ice_boundary%u(is:ie,js:je) )
    allocate( ocean_ice_boundary%v(is:ie,js:je) )
    allocate( ocean_ice_boundary%t(is:ie,js:je) )
    allocate( ocean_ice_boundary%s(is:ie,js:je) )
!frazil and sea_level are optional, if not present they should be nullified
    allocate( ocean_ice_boundary%frazil(is:ie,js:je) )
    allocate( ocean_ice_boundary%sea_level(is:ie,js:je) )
! initialize boundary fields for override experiments (mjh)
    ocean_ice_boundary%u=0.0
    ocean_ice_boundary%v=0.0
    ocean_ice_boundary%t=273.0
    ocean_ice_boundary%s=0.0
    ocean_ice_boundary%frazil=0.0
    ocean_ice_boundary%sea_level=0.0
!allocate ice_ocean_boundary
    call mpp_get_compute_domain( Ocean%domain, is, ie, js, je )
!ML ocean only requires t, q, lw, sw, fprec, calving
!AMIP ocean needs no input fields
!choice of fields will eventually be done at runtime
!via field_manager
    allocate( ice_ocean_boundary%u_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%v_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%t_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%q_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%salt_flux(is:ie,js:je) )
    allocate( ice_ocean_boundary%lw_flux  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_vis  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_vis_dir  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_vis_dif  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_dir  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_dif  (is:ie,js:je) )
    allocate( ice_ocean_boundary%lprec    (is:ie,js:je) )
    allocate( ice_ocean_boundary%fprec    (is:ie,js:je) )
    allocate( ice_ocean_boundary%runoff   (is:ie,js:je) )
    allocate( ice_ocean_boundary%calving  (is:ie,js:je) )
    allocate( ice_ocean_boundary%p        (is:ie,js:je) )
!pjp Why are the above not initialized to zero?
! initialize boundary values for override experiments
    ocean_ice_boundary%xtype = REDIST
    if( Ocean%domain.EQ.Ice%domain )ocean_ice_boundary%xtype = DIRECT
    ice_ocean_boundary%xtype = ocean_ice_boundary%xtype
!BalaCji
    cplOcnClock = mpp_clock_id( 'Ice-ocean coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    fluxIceOceanClock = mpp_clock_id( 'Flux ice to ocean', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxOceanIceClock = mpp_clock_id( 'Flux ocean to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )
!---- done ----
    do_init = .false.

  end subroutine flux_exchange_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="sfc_boundary_layer">
!  <OVERVIEW>
!   Computes explicit fluxes as well as derivatives that will be used to compute an implicit flux correction. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!  The following quantities in the land_ice_atmos_boundary_type are computed:
!
!     
!         t_surf_atm = surface temperature (used for radiation)    (K)
!         albedo_atm = surface albedo      (used for radiation)    (nondimensional)
!      rough_mom_atm = surface roughness for momentum (m)
!      land_frac_atm = fractional area of land beneath an atmospheric
!                      grid box 
!         dtaudu_atm, dtaudv_atm = derivatives of wind stress w.r.t. the
!                      lowest level wind speed  (Pa/(m/s))
!         flux_u_atm = zonal wind stress  (Pa)
!         flux_v_atm = meridional wind stress (Pa)
!         u_star_atm = friction velocity (m/s)
!         b_star_atm = buoyancy scale    (m2/s)
!
!         (u_star and b_star are defined so that u_star**2 = magnitude
!           of surface stress divided by density of air at the surface, 
!           and u_star*b_star = buoyancy flux at the surface)
!
!   </PRE>
!  </DESCRIPTION>

!  <TEMPLATE>
!   call sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" dt" TYPE="real">
!   time step. 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </INOUT>
!  <INOUT NAME="Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.   
!  </INOUT>
!
subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Boundary )

  real,                  intent(in)  :: dt
  type(time_type),       intent(in)  :: Time
  type(atmos_data_type), intent(inout)  :: Atm
  type(land_data_type),  intent(inout)  :: Land
  type(ice_data_type),   intent(inout)  :: Ice
  type(land_ice_atmos_boundary_type), intent(inout) :: Boundary

  ! ---- local vars ----------------------------------------------------------
  real, dimension(n_xgrid_sfc) :: &
       ex_albedo,     &
       ex_albedo_vis_dir,     &
       ex_albedo_nir_dir,     &
       ex_albedo_vis_dif,     &
       ex_albedo_nir_dif,     &
       ex_land_frac,  &
       ex_t_atm,      & 
       ex_q_atm,      &
       ex_p_atm,      &
       ex_u_atm, ex_v_atm,    &
       ex_gust,       &
       ex_t_surf4,    &
       ex_u_surf, ex_v_surf,  &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, &
       ex_rough_scale,&
       ex_q_star,     &
       ex_cd_q,       &
       ex_ref,        &
       ex_t_ref,      &
       ex_qs_ref,     &
       ex_del_m,      &
       ex_del_h,      &
       ex_del_q,      &
       ex_seawater

  real, dimension(size(Boundary%t,1),size(Boundary%t,2)) :: diag_atm
  real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2), size(Land%t_ca,3)) :: diag_land
  real, dimension(size(Ice%t_surf,1),size(Ice%t_surf,2),size(Ice%t_surf,3)) :: sea
  real    :: zrefm, zrefh
  logical :: used


  ! [1] check that the module was initialized
!   <ERROR MSG="must call flux_exchange_init first " STATUS="FATAL">
!      flux_exchange_init has not been called before calling sfc_boundary_layer.
!   </ERROR>
  if (do_init) call error_mesg ('flux_exchange_mod',  &
       'must call flux_exchange_init first', FATAL)
!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(sfcClock)
  ! [2] allocate storage for variables that are also used in flux_up_to_atmos
  allocate ( &
       ex_t_surf   (n_xgrid_sfc),  &
       ex_p_surf   (n_xgrid_sfc),  &
       ex_t_ca     (n_xgrid_sfc),  &
       ex_q_surf   (n_xgrid_sfc),  &
       ex_dhdt_surf(n_xgrid_sfc),  &
       ex_dedt_surf(n_xgrid_sfc),  &
       ex_dedq_surf(n_xgrid_sfc),  &
       ex_drdt_surf(n_xgrid_sfc),  &
       ex_dhdt_atm (n_xgrid_sfc),  &
       ex_dedq_atm (n_xgrid_sfc),  &
       ex_flux_t   (n_xgrid_sfc),  &
       ex_flux_q   (n_xgrid_sfc),  &
       ex_flux_lw  (n_xgrid_sfc),  &
       ex_drag_q   (n_xgrid_sfc),  &
       ex_avail    (n_xgrid_sfc),  &
       ex_f_t_delt_n(n_xgrid_sfc), &
       ex_f_q_delt_n(n_xgrid_sfc), &

! MOD these were moved from local ! so they can be passed to flux down
       ex_flux_u(n_xgrid_sfc),    &
       ex_flux_v(n_xgrid_sfc),    &
       ex_dtaudu_atm(n_xgrid_sfc),&
       ex_dtaudv_atm(n_xgrid_sfc),&

! values added for LM3
       ex_cd_t     (n_xgrid_sfc),  &
       ex_cd_m     (n_xgrid_sfc),  &
       ex_b_star   (n_xgrid_sfc),  &
       ex_u_star   (n_xgrid_sfc),  &
       ex_wind     (n_xgrid_sfc),  &
       ex_z_atm    (n_xgrid_sfc),  &

       ex_e_t_n    (n_xgrid_sfc),  &
       ex_e_q_n    (n_xgrid_sfc),  &
       ex_land     (n_xgrid_sfc)   )

  ! [3] initialize some values on exchange grid: this is actually a safeguard
  ! against using undefined values
  ex_t_surf   = 200.
  ex_u_surf   =   0.
  ex_v_surf   =   0.

  !---- do not use if relax time /= 0 ----
  ex_cd_t = 0.0
  ex_cd_m = 0.0
  ex_cd_q = 0.0
!-----------------------------------------------------------------------
!Balaji: data_override stuff moved from coupler_main
  call data_override ('ATM', 't_bot',  Atm%t_bot , Time)
  call data_override ('ATM', 'q_bot',  Atm%q_bot , Time)
  call data_override ('ATM', 'z_bot',  Atm%z_bot , Time)
  call data_override ('ATM', 'p_bot',  Atm%p_bot , Time)
  call data_override ('ATM', 'u_bot',  Atm%u_bot , Time)
  call data_override ('ATM', 'v_bot',  Atm%v_bot , Time)
  call data_override ('ATM', 'p_surf', Atm%p_surf, Time)
  call data_override ('ATM', 'gust',   Atm%gust,   Time)
  call data_override ('ICE', 't_surf',     Ice%t_surf,      Time)
  call data_override ('ICE', 'rough_mom',  Ice%rough_mom,   Time)
  call data_override ('ICE', 'rough_heat', Ice%rough_heat,  Time)
  call data_override ('ICE', 'rough_moist',Ice%rough_moist, Time)
  call data_override ('ICE', 'albedo',     Ice%albedo,      Time)
  call data_override ('ICE', 'albedo_vis_dir', Ice%albedo_vis_dir, Time)
  call data_override ('ICE', 'albedo_nir_dir', Ice%albedo_nir_dir, Time)
  call data_override ('ICE', 'albedo_vis_dif', Ice%albedo_vis_dif, Time)
  call data_override ('ICE', 'albedo_nir_dif', Ice%albedo_nir_dif, Time)
  call data_override ('ICE', 'u_surf',     Ice%u_surf,      Time)
  call data_override ('ICE', 'v_surf',     Ice%v_surf,      Time)
  call data_override ('LND', 't_surf',     Land%t_surf,     Time)
  call data_override ('LND', 't_ca',       Land%t_ca,       Time)
  call data_override ('LND', 'q_ca',       Land%q_ca,       Time)
  call data_override ('LND', 'rough_mom',  Land%rough_mom,  Time)
  call data_override ('LND', 'rough_heat', Land%rough_heat, Time)
  call data_override ('LND', 'albedo', Land%albedo,     Time)
  call data_override ('LND', 'albedo_vis_dir', Land%albedo_vis_dir,Time)
  call data_override ('LND', 'albedo_nir_dir', Land%albedo_nir_dir,Time)
  call data_override ('LND', 'albedo_vis_dif', Land%albedo_vis_dif,Time)
  call data_override ('LND', 'albedo_nir_dif', Land%albedo_nir_dif,Time)

!---- put atmosphere quantities onto exchange grid ----

  ! [4] put all the qantities we need onto exchange grid
  ! [4.1] put atmosphere quantities onto exchange grid
  call put_to_xgrid (Atm%t_bot , 'ATM', ex_t_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%q_bot , 'ATM', ex_q_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%z_bot , 'ATM', ex_z_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%p_bot , 'ATM', ex_p_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%p_surf, 'ATM', ex_p_surf, xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%gust,   'ATM', ex_gust,   xmap_sfc, remap_method=remap_method)

  ! slm, Mar 20 2002: changed order in whith the data transferred from ice and land 
  ! grids, to fill t_ca first with t_surf over ocean and then with t_ca from 
  ! land, where it is different from t_surf. It is mostly to simplify 
  ! diagnostic, since surface_flux calculations distinguish between land and 
  ! not-land anyway.

  ! [4.2] put ice quantities onto exchange grid
  ! (assume that ocean quantites are stored in no ice partition)
  ! (note: ex_avail is true at ice and ocean points)
  call put_to_xgrid (Ice%t_surf,      'OCN', ex_t_surf,      xmap_sfc)
  call put_to_xgrid (Ice%rough_mom,   'OCN', ex_rough_mom,   xmap_sfc)
  call put_to_xgrid (Ice%rough_heat,  'OCN', ex_rough_heat,  xmap_sfc)
  call put_to_xgrid (Ice%rough_moist, 'OCN', ex_rough_moist, xmap_sfc)
  call put_to_xgrid (Ice%albedo,      'OCN', ex_albedo,      xmap_sfc)
  call put_to_xgrid (Ice%albedo_vis_dir, 'OCN', ex_albedo_vis_dir, xmap_sfc)
  call put_to_xgrid (Ice%albedo_nir_dir, 'OCN', ex_albedo_nir_dir, xmap_sfc)
  call put_to_xgrid (Ice%albedo_vis_dif, 'OCN', ex_albedo_vis_dif, xmap_sfc)
  call put_to_xgrid (Ice%albedo_nir_dif, 'OCN', ex_albedo_nir_dif, xmap_sfc)
  call put_to_xgrid (Ice%u_surf,      'OCN', ex_u_surf,      xmap_sfc)
  call put_to_xgrid (Ice%v_surf,      'OCN', ex_v_surf,      xmap_sfc)
  sea = 0.0; sea(:,:,1) = 1.0;
  ex_seawater = 0.0
  call put_to_xgrid (sea,             'OCN', ex_seawater,    xmap_sfc)
  ex_t_ca = ex_t_surf ! slm, Mar 20 2002 to define values over the ocean

  ! [4.3] put land quantities onto exchange grid ----
  call some(xmap_sfc, ex_land, 'LND')

  call put_to_xgrid (Land%t_surf,     'LND', ex_t_surf,      xmap_sfc)
  call put_to_xgrid (Land%t_ca,       'LND', ex_t_ca,        xmap_sfc)
  call put_to_xgrid (Land%q_ca,       'LND', ex_q_surf,      xmap_sfc)
  call put_to_xgrid (Land%rough_mom,  'LND', ex_rough_mom,   xmap_sfc)
  call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_heat,  xmap_sfc)
  call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_moist, xmap_sfc)
  call put_to_xgrid (Land%albedo,     'LND', ex_albedo,      xmap_sfc)
  call put_to_xgrid (Land%albedo_vis_dir,     'LND', ex_albedo_vis_dir,   xmap_sfc)
  call put_to_xgrid (Land%albedo_nir_dir,     'LND', ex_albedo_nir_dir,   xmap_sfc)
  call put_to_xgrid (Land%albedo_vis_dif,     'LND', ex_albedo_vis_dif,   xmap_sfc)
  call put_to_xgrid (Land%albedo_nir_dif,     'LND', ex_albedo_nir_dif,   xmap_sfc)
  ex_rough_scale = ex_rough_mom
  call put_to_xgrid(Land%rough_scale, 'LND', ex_rough_scale, xmap_sfc)

  ex_land_frac = 0.0
  call put_logical_to_real (Land%mask,    'LND', ex_land_frac, xmap_sfc)

  ! [5] compute explicit fluxes and tendencies at all available points ---
  call some(xmap_sfc, ex_avail)
  call surface_flux (&
       ex_t_atm, ex_q_atm,  ex_u_atm, ex_v_atm,  ex_p_atm,  ex_z_atm,  &
       ex_p_surf,ex_t_surf, ex_t_ca,  ex_q_surf,                       &
       ex_u_surf, ex_v_surf,                                           &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, ex_rough_scale,    &
       ex_gust,                                                        &
       ex_flux_t, ex_flux_q, ex_flux_lw, ex_flux_u, ex_flux_v,         &
       ex_cd_m,   ex_cd_t, ex_cd_q,                                    &
       ex_wind,   ex_u_star, ex_b_star, ex_q_star,                     &
       ex_dhdt_surf, ex_dedt_surf, ex_dedq_surf,  ex_drdt_surf,        &
       ex_dhdt_atm,  ex_dedq_atm,   ex_dtaudu_atm,  ex_dtaudv_atm,     &
       dt,                                                             &
       ex_land, ex_seawater .gt. 0,  ex_avail                          )

  where (ex_avail) ex_drag_q = ex_wind*ex_cd_q
  ! [6] get mean quantities on atmosphere grid
  ! [6.1] compute t surf for radiation
  ex_t_surf4 = ex_t_surf ** 4

  ! [6.2] put relevant quantities onto atmospheric boundary
  call get_from_xgrid (Boundary%t,         'ATM', ex_t_surf4  ,  xmap_sfc)
  call get_from_xgrid (Boundary%albedo,    'ATM', ex_albedo   ,  xmap_sfc)
  call get_from_xgrid (Boundary%albedo_vis_dir,    'ATM',   &
                       ex_albedo_vis_dir   ,  xmap_sfc)
  call get_from_xgrid (Boundary%albedo_nir_dir,    'ATM',   &
                       ex_albedo_nir_dir   ,  xmap_sfc)
  call get_from_xgrid (Boundary%albedo_vis_dif,    'ATM',   &
                       ex_albedo_vis_dif   ,  xmap_sfc)
  call get_from_xgrid (Boundary%albedo_nir_dif,    'ATM',   &
                       ex_albedo_nir_dif   ,  xmap_sfc)
  call get_from_xgrid (Boundary%rough_mom, 'ATM', ex_rough_mom,  xmap_sfc)
  call get_from_xgrid (Boundary%land_frac, 'ATM', ex_land_frac,  xmap_sfc)

  call get_from_xgrid (Boundary%u_flux,    'ATM', ex_flux_u,     xmap_sfc)
  call get_from_xgrid (Boundary%v_flux,    'ATM', ex_flux_v,     xmap_sfc)
  call get_from_xgrid (Boundary%dtaudu,    'ATM', ex_dtaudu_atm, xmap_sfc)
  call get_from_xgrid (Boundary%dtaudv,    'ATM', ex_dtaudv_atm, xmap_sfc)
  call get_from_xgrid (Boundary%u_star,    'ATM', ex_u_star    , xmap_sfc)
  call get_from_xgrid (Boundary%b_star,    'ATM', ex_b_star    , xmap_sfc)
  call get_from_xgrid (Boundary%q_star,    'ATM', ex_q_star    , xmap_sfc)

  Boundary%t = Boundary%t ** 0.25
!Balaji: data_override calls moved here from coupler_main
  call data_override('ATM', 't',         Boundary%t,         Time)
  call data_override('ATM', 'albedo',    Boundary%albedo,    Time)
  call data_override('ATM', 'albedo_vis_dir',    Boundary%albedo_vis_dir,    Time)
  call data_override('ATM', 'albedo_nir_dir',    Boundary%albedo_nir_dir,    Time)
  call data_override('ATM', 'albedo_vis_dif',    Boundary%albedo_vis_dif,    Time)
  call data_override('ATM', 'albedo_nir_dif',    Boundary%albedo_nir_dif,    Time)
  call data_override('ATM', 'land_frac', Boundary%land_frac, Time)
  call data_override('ATM', 'dt_t',      Boundary%dt_t,      Time)
  call data_override('ATM', 'dt_q',      Boundary%dt_q,      Time)
  call data_override('ATM', 'u_flux',    Boundary%u_flux,    Time)
  call data_override('ATM', 'v_flux',    Boundary%v_flux,    Time)
  call data_override('ATM', 'dtaudu',    Boundary%dtaudu,    Time)
  call data_override('ATM', 'dtaudv',    Boundary%dtaudv,    Time)
  call data_override('ATM', 'u_star',    Boundary%u_star,    Time)
  call data_override('ATM', 'b_star',    Boundary%b_star,    Time)
! call data_override('ATM', 'q_star',    Boundary%q_star,    Time)
  call data_override('ATM', 'rough_mom', Boundary%rough_mom, Time)

  ! [6.3] save atmos albedo fix and old albedo (for downward SW flux calculations)
  ! on exchange grid
  ! allocate ( ex_old_albedo(n_xgrid_sfc)  )
  ! ex_old_albedo = ex_albedo
  
!!  STILL NEEDED   ????
!! IS THIS CORRECT ??
  allocate ( ex_albedo_fix(n_xgrid_sfc) )
  call put_to_xgrid (Boundary%albedo, 'ATM',  ex_albedo_fix, xmap_sfc)
  ex_albedo_fix = (1.0-ex_albedo) / (1.0-ex_albedo_fix)

  allocate ( ex_albedo_vis_dir_fix(n_xgrid_sfc) )
  call put_to_xgrid (Boundary%albedo_vis_dir, 'ATM',  &
           ex_albedo_vis_dir_fix, xmap_sfc)
  ex_albedo_vis_dir_fix = (1.0-ex_albedo_vis_dir) /  &
(1.0-ex_albedo_vis_dir_fix)
  allocate ( ex_albedo_nir_dir_fix(n_xgrid_sfc) )
  call put_to_xgrid (Boundary%albedo_nir_dir, 'ATM', &
 ex_albedo_nir_dir_fix, xmap_sfc)
  ex_albedo_nir_dir_fix = (1.0-ex_albedo_nir_dir) /  &
(1.0-ex_albedo_nir_dir_fix)
  allocate ( ex_albedo_vis_dif_fix(n_xgrid_sfc) )
  call put_to_xgrid (Boundary%albedo_vis_dif, 'ATM',   &
        ex_albedo_vis_dif_fix, xmap_sfc)
   ex_albedo_vis_dif_fix = (1.0-ex_albedo_vis_dif) /   &
       (1.0-ex_albedo_vis_dif_fix)
  allocate ( ex_albedo_nir_dif_fix(n_xgrid_sfc) )
  call put_to_xgrid (Boundary%albedo_nir_dif, 'ATM',  &
 ex_albedo_nir_dif_fix, xmap_sfc)
  ex_albedo_nir_dif_fix = (1.0-ex_albedo_nir_dif) /   &
       (1.0-ex_albedo_nir_dif_fix)

  !=======================================================================
  ! [7] diagnostics section

  !------- save static fields first time only ------
  if (first_static) then

     !------- land fraction ------
     if ( id_land_mask > 0 ) then
        used = send_data ( id_land_mask, Boundary%land_frac, Time )
     endif

     first_static = .false.
  endif

  !------- drag coeff moisture -----------
  if ( id_wind > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_wind, xmap_sfc)
     used = send_data ( id_wind, diag_atm, Time )
  endif
  !------- drag coeff moisture -----------
  if ( id_drag_moist > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_q, xmap_sfc)
     used = send_data ( id_drag_moist, diag_atm, Time )
  endif

  !------- drag coeff heat -----------
  if ( id_drag_heat > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_t, xmap_sfc)
     used = send_data ( id_drag_heat, diag_atm, Time )
  endif
  
  !------- drag coeff momemtum -----------
  if ( id_drag_mom > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_m, xmap_sfc)
     used = send_data ( id_drag_mom, diag_atm, Time )
  endif
  
  !------- roughness moisture -----------
  if ( id_rough_moist > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_rough_moist, xmap_sfc)
     used = send_data ( id_rough_moist, diag_atm, Time )
  endif
  
  !------- roughness heat -----------
  if ( id_rough_heat > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_rough_heat, xmap_sfc)
     used = send_data ( id_rough_heat, diag_atm, Time )
  endif
  
  !------- roughness momemtum -----------
  if ( id_rough_mom > 0 ) then
     used = send_data ( id_rough_mom, Boundary%rough_mom, Time )
  endif
  
  !------- friction velocity -----------
  if ( id_u_star > 0 ) then
     used = send_data ( id_u_star, Boundary%u_star, Time )
  endif
  
  !------- bouyancy -----------
  if ( id_b_star > 0 ) then
     used = send_data ( id_b_star, Boundary%b_star, Time )
  endif

  !------- moisture scale -----------
  if ( id_q_star > 0 ) then
     used = send_data ( id_q_star, Boundary%q_star, Time )
  endif

  !-----------------------------------------------------------------------
  !------ diagnostics for fields at bottom atmospheric level ------
  
  if ( id_t_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_atm, xmap_sfc)
     used = send_data ( id_t_atm, diag_atm, Time )
  endif
  
  if ( id_u_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_u_atm, xmap_sfc)
     used = send_data ( id_u_atm, diag_atm, Time )
  endif
  
  if ( id_v_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_v_atm, xmap_sfc)
     used = send_data ( id_v_atm, diag_atm, Time )
  endif
  
  ! + slm, Mar 25, 2002
  if ( id_q_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_q_atm, xmap_sfc)
     used = send_data ( id_q_atm, diag_atm, Time )
  endif
  ! - slm, Mar 25, 2002
  if ( id_p_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_p_atm, xmap_sfc)
     used = send_data ( id_p_atm, diag_atm, Time )
  endif
  if ( id_z_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_z_atm, xmap_sfc)
     used = send_data ( id_z_atm, diag_atm, Time )
  endif
  if ( id_gust > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_gust, xmap_sfc)
     used = send_data ( id_gust, diag_atm, Time )
  endif

  !-----------------------------------------------------------------------
  !--------- diagnostics for fields at reference level ---------
  
  if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
       id_u_ref > 0 .or. id_v_ref  > 0 .or. &
       id_q_ref > 0 .or. id_q_ref_land > 0 .or. &
       id_t_ref_land > 0 .or. id_rh_ref_land > 0 .or. &
       id_u_ref_land > 0 .or. id_v_ref_land  > 0 ) then
     
     zrefm = z_ref_mom
     zrefh = z_ref_heat
     !      ---- optimize calculation ----
     if ( id_t_ref <= 0 ) zrefh = zrefm
     
     call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
          ex_rough_heat, ex_rough_moist,          &
          ex_u_star, ex_b_star, ex_q_star,        &
          ex_del_m, ex_del_h, ex_del_q, ex_avail  )

     !    ------- reference relative humidity -----------
     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
          id_q_ref > 0 .or. id_q_ref_land >0 ) then
        ex_ref   = ex_q_surf + (ex_q_atm-ex_q_surf) * ex_del_q
        if(id_q_ref > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data(id_q_ref,diag_atm,Time)
        endif
        if(id_q_ref_land > 0) then
           call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data(id_q_ref_land, diag_land, &
                Land%tile_size, Time, mask=Land%mask)
        endif
        where (ex_avail) &
           ex_t_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        call escomp (ex_t_ref, ex_qs_ref)
        ex_qs_ref = d622*ex_qs_ref/(ex_p_surf-d378*ex_qs_ref)
        ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)

        if ( id_rh_ref_land > 0 ) then
           call get_from_xgrid (diag_land,'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data ( id_rh_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if(id_rh_ref > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_rh_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference temp -----------
     if ( id_t_ref > 0 .or. id_t_ref_land > 0 ) then
        ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        if (id_t_ref_land > 0) then
           call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data ( id_t_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_t_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_t_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference u comp -----------
     if ( id_u_ref > 0 .or. id_u_ref_land > 0) then
        ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
        if ( id_u_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_tile_averaged_data ( id_u_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_u_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_u_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference v comp -----------
     if ( id_v_ref > 0 .or. id_v_ref_land > 0 ) then
        ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
        if ( id_v_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_tile_averaged_data ( id_v_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_v_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_v_ref, diag_atm, Time )
        endif
     endif

     !    ------- interp factor for heat ------
     if ( id_del_h > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_h, xmap_sfc)
        used = send_data ( id_del_h, diag_atm, Time )
     endif

     !    ------- interp factor for momentum ------
     if ( id_del_m > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_m, xmap_sfc)
        used = send_data ( id_del_m, diag_atm, Time )
     endif

     !    ------- interp factor for moisture ------
     if ( id_del_q > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_q, xmap_sfc)
        used = send_data ( id_del_q, diag_atm, Time )
     endif

  endif
  ! topographic roughness scale
  if(id_rough_scale>0) then
     call get_from_xgrid (diag_atm, 'ATM',&
          (log(ex_z_atm/ex_rough_mom+1)/log(ex_z_atm/ex_rough_scale+1))**2, xmap_sfc)
     used = send_data(id_rough_scale, diag_atm, Time)
  endif

!Balaji
  call mpp_clock_end(sfcClock)
  call mpp_clock_end(cplClock)

!=======================================================================

end subroutine sfc_boundary_layer
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="flux_down_from_atmos">
!  <OVERVIEW>
!   Returns fluxes and derivatives corrected for the implicit treatment of atmospheric 
!   diffusive fluxes, as well as the increments in the temperature and specific humidity 
!   of the lowest atmospheric layer due to all explicit processes as well as the diffusive 
!   fluxes through the top of this layer. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following elements from Atmos_boundary are used as input: 
!
!        flux_u_atm = zonal wind stress (Pa)  
!        flux_v_atm = meridional wind stress (Pa)
!
!
!    The following elements of Land_boundary are output: 
!
!       flux_t_land = sensible heat flux (W/m2)
!       flux_q_land = specific humidity flux (Kg/(m2 s)
!      flux_lw_land = net longwave flux (W/m2), uncorrected for
!                     changes in surface temperature
!      flux_sw_land = net shortwave flux (W/m2)
!         dhdt_land = derivative of sensible heat flux w.r.t.
!                     surface temperature (on land model grid)  (W/(m2 K)
!         dedt_land = derivative of specific humidity flux w.r.t.
!                     surface temperature (on land model grid)  (Kg/(m2 s K)
!         drdt_land = derivative of upward longwave flux w.r.t.
!                     surface temperature (on land model grid) (W/(m2 K)
!        lprec_land = liquid precipitation, mass for one time step
!                      (Kg/m2)
!        fprec_land = frozen precipitation, mass for one time step
!                      (Kg/m2)
!
!
!    The following elements of Ice_boundary are output: 
!
!        flux_u_ice = zonal wind stress (Pa)
!        flux_v_ice = meridional wind stress (Pa)
!        coszen_ice = cosine of the zenith angle
!
!   </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_down_from_atmos (Time, Atm, Land, Ice, &
!		Atmos_boundary, Land_boundary, Ice_boundary )
!		
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.
!  </IN>
!  <INOUT NAME="Land_boundary" TYPE="atmos_land_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to land.
!  </INOUT>
!  <INOUT NAME="Ice_boundary" TYPE="atmos_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to ice.
!  </INOUT>
!
subroutine flux_down_from_atmos (Time, Atm, Land, Ice, &
     Atmos_boundary, Land_boundary, Ice_boundary )

  type(time_type),       intent(in) :: Time
  type(atmos_data_type), intent(inout) :: Atm
  type(land_data_type),  intent(in) :: Land
  type(ice_data_type),   intent(in) :: Ice
  type(land_ice_atmos_boundary_type),intent(in) :: Atmos_boundary
  type(atmos_land_boundary_type),    intent(inout):: Land_boundary
  type(atmos_ice_boundary_type),     intent(inout):: Ice_boundary

  real, dimension(n_xgrid_sfc) :: ex_flux_sw, ex_flux_lwd, &
       ex_flux_sw_dir,  &
                                    ex_flux_sw_dif,  &
      ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
      ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
       ex_flux_sw_vis, &
       ex_flux_sw_vis_dir, &
       ex_flux_sw_vis_dif, &
       ex_lprec, ex_fprec,      &
       ex_u_star_smooth,        &
       ex_coszen, ex_ice_frac

  real, dimension(n_xgrid_sfc) :: ex_gamma  , ex_dtmass,  &
       ex_delta_t, ex_delta_q, ex_delta_u, ex_delta_v, &
       ex_dflux_t, ex_dflux_q

  real :: ice_frac (size(Ice_boundary%dhdt,1),size(Ice_boundary%dhdt,2),size(Ice_boundary%dhdt,3))
  real :: diag_atm (size(Atmos_boundary%u_flux,1),size(Atmos_boundary%u_flux,2))

  real    :: cp_inv
  logical :: used
  logical :: ov

!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxAtmDnClock)
  ov = .FALSE.
!-----------------------------------------------------------------------
!Balaji: data_override calls moved here from coupler_main            
  call data_override ('ATM', 'flux_sw',  Atm%flux_sw, Time)
  call data_override ('ATM', 'flux_sw_dir',  Atm%flux_sw_dir, Time)
  call data_override ('ATM', 'flux_sw_dif',  Atm%flux_sw_dif, Time)
  call data_override ('ATM', 'flux_sw_down_vis_dir',  Atm%flux_sw_down_vis_dir, Time)
  call data_override ('ATM', 'flux_sw_down_vis_dif',  Atm%flux_sw_down_vis_dif, Time)
  call data_override ('ATM', 'flux_sw_down_total_dir',  Atm%flux_sw_down_total_dir, Time)
  call data_override ('ATM', 'flux_sw_down_total_dif',  Atm%flux_sw_down_total_dif, Time)
  call data_override ('ATM', 'flux_sw_vis',  Atm%flux_sw_vis, Time)
  call data_override ('ATM', 'flux_sw_vis_dir',  Atm%flux_sw_vis_dir, Time)
  call data_override ('ATM', 'flux_sw_vis_dif',  Atm%flux_sw_vis_dif, Time)
  call data_override ('ATM', 'flux_lw',  Atm%flux_lw, Time)
  call data_override ('ATM', 'lprec',    Atm%lprec,   Time)
  call data_override ('ATM', 'fprec',    Atm%fprec,   Time)
  call data_override ('ATM', 'coszen',   Atm%coszen,  Time)
  call data_override ('ATM', 'dtmass',   Atm%Surf_Diff%dtmass, Time)
  call data_override ('ATM', 'delta_t',  Atm%Surf_Diff%delta_t, Time)
  call data_override ('ATM', 'delta_q',  Atm%Surf_Diff%delta_q, Time)
  call data_override ('ATM', 'dflux_t',  Atm%Surf_Diff%dflux_t, Time)
  call data_override ('ATM', 'dflux_q',  Atm%Surf_Diff%dflux_q, Time)
!---- put atmosphere quantities onto exchange grid ----

  call put_to_xgrid (Atm%flux_sw, 'ATM', ex_flux_sw, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_vis, 'ATM', ex_flux_sw_vis, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_dir, 'ATM', ex_flux_sw_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_vis_dir, 'ATM', ex_flux_sw_vis_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_dif, 'ATM', ex_flux_sw_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_vis_dif, 'ATM', ex_flux_sw_vis_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_vis_dir, 'ATM', ex_flux_sw_down_vis_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_total_dir, 'ATM', ex_flux_sw_down_total_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_vis_dif, 'ATM', ex_flux_sw_down_vis_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_total_dif, 'ATM', ex_flux_sw_down_total_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_lw, 'ATM', ex_flux_lwd, xmap_sfc, remap_method=remap_method)
  !  ccc = conservation_check(Atm%lprec, 'ATM', xmap_sfc)
  !  if (mpp_pe()== mpp_root_pe()) print *,'LPREC', ccc

  call put_to_xgrid (Atm%lprec,   'ATM', ex_lprec, xmap_sfc)
  call put_to_xgrid (Atm%fprec,   'ATM', ex_fprec, xmap_sfc)

  call put_to_xgrid (Atm%coszen,  'ATM', ex_coszen, xmap_sfc)


  if(ex_u_star_smooth_bug) then
     call put_to_xgrid (Atmos_boundary%u_star, 'ATM', ex_u_star_smooth, xmap_sfc, remap_method=remap_method)
     ex_u_star = ex_u_star_smooth
  endif


! MOD changed the following two lines to put Atmos%surf_diff%delta_u and v
! on exchange grid instead of the stresses themselves so that only the 
! implicit corrections are filtered through the atmospheric grid not the
! stresses themselves
  call put_to_xgrid (Atm%Surf_Diff%delta_u, 'ATM', ex_delta_u, xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%Surf_Diff%delta_v, 'ATM', ex_delta_v, xmap_sfc, remap_method=remap_method)

  ! MOD update stresses using atmos delta's but derivatives on exchange grid
  ex_flux_u = ex_flux_u + ex_delta_u*ex_dtaudu_atm
  ex_flux_v = ex_flux_v + ex_delta_v*ex_dtaudv_atm

!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----

  ex_flux_sw = ex_flux_sw * ex_albedo_fix

! NEEDED ??
!! IS THIS CORRECT ????
!! CHECK FOR PROPER ALBEDO TO USE in each case
  ex_flux_sw_vis = ex_flux_sw_vis * ex_albedo_vis_dir_fix
  ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_vis_dir_fix
  ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_vis_dif_fix
  ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix
  ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix
!  ex_flux_sw_down_vis_dir = ex_flux_sw_down_vis_dir * ex_albedo_fix
!  ex_flux_sw_down_total_dir = ex_flux_sw_down_total_dir * ex_albedo_fix
!  ex_flux_sw_down_vis_dif = ex_flux_sw_down_vis_dif * ex_albedo_fix
!  ex_flux_sw_down_total_dif = ex_flux_sw_down_total_dif * ex_albedo_fix
 
  deallocate ( ex_albedo_fix )
  deallocate ( ex_albedo_vis_dir_fix )
  deallocate ( ex_albedo_nir_dir_fix )
  deallocate ( ex_albedo_vis_dif_fix )
  deallocate ( ex_albedo_nir_dif_fix )
!----- compute net longwave flux (down-up) -----
  ! (note: lw up already in ex_flux_lw)

  ex_flux_lw = ex_flux_lwd - ex_flux_lw

!-----------------------------------------------------------------------
!----- adjust fluxes for implicit dependence on atmosphere ----


  call put_to_xgrid (Atm%Surf_Diff%dtmass , 'ATM', ex_dtmass , xmap_sfc )
  call put_to_xgrid (Atm%Surf_Diff%delta_t, 'ATM', ex_delta_t, xmap_sfc )
  call put_to_xgrid (Atm%Surf_Diff%delta_q, 'ATM', ex_delta_q, xmap_sfc )
  call put_to_xgrid (Atm%Surf_Diff%dflux_t, 'ATM', ex_dflux_t, xmap_sfc )
  call put_to_xgrid (Atm%Surf_Diff%dflux_q, 'ATM', ex_dflux_q, xmap_sfc )

  cp_inv = 1.0/cp_air

  where(ex_avail)

     ! temperature

     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_t + ex_dhdt_atm*cp_inv))
     ex_e_t_n      =  ex_dtmass*ex_dhdt_surf*cp_inv*ex_gamma
     ex_f_t_delt_n = (ex_delta_t + ex_dtmass * ex_flux_t*cp_inv) * ex_gamma    
     
     ex_flux_t     =  ex_flux_t        + ex_dhdt_atm * ex_f_t_delt_n 
     ex_dhdt_surf  =  ex_dhdt_surf     + ex_dhdt_atm * ex_e_t_n   

     ! moisture
     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_q + ex_dedq_atm))

! here it looks like two derivatives with different units are added together,
! but in fact they are not: ex_dedt_surf and ex_dedq_surf defined in complimentary
! regions of exchange grid, so that if one of them is not zero the other is, and
! vice versa.
     ex_e_q_n      =  ex_dtmass*(ex_dedt_surf+ex_dedq_surf) * ex_gamma

     ex_f_q_delt_n = (ex_delta_q  + ex_dtmass * ex_flux_q) * ex_gamma    
     
     ex_flux_q     =  ex_flux_q    + ex_dedq_atm * ex_f_q_delt_n 
     ex_dedt_surf  =  ex_dedt_surf + ex_dedq_atm * ex_e_q_n
     ex_dedq_surf  =  ex_dedq_surf + ex_dedq_atm * ex_e_q_n
  endwhere

!-----------------------------------------------------------------------
!---- output fields on the land grid -------

  call get_from_xgrid (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Land_boundary%q_flux,  'LND', ex_flux_q,    xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_vis_dir, 'LND', ex_flux_sw_down_vis_dir,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_total_dir, 'LND', ex_flux_sw_down_total_dir,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_vis_dif, 'LND', ex_flux_sw_down_vis_dif,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_total_dif, 'LND', ex_flux_sw_down_total_dif,   xmap_sfc)
  call get_from_xgrid (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
  call get_from_xgrid (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%dedt,    'LND', ex_dedt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%dedq,    'LND', ex_dedq_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

  if(associated(Land_boundary%drag_q)) then
     call get_from_xgrid (Land_boundary%drag_q, 'LND', ex_drag_q,    xmap_sfc)
     call data_override('LND', 'drag_q', Land_boundary%drag_q,  Time )
  endif
  if(associated(Land_boundary%lwdn_flux)) then
     call get_from_xgrid (Land_boundary%lwdn_flux, 'LND', ex_flux_lwd, xmap_sfc)
     call data_override('LND', 'lwdn_flux', Land_boundary%lwdn_flux, Time )
  endif
  if(associated(Land_boundary%cd_m)) then
     call get_from_xgrid (Land_boundary%cd_m, 'LND', ex_cd_m, xmap_sfc)
     call data_override('LND', 'cd_m', Land_boundary%cd_m, Time )
  endif
  if(associated(Land_boundary%cd_t)) then
     call get_from_xgrid (Land_boundary%cd_t, 'LND', ex_cd_t, xmap_sfc)
     call data_override('LND', 'cd_t', Land_boundary%cd_t, Time )
  endif
  if(associated(Land_boundary%bstar)) then
     call get_from_xgrid (Land_boundary%bstar, 'LND', ex_b_star, xmap_sfc)
     call data_override('LND', 'bstar',  Land_boundary%bstar, Time )
  endif
  if(associated(Land_boundary%ustar)) then
     call get_from_xgrid (Land_boundary%ustar, 'LND', ex_u_star, xmap_sfc)
     call data_override('LND', 'ustar',  Land_boundary%ustar, Time )
  endif
  if(associated(Land_boundary%wind)) then
     call get_from_xgrid (Land_boundary%wind, 'LND', ex_wind, xmap_sfc)
     call data_override('LND', 'wind',  Land_boundary%wind, Time )
  endif
  if(associated(Land_boundary%z_bot)) then
     call get_from_xgrid (Land_boundary%z_bot, 'LND', ex_z_atm, xmap_sfc)
     call data_override('LND', 'bstar',  Land_boundary%bstar, Time )
  endif

!  current time is Time: is that ok? not available in land_data_type
!Balaji: data_override calls moved here from coupler_main
  call data_override('LND', 't_flux',  Land_boundary%t_flux,  Time )
  call data_override('LND', 'q_flux',  Land_boundary%q_flux,  Time )
  call data_override('LND', 'lw_flux', Land_boundary%lw_flux, Time )
  call data_override('LND', 'sw_flux', Land_boundary%sw_flux, Time )
  call data_override('LND', 'sw_flux_down_vis_dir', Land_boundary%sw_flux_down_vis_dir, Time )
  call data_override('LND', 'sw_flux_down_total_dir', Land_boundary%sw_flux_down_total_dir, Time )
  call data_override('LND', 'sw_flux_down_vis_dif', Land_boundary%sw_flux_down_vis_dif, Time )
  call data_override('LND', 'sw_flux_down_total_dif', Land_boundary%sw_flux_down_total_dif, Time )
  
  call data_override('LND', 'lprec',   Land_boundary%lprec,   Time )
  call data_override('LND', 'fprec',   Land_boundary%fprec,   Time )
  call data_override('LND', 'dhdt',    Land_boundary%dhdt,    Time )
  call data_override('LND', 'dedt',    Land_boundary%dedt,    Time )
  call data_override('LND', 'dedq',    Land_boundary%dedq,    Time )
  call data_override('LND', 'drdt',    Land_boundary%drdt,    Time )
  call data_override('LND', 'p_surf',  Land_boundary%p_surf,  Time )

!-----------------------------------------------------------------------
!---- output fields on the ice grid -------

  call get_from_xgrid (Ice_boundary%t_flux,   'OCN', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%q_flux,   'OCN', ex_flux_q,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux,  'OCN', ex_flux_sw,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_vis,  'OCN', ex_flux_sw_vis,xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_dir,  'OCN', ex_flux_sw_dir,xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_vis_dir,  'OCN', ex_flux_sw_vis_dir,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_dif,  'OCN', ex_flux_sw_dif,xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_vis_dif,  'OCN', ex_flux_sw_vis_dif,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%lw_flux,  'OCN', ex_flux_lw,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%dhdt,     'OCN', ex_dhdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%dedt,     'OCN', ex_dedt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%drdt,     'OCN', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%lprec,    'OCN', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Ice_boundary%fprec,    'OCN', ex_fprec,     xmap_sfc)
  call get_from_xgrid (Ice_boundary%u_flux,   'OCN', ex_flux_u,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%v_flux,   'OCN', ex_flux_v,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%u_star,   'OCN', ex_u_star,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%coszen,   'OCN', ex_coszen,    xmap_sfc)
!Balaji: data_override calls moved here from coupler_main
  call data_override('ICE', 'u_flux', Ice_boundary%u_flux,  Time)
  call data_override('ICE', 'v_flux', Ice_boundary%v_flux,  Time)
  call data_override('ICE', 't_flux', Ice_boundary%t_flux,  Time)
  call data_override('ICE', 'q_flux', Ice_boundary%q_flux,  Time)
  call data_override('ICE', 'lw_flux',Ice_boundary%lw_flux, Time)
  call data_override('ICE', 'lw_flux_dn',Ice_boundary%lw_flux, Time, override=ov)
  if (ov) then
    Ice_boundary%lw_flux = Ice_boundary%lw_flux - stefan*Ice%t_surf**4
  endif
  call data_override('ICE', 'sw_flux',Ice_boundary%sw_flux, Time)
  call data_override('ICE', 'sw_flux_vis',Ice_boundary%sw_flux_vis, Time)
  call data_override('ICE', 'sw_flux_vis_dir',Ice_boundary%sw_flux_vis_dir, Time)
  call data_override('ICE', 'sw_flux_dir',Ice_boundary%sw_flux_dir, Time, override=ov)
  call data_override('ICE', 'sw_flux_vis_dif',Ice_boundary%sw_flux_vis_dif, Time)
  call data_override('ICE', 'sw_flux_dif',Ice_boundary%sw_flux_dif, Time, override=ov)
  call data_override('ICE', 'sw_flux_dn',Ice_boundary%sw_flux, Time, override=ov)
!! ANY CHANGE NEEDED HERE FOR sw_flux_vis ??
! this converts net fluxes to downward fluxes
  if (ov) then
    Ice_boundary%sw_flux = Ice_boundary%sw_flux*(1-Ice%albedo)
! ?? NEEDED:    Ice_boundary%sw_flux_vis = Ice_boundary%sw_flux_vis*(1-Ice%albedo_vis_dir or dif, must choose ??)
! ?? NEEDED:    Ice_boundary%sw_flux_vis_dir = Ice_boundary%sw_flux_vis_dir*(1-Ice%albedo_vis_dir)
! ?? NEEDED:    Ice_boundary%sw_flux_vis_dif = Ice_boundary%sw_flux_vis_dif*(1-Ice%albedo_vis_dif)
! ?? NEEDED:    Ice_boundary%sw_flux_dir = Ice_boundary%sw_flux_dir*(1-Ice%albedo_( select vis or nir) dir)
! ?? NEEDED:    Ice_boundary%sw_flux_dif = Ice_boundary%sw_flux_dif*(1-Ice%albedo_ ( select vis or nir ) dif)
  endif
  call data_override('ICE', 'lprec',  Ice_boundary%lprec,   Time)
  call data_override('ICE', 'fprec',  Ice_boundary%fprec,   Time)
  call data_override('ICE', 'dhdt',   Ice_boundary%dhdt,    Time)
  call data_override('ICE', 'dedt',   Ice_boundary%dedt,    Time)
  call data_override('ICE', 'drdt',   Ice_boundary%drdt,    Time)
  call data_override('ICE', 'coszen', Ice_boundary%coszen,  Time)
  call data_override('ICE', 'p',      Ice_boundary%p,       Time)

  deallocate ( ex_flux_u, ex_flux_v, ex_dtaudu_atm, ex_dtaudv_atm)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- zonal wind stress -----------
  if ( id_u_flux > 0 ) then
     used = send_data ( id_u_flux, Atmos_boundary%u_flux, Time )
  endif

  !------- meridional wind stress -----------
  if ( id_v_flux > 0 ) then
     used = send_data ( id_v_flux, Atmos_boundary%v_flux, Time )
  endif

!Balaji
  call mpp_clock_end(fluxAtmDnClock)
  call mpp_clock_end(cplClock)
!=======================================================================

  end subroutine flux_down_from_atmos
! </SUBROUTINE>

!#######################################################################
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! flux_land_to_ice - translate runoff from land to ice grids                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! <SUBROUTINE NAME="flux_land_to_ice">
!  <OVERVIEW>
!   Conservative transfer of water and snow discharge from the land model to sea ice/ocean model. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following elements are transferred from the Land to the Land_ice_boundary: 
!
!        discharge --> runoff (kg/m2)
!        discharge_snow --> calving (kg/m2)
!
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_land_to_ice(Time, Land, Ice, Boundary )
!		
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <INOUT NAME="Boundary" TYPE="land_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from land to ice.
!  </INOUT>
!
subroutine flux_land_to_ice( Time, Land, Ice, Boundary )
  type(time_type),               intent(in) :: Time
  type(land_data_type),          intent(in) :: Land
  type(ice_data_type),           intent(in) :: Ice
!real, dimension(:,:), intent(out) :: runoff_ice, calving_ice
  type(land_ice_boundary_type),  intent(inout):: Boundary
  
  real, dimension(n_xgrid_runoff) :: ex_runoff, ex_calving
  real, dimension(size(Boundary%runoff,1),size(Boundary%runoff,2),1) :: ice_buf
!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxLandIceClock)

  ! ccc = conservation_check(Land%discharge, 'LND', xmap_runoff)
  ! if (mpp_pe()==mpp_root_pe()) print *,'RUNOFF', ccc

  call put_to_xgrid ( Land%discharge,      'LND', ex_runoff,  xmap_runoff)
  call put_to_xgrid ( Land%discharge_snow, 'LND', ex_calving, xmap_runoff)
  call get_from_xgrid (ice_buf, 'OCN', ex_runoff,  xmap_runoff)
  Boundary%runoff = ice_buf(:,:,1);
  call get_from_xgrid (ice_buf, 'OCN', ex_calving, xmap_runoff)
  Boundary%calving = ice_buf(:,:,1);
!Balaji
  call data_override('ICE', 'runoff' , Boundary%runoff , Time)
  call data_override('ICE', 'calving', Boundary%calving, Time)
  call mpp_clock_end(fluxLandIceClock)
  call mpp_clock_end(cplClock)

end subroutine flux_land_to_ice
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_ice_to_ocean">
!  <OVERVIEW>
!   Takes the ice model state (fluxes at the bottom of the ice) and interpolates it to the ocean model grid. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!   The following quantities are transferred from the Ice to the ice_ocean_boundary_type: 
!
!       flux_u = zonal wind stress (Pa)
!       flux_v = meridional wind stress (Pa)
!       flux_t = sensible heat flux (W/m2)
!       flux_q = specific humidity flux (Kg/m2/s)
!    flux_salt = salt flux (Kg/m2/s)
!      flux_sw = net (down-up) shortwave flux (W/m2)
!      flux_lw = net (down-up) longwave flux (W/m2)
!        lprec = mass of liquid precipitation since last
!                      time step (Kg/m2)
!        fprec = mass of frozen precipitation since last
!                time step (Kg/m2)
!       runoff = mass (?) of runoff since last time step
!                       (Kg/m2)
!       p_surf = surface pressure (Pa)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ice_to_ocean ( Time, Ice, Ocean, Boundary )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME=" Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Ocean" TYPE="ocean_data_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <INOUT NAME="Boundary" TYPE="ice_ocean_boundary_type">
!   A derived data type to specify properties and fluxes passed from ice to ocean.
!  </INOUT>
!
subroutine flux_ice_to_ocean ( Time, Ice, Ocean, Boundary )
  type(time_type),        intent(in) :: Time
  type(ice_data_type),   intent(in)  :: Ice
  type(ocean_data_type), intent(in)  :: Ocean
!  real, dimension(:,:),   intent(out) :: flux_u_ocean,  flux_v_ocean,  &
!                                         flux_t_ocean,  flux_q_ocean,  &
!                                         flux_sw_ocean, flux_lw_ocean, &
!                                         lprec_ocean,   fprec_ocean,   &
!                                         runoff_ocean,  calving_ocean, &
!                                         flux_salt_ocean, p_surf_ocean
  type(ice_ocean_boundary_type), intent(inout) :: Boundary
  integer(8) :: chksum

!Balaji
  call mpp_clock_begin(cplOcnClock)
  call mpp_clock_begin(fluxIceOceanClock)
  select case (Boundary%xtype)
  case(DIRECT)
     !same grid and domain decomp for ocean and ice    
     if( ASSOCIATED(Boundary%u_flux   ) )Boundary%u_flux    = Ice%flux_u
     if( ASSOCIATED(Boundary%v_flux   ) )Boundary%v_flux    = Ice%flux_v
     if( ASSOCIATED(Boundary%t_flux   ) )Boundary%t_flux    = Ice%flux_t
     if( ASSOCIATED(Boundary%q_flux   ) )Boundary%q_flux    = Ice%flux_q
     if( ASSOCIATED(Boundary%salt_flux) )Boundary%salt_flux = Ice%flux_salt
     if( ASSOCIATED(Boundary%sw_flux  ) )Boundary%sw_flux   = Ice%flux_sw
     if( ASSOCIATED(Boundary%sw_flux_vis  ) )Boundary%sw_flux_vis   = Ice%flux_sw_vis
     if( ASSOCIATED(Boundary%sw_flux_dir  ) )Boundary%sw_flux_dir   = Ice%flux_sw_dir
     if( ASSOCIATED(Boundary%sw_flux_dif  ) )Boundary%sw_flux_dif   = Ice%flux_sw_dif
     if( ASSOCIATED(Boundary%sw_flux_vis_dir  ) )Boundary%sw_flux_vis_dir   = Ice%flux_sw_vis_dir
     if( ASSOCIATED(Boundary%sw_flux_vis_dif  ) )Boundary%sw_flux_vis_dif   = Ice%flux_sw_vis_dif
     if( ASSOCIATED(Boundary%lw_flux  ) )Boundary%lw_flux   = Ice%flux_lw
     if( ASSOCIATED(Boundary%lprec    ) )Boundary%lprec     = Ice%lprec
     if( ASSOCIATED(Boundary%fprec    ) )Boundary%fprec     = Ice%fprec
     if( ASSOCIATED(Boundary%runoff   ) )Boundary%runoff    = Ice%runoff
     if( ASSOCIATED(Boundary%calving  ) )Boundary%calving   = Ice%calving
     if( ASSOCIATED(Boundary%p        ) )Boundary%p         = Ice%p_surf
  case(REDIST)
     !same grid, different domain decomp for ocean and ice    
     if (ASSOCIATED(Boundary%u_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_u, Ocean%Domain, Boundary%u_flux)
     if (ASSOCIATED(Boundary%v_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_v, Ocean%Domain, Boundary%v_flux)
     if (ASSOCIATED(Boundary%t_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_t, Ocean%Domain, Boundary%t_flux)
     if (ASSOCIATED(Boundary%q_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_q, Ocean%Domain, Boundary%q_flux)
     if (ASSOCIATED(Boundary%salt_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_salt, Ocean%Domain, Boundary%salt_flux)
     if (ASSOCIATED(Boundary%sw_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw, Ocean%Domain, Boundary%sw_flux)
     if (ASSOCIATED(Boundary%sw_flux_vis)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_vis, Ocean%Domain, Boundary%sw_flux_vis)
     if (ASSOCIATED(Boundary%sw_flux_dir)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_dir, Ocean%Domain, Boundary%sw_flux_dir)
     if (ASSOCIATED(Boundary%sw_flux_dif)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_dif, Ocean%Domain, Boundary%sw_flux_dif)
     if (ASSOCIATED(Boundary%sw_flux_vis_dir)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_vis_dir, Ocean%Domain, Boundary%sw_flux_vis_dir)
     if (ASSOCIATED(Boundary%sw_flux_vis_dif)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_vis_dif, Ocean%Domain, Boundary%sw_flux_vis_dif)
     if (ASSOCIATED(Boundary%lw_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_lw, Ocean%Domain, Boundary%lw_flux)
     if (ASSOCIATED(Boundary%lprec)) &
          call mpp_redistribute(Ice%Domain, Ice%lprec, Ocean%Domain, Boundary%lprec)
     if (ASSOCIATED(Boundary%fprec)) &
          call mpp_redistribute(Ice%Domain, Ice%fprec, Ocean%Domain, Boundary%fprec)
     if (ASSOCIATED(Boundary%runoff)) &
          call mpp_redistribute(Ice%Domain, Ice%runoff, Ocean%Domain, Boundary%runoff)
     if (ASSOCIATED(Boundary%calving)) &
          call mpp_redistribute(Ice%Domain, Ice%calving, Ocean%Domain, Boundary%calving)
     if (ASSOCIATED(Boundary%p)) &
          call mpp_redistribute(Ice%Domain, Ice%p_surf, Ocean%Domain, Boundary%p)
  case DEFAULT
!   <ERROR MSG="Boundary%xtype must be DIRECT or REDIST." STATUS="FATAL">
!      The value of variable xtype of ice_ocean_boundary_type data must be DIRECT or REDIST.
!   </ERROR>
     call mpp_error( FATAL, 'FLUX_ICE_TO_OCEAN: Boundary%xtype must be DIRECT or REDIST.' )
  end select
!Balaji: moved data_override calls here from coupler_main
  call data_override('OCN', 'u_flux',    Boundary%u_flux   , Time )
  call data_override('OCN', 'v_flux',    Boundary%v_flux   , Time )
  call data_override('OCN', 't_flux',    Boundary%t_flux   , Time )
  call data_override('OCN', 'q_flux',    Boundary%q_flux   , Time )
  call data_override('OCN', 'salt_flux', Boundary%salt_flux, Time )
  call data_override('OCN', 'lw_flux',   Boundary%lw_flux  , Time )
  call data_override('OCN', 'sw_flux',   Boundary%sw_flux  , Time )
  call data_override('OCN', 'sw_flux_vis',   Boundary%sw_flux_vis  , Time )
  call data_override('OCN', 'sw_flux_dir',   Boundary%sw_flux_dir  , Time )
  call data_override('OCN', 'sw_flux_dif',   Boundary%sw_flux_dif  , Time )
  call data_override('OCN', 'sw_flux_vis_dir',   Boundary%sw_flux_vis_dir  , Time )
  call data_override('OCN', 'sw_flux_vis_dif',   Boundary%sw_flux_vis_dif  , Time )
  call data_override('OCN', 'lprec',     Boundary%lprec    , Time )
  call data_override('OCN', 'fprec',     Boundary%fprec    , Time )
  call data_override('OCN', 'runoff',    Boundary%runoff   , Time )
  call data_override('OCN', 'calving',   Boundary%calving  , Time )
  call data_override('OCN', 'p',         Boundary%p        , Time )
!Balaji
  call mpp_clock_end(fluxIceOceanClock)
  call mpp_clock_end(cplOcnClock)
!-----------------------------------------------------------------------

  end subroutine flux_ice_to_ocean
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_ocean_to_ice">
!  <OVERVIEW>
!   Takes the ocean model state and interpolates it onto the bottom of the ice. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following quantities are transferred from the Ocean to the ocean_ice_boundary_type: 
!
!        t_surf = surface temperature (deg K)
!        frazil = frazil (???)
!        u_surf = zonal ocean current/ice motion (m/s)
!        v_surf = meridional ocean current/ice motion (m/s
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ocean_to_ice ( Time, Ocean, Ice, Boundary)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME=" Ocean" TYPE="ocean_data_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <INOUT NAME="Boundary" TYPE="ocean_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from ocean to ice.
!  </INOUT>
!
subroutine flux_ocean_to_ice ( Time, Ocean, Ice, Boundary )
  type(time_type),       intent(in)  :: Time
  type(ocean_data_type), intent(in)  :: Ocean
  type(ice_data_type),   intent(in)  :: Ice
!  real, dimension(:,:),   intent(out) :: t_surf_ice, u_surf_ice, v_surf_ice, &
!                                         frazil_ice, s_surf_ice, sea_lev_ice
  type(ocean_ice_boundary_type), intent(inout) :: Boundary
  real, dimension(size(Boundary%t,1),size(Boundary%t,2),size(Ice%part_size,3)) &
       :: ice_frac
  real, dimension(:), allocatable :: ex_ice_frac
  real, dimension(ni_atm, nj_atm) :: diag_atm
  logical :: used

!Balaji
  call mpp_clock_begin(cplOcnClock)
  call mpp_clock_begin(fluxOceanIceClock)

  select case (Boundary%xtype)
  case(DIRECT)
     !same grid and domain decomp for ocean and ice    
     if( ASSOCIATED(Boundary%u) )Boundary%u = Ocean%u_surf
     if( ASSOCIATED(Boundary%v) )Boundary%v = Ocean%v_surf
     if( ASSOCIATED(Boundary%t) )Boundary%t = Ocean%t_surf
     if( ASSOCIATED(Boundary%s) )Boundary%s = Ocean%s_surf
     if( ASSOCIATED(Boundary%frazil) )Boundary%frazil = Ocean%frazil
     if( ASSOCIATED(Boundary%sea_level) )Boundary%sea_level = Ocean%sea_lev
  case(REDIST)
     !same grid, different domain decomp for ocean and ice    
     if( ASSOCIATED(Boundary%u) )call mpp_redistribute(Ocean%Domain, Ocean%u_surf, Ice%Domain, Boundary%u)
     if( ASSOCIATED(Boundary%v) )call mpp_redistribute(Ocean%Domain, Ocean%v_surf, Ice%Domain, Boundary%v)
     if( ASSOCIATED(Boundary%t) )call mpp_redistribute(Ocean%Domain, Ocean%t_surf, Ice%Domain, Boundary%t)
     if( ASSOCIATED(Boundary%s) )call mpp_redistribute(Ocean%Domain, Ocean%s_surf, Ice%Domain, Boundary%s)
     if( ASSOCIATED(Boundary%frazil) )call mpp_redistribute(Ocean%Domain, Ocean%frazil, Ice%Domain, Boundary%frazil)
     if( ASSOCIATED(Boundary%sea_level) )call mpp_redistribute(Ocean%Domain, Ocean%sea_lev, Ice%Domain, Boundary%sea_level)
  case DEFAULT
!   <ERROR MSG="Boundary%xtype must be DIRECT or REDIST." STATUS="FATAL">
!     The value of variable xtype of ice_ocean_boundary_type data must be DIRECT or REDIST.
!   </ERROR>
     call mpp_error( FATAL, 'FLUX_OCEAN_TO_ICE: Boundary%xtype must be DIRECT or REDIST.' )
  end select
!Balaji: data_override moved here from coupler_main
  call data_override('ICE', 'u',         Boundary%u,         Time)
  call data_override('ICE', 'v',         Boundary%v,         Time)
  call data_override('ICE', 't',         Boundary%t,         Time)
  call data_override('ICE', 's',         Boundary%s,         Time)
  call data_override('ICE', 'frazil',    Boundary%frazil,    Time)
  call data_override('ICE', 'sea_level', Boundary%sea_level, Time)
  if ( id_ice_mask > 0 ) then
     allocate ( ex_ice_frac(n_xgrid_sfc) )
     ice_frac        = 1.
     ice_frac(:,:,1) = 0.
     ex_ice_frac     = 0.
     call put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
     call get_from_xgrid (diag_atm, 'ATM', ex_ice_frac, xmap_sfc)
     used = send_data ( id_ice_mask, diag_atm, Time )
     deallocate ( ex_ice_frac )
  endif

!Balaji
  call mpp_clock_end(fluxOceanIceClock)
  call mpp_clock_end(cplOcnClock)
!-----------------------------------------------------------------------

  end subroutine flux_ocean_to_ice
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="generate_sfc_xgrid">
!  <OVERVIEW>
!   Optimizes the exchange grids by eliminating land and ice partitions with no data. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   Optimizes the exchange grids by eliminating land and ice partitions with no data. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call generate_sfc_xgrid( Land, Ice )
!		
!  </TEMPLATE>
!  <IN NAME=" Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!  A derived data type to specify ice boundary data.
!  </IN>
!
subroutine generate_sfc_xgrid( Land, Ice )
! subroutine to regenerate exchange grid eliminating side 2 tiles with 0 frac area
    type(land_data_type), intent(in) :: Land
    type(ice_data_type),  intent(in) :: Ice

    integer :: isc, iec, jsc, jec

!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(regenClock)

  call mpp_get_compute_domain(Ice%Domain, isc, iec, jsc, jec)

  call set_frac_area (Ice%part_size(isc:iec,jsc:jec,:) , 'OCN', xmap_sfc)
  call set_frac_area (Land%tile_size, 'LND', xmap_sfc)
  n_xgrid_sfc = max(xgrid_count(xmap_sfc),1)

!Balaji
  call mpp_clock_end(regenClock)
  call mpp_clock_end(cplClock)
  return
end subroutine generate_sfc_xgrid
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_up_to_atmos">
!  <OVERVIEW>
!   Corrects the fluxes for consistency with the new surface temperatures in land 
!   and ice models.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Corrects the fluxes for consistency with the new surface temperatures in land 
!   and ice models. Final increments for temperature and specific humidity in the 
!   lowest atmospheric layer are computed and returned to the atmospheric model
!   so that it can finalize the increments in the rest of the atmosphere. 
!  <PRE>
!
!   The following elements of the land_ice_atmos_boundary_type are computed:
!        dt_t  = temperature change at the lowest
!                 atmospheric level (deg k)
!        dt_q  = specific humidity change at the lowest
!                 atmospheric level (kg/kg)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_up_to_atmos ( Time, Land, Ice, Boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   Current time.
!  </IN>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </INOUT>
!  <INOUT NAME="Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice. 
!  </INOUT>
!
subroutine flux_up_to_atmos ( Time, Land, Ice, Boundary )

  type(time_type),      intent(in)  :: Time
  type(land_data_type), intent(inout)  :: Land
  type(ice_data_type),  intent(inout)  :: Ice
! real, dimension(:,:),   intent(out) :: dt_t_atm, dt_q_atm
  type(land_ice_atmos_boundary_type), intent(inout) :: Boundary
  real, dimension(n_xgrid_sfc) :: ex_t_surf_new, ex_dt_t_surf,  &
       ex_dt_t, ex_dt_q, &
       ex_delta_t_n,  ex_delta_q_n, &
       ! + slm, Mar 20 2002
       ex_t_ca_new,   ex_dt_t_ca,   &
       ex_q_surf_new, ex_dt_q_surf
       ! - slm, Mar 20 2002

  real, dimension(size(Boundary%dt_t,1),size(Boundary%dt_t,2)) :: diag_atm, &
       evap_atm
  logical :: used

!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxAtmUpClock)
  !-----------------------------------------------------------------------
!Balaji: data_override calls moved here from coupler_main
  call data_override ( 'ICE', 't_surf', Ice%t_surf,  Time)
  call data_override ( 'LND', 't_ca',   Land%t_ca,   Time)
  call data_override ( 'LND', 't_surf', Land%t_surf, Time)
  call data_override ( 'LND', 'q_ca',   Land%q_ca,   Time)

  !----- compute surface temperature change -----

  ex_t_surf_new = 200.0

  call put_to_xgrid (Ice%t_surf,  'OCN', ex_t_surf_new, xmap_sfc)
  ex_t_ca_new = ex_t_surf_new  ! since it is the same thing over oceans
  call put_to_xgrid (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
  call put_to_xgrid (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)

  call escomp(ex_t_ca_new, ex_q_surf_new)
  ex_q_surf_new  = d622*ex_q_surf_new/(ex_p_surf-d378*ex_q_surf_new) 
  call put_to_xgrid (Land%q_ca, 'LND', ex_q_surf_new, xmap_sfc)

  where (ex_avail)
     ex_dt_t_ca   = ex_t_ca_new   - ex_t_ca   ! changes in near-surface T
     ex_dt_t_surf = ex_t_surf_new - ex_t_surf ! changes in radiative T
     ex_dt_q_surf = ex_q_surf_new - ex_q_surf ! changes in near-surface q
  endwhere

  !-----------------------------------------------------------------------
  !-----  adjust fluxes and atmospheric increments for 
  !-----  implicit dependence on surface temperature -----

  ex_delta_t_n = 0.0
  ex_delta_q_n = 0.0

  where(ex_avail)
     ex_flux_t     = ex_flux_t  + ex_dt_t_ca   * ex_dhdt_surf
     ex_flux_lw    = ex_flux_lw - ex_dt_t_surf * ex_drdt_surf
     ex_delta_t_n  = ex_f_t_delt_n  + ex_dt_t_ca*ex_e_t_n

     ! Note that in the expressions below ex_e_q_n used to relate changes
     ! of humidity on the lower atmos level to changes of sfc temperature or 
     ! near-sfc humidity; it may seem that the units are mixed up. However
     ! it is not the case since for each exchange grid cell ex_e_q_n is
     ! defined either in terms of sfc temperature or near-sfc humidity
     
     where (ex_land) 
        ex_delta_q_n  = ex_f_q_delt_n + ex_dt_q_surf * ex_e_q_n
        ex_flux_q     = ex_flux_q     + ex_dt_q_surf * ex_dedq_surf
     elsewhere
        ! note that in this region (over ocean) ex_dt_t_surf == ex_dt_t_ca
        ex_delta_q_n  = ex_f_q_delt_n + ex_dt_t_surf * ex_e_q_n
        ex_flux_q     = ex_flux_q     + ex_dt_t_surf * ex_dedt_surf
     endwhere
  endwhere

  !-----------------------------------------------------------------------
  !---- get mean quantites on atmospheric grid ----

  call get_from_xgrid (Boundary%dt_t, 'ATM', ex_delta_t_n, xmap_sfc)
  call get_from_xgrid (Boundary%dt_q, 'ATM', ex_delta_q_n, xmap_sfc)

  !  ---- always get evaporation for diagnostic purposes ----

  call get_from_xgrid (evap_atm,      'ATM', ex_flux_q,    xmap_sfc)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- new surface temperature -----------
  if ( id_t_surf > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
     used = send_data ( id_t_surf, diag_atm, Time )
  endif


  ! + slm, Mar 27 2002
  ! ------ new canopy temperature --------
  !   NOTE, that in the particular case of LM2 t_ca is identical to t_surf,
  !   but this will be changed in future version of the land madel
  if ( id_t_ca > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_ca_new, xmap_sfc)
     used = send_data ( id_t_ca, diag_atm, Time )
  endif

  ! ---- new surface humidity ---------
  if ( id_q_surf > 0 ) then
     call get_from_xgrid(diag_atm, 'ATM', ex_q_surf_new, xmap_sfc )
     used = send_data ( id_q_surf, diag_atm, Time )
  endif
  ! - slm, Mar 27 2002

  !------- sensible heat flux -----------
  if ( id_t_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_t, xmap_sfc)
     used = send_data ( id_t_flux, diag_atm, Time )
  endif

!------- latent heat flux (see below) -----------

!  if ( id_q_flux > 0 ) then
!!!!  latent  = hlv
!!!!  where (ice_or_snow) latent = hlf + hlv
!     call get_from_xgrid (diag_atm, 'ATM', ex_flux_q, xmap_sfc)
!     used = send_data ( id_q_flux, diag_atm, Time )
!  endif

  !------- net longwave flux -----------
  if ( id_r_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_lw, xmap_sfc)
     used = send_data ( id_r_flux, diag_atm, Time )
  endif

  !------- evaporation rate -----------
  
  if ( id_q_flux > 0 ) then
     used = send_data ( id_q_flux, evap_atm, Time )
  endif

  !-----------------------------------------------------------------------
  !---- accumulate global integral of evaporation (mm/day) -----

  call sum_diag_integral_field ('evap', evap_atm*86400.)

  !=======================================================================
  !---- deallocate module storage ----
  deallocate ( &
       ex_t_surf   ,  &
       ex_p_surf   ,  &
       ex_t_ca     ,  &
       ex_q_surf   ,  &
       ex_dhdt_surf,  &
       ex_dedt_surf,  &
       ex_dedq_surf,  &
       ex_drdt_surf,  &
       ex_dhdt_atm ,  &
       ex_dedq_atm ,  &
       ex_flux_t   ,  &
       ex_flux_q   ,  &
       ex_flux_lw  ,  &
       ex_drag_q   ,  &
       ex_avail    ,  &
       ex_f_t_delt_n, &
       ex_f_q_delt_n, &
       ex_e_t_n    ,  &
       ex_e_q_n    ,  &
       ! values added for LM3
       ex_cd_t     ,  &
       ex_cd_m     ,  &
       ex_b_star   ,  &
       ex_u_star   ,  &
       ex_wind     ,  &
       ex_z_atm    ,  &

       ex_land        )


!Balaji
  call mpp_clock_end(fluxAtmUpClock)
  call mpp_clock_end(cplClock)

!-----------------------------------------------------------------------

end subroutine flux_up_to_atmos
! </SUBROUTINE>

!#######################################################################

subroutine put_logical_to_real (mask, id, ex_mask, xmap)

  logical         , intent(in)    :: mask(:,:,:)
  character(len=3), intent(in)    :: id
  real            , intent(inout) :: ex_mask(:)
  type(xmap_type), intent(inout) :: xmap

  !-----------------------------------------------------------------------
  !    puts land or ice model masks (with partitions) onto the
  !    exchange grid as a real array (1.=true, 0.=false)
  !-----------------------------------------------------------------------

  real, dimension(size(mask,1),size(mask,2),size(mask,3)) :: rmask
  
  where (mask)
     rmask = 1.0
  elsewhere
     rmask = 0.0
  endwhere

  call put_to_xgrid(rmask, id, ex_mask, xmap)

end subroutine put_logical_to_real

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes, land_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)
  integer,         intent(in) :: land_axes(2)

  integer :: iref
  character(len=6) :: label_zm, label_zh
  real, dimension(2) :: trange = (/  100., 400. /), &
       vrange = (/ -400., 400. /), &
       frange = (/ -0.01, 1.01 /)
!-----------------------------------------------------------------------
!  initializes diagnostic fields that may be output from this module
!  (the id numbers may be referenced anywhere in this module)
!-----------------------------------------------------------------------

  !------ labels for diagnostics -------
  !  (z_ref_mom, z_ref_heat are namelist variables)

  iref = int(z_ref_mom+0.5)
  if ( real(iref) == z_ref_mom ) then
     write (label_zm,105) iref
     if (iref < 10) write (label_zm,100) iref
  else
     write (label_zm,110) z_ref_mom
  endif

  iref = int(z_ref_heat+0.5)
  if ( real(iref) == z_ref_heat ) then
     write (label_zh,105) iref
     if (iref < 10) write (label_zh,100) iref
  else
     write (label_zh,110) z_ref_heat
  endif

100 format (i1,' m',3x)
105 format (i2,' m',2x)
110 format (f4.1,' m')

  !--------- initialize static diagnostic fields --------------------

  id_land_mask = &
       register_static_field ( mod_name, 'land_mask', atmos_axes,  &
       'fractional amount of land', 'none', &
       range=frange )
  
  !--------- initialize diagnostic fields --------------------

  id_ice_mask = &
       register_diag_field ( mod_name, 'ice_mask', atmos_axes, Time, &
       'fractional amount of sea ice', 'none',  &
       range=frange )
  
  id_wind = &
       register_diag_field ( mod_name, 'wind', atmos_axes, Time, &
       'wind speed for flux calculations', 'm/s', &
       range=(/0.,vrange(2)/) )
  
  id_drag_moist = &
       register_diag_field ( mod_name, 'drag_moist', atmos_axes, Time, &
       'drag coeff for moisture',    'none'     )
  
  id_drag_heat  = &
       register_diag_field ( mod_name, 'drag_heat', atmos_axes, Time, &
       'drag coeff for heat',    'none'     )
  
  id_drag_mom   = &
       register_diag_field ( mod_name, 'drag_mom',  atmos_axes, Time, &
       'drag coeff for momentum',     'none'     )
  
  id_rough_moist = &
       register_diag_field ( mod_name, 'rough_moist', atmos_axes, Time, &
       'surface roughness for moisture',  'm'  )

  id_rough_heat = &
       register_diag_field ( mod_name, 'rough_heat', atmos_axes, Time, &
       'surface roughness for heat',  'm'  )

  id_rough_mom  = &
       register_diag_field ( mod_name, 'rough_mom',  atmos_axes, Time, &
       'surface roughness for momentum',  'm'  )

  id_u_star     = &
       register_diag_field ( mod_name, 'u_star',     atmos_axes, Time, &
       'friction velocity',   'm/s'   )

  id_b_star     = &
       register_diag_field ( mod_name, 'b_star',     atmos_axes, Time, &
       'buoyancy scale',      'm/s2'   )

  id_q_star     = &
       register_diag_field ( mod_name, 'q_star',     atmos_axes, Time, &
       'moisture scale',      'kg water/kg air'   )

  id_u_flux     = &
       register_diag_field ( mod_name, 'tau_x',      atmos_axes, Time, &
       'zonal wind stress',     'pa'   )

  id_v_flux     = &
       register_diag_field ( mod_name, 'tau_y',      atmos_axes, Time, &
       'meridional wind stress',     'pa'   )

  id_t_surf     = &
       register_diag_field ( mod_name, 't_surf',     atmos_axes, Time, &
       'surface temperature',    'deg_k', &
       range=trange    )

  ! + slm, Mar 25, 2002 -- add diagnositcs for t_ca, q_ca, and q_atm
  id_t_ca       = &
       register_diag_field ( mod_name, 't_ca',     atmos_axes, Time, &
       'canopy air temperature',    'deg_k', &
       range=trange    )

  id_q_atm     = &
       register_diag_field ( mod_name, 'q_atm',     atmos_axes, Time, &
       'specific humidity at btm level',    'kg/kg')

  id_q_surf     = &
       register_diag_field ( mod_name, 'q_surf',     atmos_axes, Time, &
       'surface specific humidity',    'kg/kg')
  ! - slm, Mar 25, 2002
  id_z_atm      = &
       register_diag_field ( mod_name, 'z_atm',     atmos_axes, Time, &
       'height of btm level',    'm')

  id_p_atm      = &
       register_diag_field ( mod_name, 'p_atm',     atmos_axes, Time, &
       'pressure at btm level',    'pa')

  id_gust       = &
       register_diag_field ( mod_name, 'gust',     atmos_axes, Time, &
       'gust scale',    'm/s')

  id_t_flux     = &
       register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
       'sensible heat flux',     'w/m2'    )

  id_q_flux     = &
       register_diag_field ( mod_name, 'evap',       atmos_axes, Time, &
       'evaporation rate',        'kg/m2/s'  )

  id_r_flux     = &
       register_diag_field ( mod_name, 'lwflx',      atmos_axes, Time, &
       'net (down-up) longwave flux',   'w/m2'    )

  id_t_atm      = &
       register_diag_field ( mod_name, 't_atm',      atmos_axes, Time, &
       'temperature at btm level',    'deg_k', &
       range=trange     )

  id_u_atm      = &
       register_diag_field ( mod_name, 'u_atm',      atmos_axes, Time, &
       'u wind component at btm level',  'm/s', &
       range=vrange    )

  id_v_atm      = &
       register_diag_field ( mod_name, 'v_atm',      atmos_axes, Time, &
       'v wind component at btm level',  'm/s', &
       range=vrange    )

  id_t_ref      = &
       register_diag_field ( mod_name, 't_ref',      atmos_axes, Time, &
       'temperature at '//label_zh, 'deg_k' , &
       range=trange      )

  id_rh_ref     = &
       register_diag_field ( mod_name, 'rh_ref',     atmos_axes, Time,   &
       'relative humidity at '//label_zh, 'percent' )

  id_u_ref      = &
       register_diag_field ( mod_name, 'u_ref',      atmos_axes, Time, &
       'zonal wind component at '//label_zm,  'm/s', &
       range=vrange )

  id_v_ref      = &
       register_diag_field ( mod_name, 'v_ref',      atmos_axes, Time,     &
       'meridional wind component at '//label_zm, 'm/s', &
       range=vrange )

  id_del_h      = &
       register_diag_field ( mod_name, 'del_h',      atmos_axes, Time,  &
       'ref height interp factor for heat', 'none' )
  id_del_m      = &
       register_diag_field ( mod_name, 'del_m',      atmos_axes, Time,     &
       'ref height interp factor for momentum','none' )
  id_del_q      = &
       register_diag_field ( mod_name, 'del_q',      atmos_axes, Time,     &
       'ref height interp factor for moisture','none' )

  ! + slm Jun 02, 2002 -- diagnostics of reference values over the land
  id_t_ref_land = &
       register_diag_field ( mod_name, 't_ref_land', Land_axes, Time, &
       'temperature at '//trim(label_zh)//' over land', 'deg_k' , &
       range=trange, missing_value =  -100.0)
  id_rh_ref_land= &
       register_diag_field ( mod_name, 'rh_ref_land', Land_axes, Time,   &
       'relative humidity at '//trim(label_zh)//' over land', 'percent',       &
       missing_value=-999.0)
  id_u_ref_land = &
       register_diag_field ( mod_name, 'u_ref_land',  Land_axes, Time, &
       'zonal wind component at '//trim(label_zm)//' over land',  'm/s', &
       range=vrange, missing_value=-999.0 )
  id_v_ref_land = &
       register_diag_field ( mod_name, 'v_ref_land',  Land_axes, Time,     &
       'meridional wind component at '//trim(label_zm)//' over land', 'm/s', &
       range=vrange, missing_value = -999.0 )
  ! - slm Jun 02, 2002
  id_q_ref = &
       register_diag_field ( mod_name, 'q_ref', atmos_axes, Time,     &
       'specific humidity at '//trim(label_zh), 'kg/kg', missing_value=-1.0)
  id_q_ref_land = &
       register_diag_field ( mod_name, 'q_ref_land', Land_axes, Time, &
       'specific humidity at '//trim(label_zh)//' over land', 'kg/kg',          &
       missing_value=-1.0)

  id_rough_scale = &
       register_diag_field ( mod_name, 'rough_scale', atmos_axes, Time, &
       'topographic scaling factor for momentum drag','1' )
!-----------------------------------------------------------------------

  end subroutine diag_field_init

! <DIAGFIELDS>
!   <NETCDF NAME="land_mask" UNITS="none">
!     fractional amount of land
!   </NETCDF>
!   <NETCDF NAME="wind" UNITS="m/s">
!     wind speed for flux calculations
!   </NETCDF>
!   <NETCDF NAME="drag_moist" UNITS="none">
!     drag coeff for moisture
!   </NETCDF>
!   <NETCDF NAME="drag_heat" UNITS="none">
!     drag coeff for heat
!   </NETCDF>
!   <NETCDF NAME="drag_mom" UNITS="none">
!     drag coeff for momentum
!   </NETCDF>
!   <NETCDF NAME="rough_moist" UNITS="m">
!     surface roughness for moisture
!   </NETCDF>
!   <NETCDF NAME="rough_heat" UNITS="m">
!     surface roughness for heat
!   </NETCDF>
!   <NETCDF NAME="rough_mom" UNITS="m">
!     surface roughness for momentum
!   </NETCDF>
!   <NETCDF NAME="u_star" UNITS="m/s">
!     friction velocity
!   </NETCDF>
!   <NETCDF NAME="b_star" UNITS="m/s">
!     buoyancy scale
!   </NETCDF>
!   <NETCDF NAME="q_star" UNITS="kg water/kg air">
!     moisture scale
!   </NETCDF>
!   <NETCDF NAME="t_atm" UNITS="deg_k">
!     temperature at btm level
!   </NETCDF>
!   <NETCDF NAME="u_atm" UNITS="m/s">
!     u wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="v_atm" UNITS="m/s">
!     v wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="q_atm" UNITS="kg/kg">
!     specific humidity at btm level
!   </NETCDF>
!   <NETCDF NAME="p_atm" UNITS="pa">
!     pressure at btm level
!   </NETCDF>
!   <NETCDF NAME="z_atm" UNITS="m">
!     height of btm level
!   </NETCDF>
!   <NETCDF NAME="gust" UNITS="m/s">
!     gust scale 
!   </NETCDF>
!   <NETCDF NAME="rh_ref" UNITS="percent">
!     relative humidity at ref height
!   </NETCDF>
!   <NETCDF NAME="t_ref" UNITS="deg_k">
!    temperature at ref height
!   </NETCDF>
!   <NETCDF NAME="u_ref" UNITS="m/s">
!    zonal wind component at ref height
!   </NETCDF>
!   <NETCDF NAME="v_ref" UNITS="m/s">
!    meridional wind component at ref height 
!   </NETCDF>
!   <NETCDF NAME="del_h" UNITS="none">
!    ref height interp factor for heat 
!   </NETCDF>
!   <NETCDF NAME="del_m" UNITS="none">
!    ref height interp factor for momentum 
!   </NETCDF>
!   <NETCDF NAME="del_q" UNITS="none">
!    ref height interp factor for moisture
!   </NETCDF>
!   <NETCDF NAME="tau_x" UNITS="pa">
!    zonal wind stress
!   </NETCDF>
!   <NETCDF NAME="tau_y" UNITS="pa">
!    meridional wind stress
!   </NETCDF>
!   <NETCDF NAME="ice_mask" UNITS="none">
!    fractional amount of sea ice 
!   </NETCDF>
!   <NETCDF NAME="t_surf" UNITS="deg_k">
!     surface temperature
!   </NETCDF>
!   <NETCDF NAME="t_ca" UNITS="deg_k">
!     canopy air temperature
!   </NETCDF>
!   <NETCDF NAME="q_surf" UNITS="kg/kg">
!     surface specific humidity 
!   </NETCDF>
!   <NETCDF NAME="shflx" UNITS="w/m2">
!     sensible heat flux
!   </NETCDF>
!   <NETCDF NAME="evap" UNITS="kg/m2/s">
!     evaporation rate 
!   </NETCDF>
!   <NETCDF NAME="lwflx" UNITS="w/m2">
!    net (down-up) longwave flux 
!   </NETCDF>

! </DIAGFIELDS>

! <INFO>


!   <NOTE>
!   <PRE>
!
!  MAIN PROGRAM EXAMPLE
!  --------------------
!
!       DO slow time steps (ocean)
!
!           call flux_ocean_to_ice
!
!           call ICE_SLOW_UP
!
!
!           DO fast time steps (atmos)
!
!                call sfc_boundary_layer
!
!                call ATMOS_DOWN
!
!                call flux_down_from_atmos
!
!                call LAND_FAST
!
!                call ICE_FAST
!
!                call flux_up_to_atmos
!
!                call ATMOS_UP
!
!           END DO
!
!           call ICE_SLOW_DN
!
!           call flux_ice_to_ocean
!
!           call OCEAN
!
!      END DO
!
!   LAND_FAST and ICE_FAST must update the surface temperature
!
! =======================================================================
!
! REQUIRED VARIABLES IN DEFINED DATA TYPES FOR COMPONENT MODELS
! --------------------------------------------------------------
!
! type (atmos_boundary_data_type) :: Atm
! type (surf_diff_type) :: Atm%Surf_Diff
!
! real, dimension(:)
!
!    Atm%lon_bnd   longitude axis grid box boundaries in radians
!                  must be monotonic
!    Atm%lat_bnd   latitude axis grid box boundaries in radians
!                  must be monotonic
!
! real, dimension(:,:)
!
!    Atm%t_bot     temperature at lowest model level
!    Atm%q_bot     specific humidity at lowest model level
!    Atm%z_bot     height above the surface for the lowest model level (m)
!    Atm%p_bot     pressure at lowest model level (pa)
!    Atm%u_bot     zonal wind component at lowest model level (m/s)
!    Atm%v_bot     meridional wind component at lowest model level (m/s)
!    Atm%p_surf    surface pressure (pa)
!    Atm%gust      gustiness factor (m/s)
!    Atm%flux_sw   net shortwave flux at the surface
!    Atm%flux_lw   downward longwave flux at the surface
!    Atm%lprec     liquid precipitation (kg/m2)
!    Atm%fprec     water equivalent frozen precipitation (kg/m2)
!    Atm%coszen    cosine of the zenith angle
!
!   (the following five fields are gathered into a data type for convenience in passing
!   this information through the different levels of the atmospheric model --
!   these fields are rlated to the simultaneous implicit time steps in the
!   atmosphere and surface models -- they are described more fully in
!   flux_exchange.tech.ps and
!   in the documntation for vert_diff_mod
!
!
!    Atm%Surf_Diff%dtmass   = dt/mass where dt = atmospheric time step ((i+1) = (i-1) for leapfrog) (s)
!                           mass = mass per unit area of lowest atmosphehic layer  (Kg/m2))
!    Atm%Surf_Diff%delta_t  increment ((i+1) = (i-1) for leapfrog) in temperature of
!                           lowest atmospheric layer  (K)
!    Atm%Surf_Diff%delta_q  increment ((i+1) = (i-1) for leapfrog) in specific humidity of
!                           lowest atmospheric layer (nondimensional -- Kg/Kg)
!    Atm%Surf_Diff%dflux_t  derivative of implicit part of downward temperature flux at top of lowest
!                           atmospheric layer with respect to temperature
!                           of lowest atmospheric layer (Kg/(m2 s))
!    Atm%Surf_Diff%dflux_q  derivative of implicit part of downward moisture flux at top of lowest
!                           atmospheric layer with respect to specific humidity of
!                           of lowest atmospheric layer (Kg/(m2 s))
!
!
! integer, dimension(4)
!
!    Atm%axes      Axis identifiers returned by diag_axis_init for the
!                  atmospheric model axes: X, Y, Z_full, Z_half.
!
! -----------------------------------------------
!
! type (land_boundary_data_type) :: Land
!
! real, dimension(:)
!
!    Land%lon_bnd     longitude axis grid box boundaries in radians
!                     must be monotonic
!    Land%lat_bnd     latitude axis grid box boundaries in radians
!                     must be monotonic
!
! logical, dimension(:,:,:)
!
!    Land%mask        land/sea mask (true for land)
!    Land%glacier     glacier mask  (true for glacier)
!
! real, dimension(:,:,:)
!
!    Land%tile_size   fractional area of each tile (partition)
!
!    Land%t_surf      surface temperature (deg k)
!    Land%albedo      surface albedo (fraction)
!    Land%rough_mom   surface roughness for momentum (m)
!    Land%rough_heat  surface roughness for heat/moisture (m)
!    Land%stomatal    stomatal resistance
!    Land%snow        snow depth (water equivalent) (kg/m2)
!    Land%water       water depth of the uppermost bucket (kg/m2)
!    Land%max_water   maximum water depth allowed in the uppermost bucket (kg/m2)
!
! -----------------------------------------------
!
!
! type (ice_boundary_data_type) :: Ice
!
! real, dimension(:)
!
!    Ice%lon_bnd       longitude axis grid box boundaries for temperature points
!                      in radians (must be monotonic)
!    Ice%lat_bnd       latitude axis grid box boundaries for temperature points
!                      in radians (must be monotonic)
!    Ice%lon_bnd_uv    longitude axis grid box boundaries for momentum points
!                      in radians (must be monotonic)
!    Ice%lat_bnd_uv    latitude axis grid box boundaries for momentum points
!                      in radians (must be monotonic)
!
! logical, dimension(:,:,:)
!
!    Ice%mask          ocean/land mask for temperature points
!                        (true for ocean, with or without ice)
!    Ice%mask_uv       ocean/land mask for momentum points
!                        (true for ocean, with or without ice)
!    Ice%ice_mask      optional ice mask (true for ice)
!
! real, dimension(:,:,:)
!
!    Ice%part_size     fractional area of each partition of a temperature grid box
!    Ice%part_size_uv  fractional area of each partition of a momentum grid box
!
!    the following fields are located on the ice top grid
!
!    Ice%t_surf        surface temperature (deg k)
!    Ice%albedo        surface albedo (fraction)
!    Ice%rough_mom     surface roughness for momentum (m)
!    Ice%rough_heat    surface roughness for heat/moisture (m)
!    Ice%u_surf        zonal (ocean/ice) current at the surface (m/s)
!    Ice%v_surf        meridional (ocean/ice) current at the surface (m/s)
!
!    the following fields are located on the ice bottom grid
!
!    Ice%flux_u        zonal wind stress (Pa)
!    Ice%flux_v        meridional wind stress (Pa)
!    Ice%flux_t        sensible heat flux (w/m2)
!    Ice%flux_q        specific humidity flux (kg/m2/s)
!    Ice%flux_sw       net (down-up) shortwave flux (w/m2)
!    Ice%flux_lw       net (down-up) longwave flux (w/m2)
!    Ice%lprec         mass of liquid precipitation since last time step (Kg/m2)
!    Ice%fprec         mass of frozen precipitation since last time step (Kg/m2)
!    Ice%runoff        mass of runoff water since last time step (Kg/m2)
!
! -----------------------------------------------
!
! type (ocean_boundary_data_type) :: Ocean
!
! real, dimension(:)
!
!    Ocean%Data%lon_bnd      longitude axis grid box boundaries for temperature
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lat_bnd      latitude axis grid box boundaries for temperature
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lon_bnd_uv   longitude axis grid box boundaries for momentum
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lat_bnd_uv   latitude axis grid box boundaries for momentum
!                            points on the ocean DATA GRID (radians)
!
!    Ocean%Ocean%lon_bnd     longitude axis grid box boundaries for temperature
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lat_bnd     latitude axis grid box boundaries for temperature
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lon_bnd_uv  longitude axis grid box boundaries for momentum
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lat_bnd_uv  latitude axis grid box boundaries for momentum
!                            points on the ocean MODEL GRID (radians)
!
!      Note: The data values in all longitude and latitude grid box boundary
!            array must be monotonic.
!
! logical, dimension(:,:)
!
!    Ocean%Data%mask       ocean/land mask for temperature points on the ocean
!                          DATA GRID (true for ocean)
!    Ocean%Data%mask_uv    ocean/land mask for momentum points on the ocean
!                          DATA GRID (true for ocean)
!
!    Ocean%Ocean%mask      ocean/land mask for temperature points on the ocean
!                          MODEL GRID (true for ocean)
!    Ocean%Ocean%mask_uv   ocean/land mask for momentum points on the ocean
!                          MODEL GRID (true for ocean)
!
! real, dimension(:,:)
!
!    Ocean%t_surf_data  surface temperature on the ocean DATA GRID (deg k)
!
!    Ocean%t_surf       surface temperature on the ocean MODEL GRID (deg k)
!    Ocean%u_surf       zonal ocean current at the surface on the ocean
!                       MODEL GRID (m/s)
!    Ocean%v_surf       meridional ocean current at the surface on the
!                       ocean MODEL GRID (m/s)
!    Ocean%frazil       frazil at temperature points on the ocean MODEL GRID
!
!   </PRE>
!   </NOTE>
! </INFO>

end module flux_exchange_mod

