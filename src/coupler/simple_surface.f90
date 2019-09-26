
module simple_surface_mod

!use atmos_coupled_mod, only: atmos_boundary_data_type
use atmos_model_mod, only: atmos_data_type

use surface_flux_mod, only: surface_flux


use diag_integral_mod, only: diag_integral_field_init, &
                             sum_diag_integral_field

use           fms_mod, only: file_exist, open_namelist_file, check_nml_error, &
                             error_mesg, FATAL, mpp_pe, mpp_root_pe, &
                             open_file, close_file, read_data, write_data, &
                             write_version_number, stdlog, set_domain

use  diag_manager_mod, only: register_diag_field,  &
                             register_static_field, send_data

use  time_manager_mod, only: time_type,get_time

use      constants_mod, only: rdgas, rvgas, cp_air, hlv, hlf

use ocean_rough_mod, only: compute_ocean_roughness
! mj know about surface topography
use spectral_dynamics_mod,only: get_surf_geopotential
use topography_mod, only:       get_ocean_mask
! mj read SSTs
use interpolator_mod, only: interpolate_type,interpolator_init&
     &,CONSTANT,interpolator
!mj q-flux
use qflux_mod, only: qflux_init,qflux,warmpool
!mj local surface heating
use physics_driver_mod, only: do_local_heating
use local_heating_mod, only: horizontal_heating,ngauss,hamp,pcenter

implicit none
private

public :: simple_surface_init,   &
          compute_flux,          &
          update_simple_surface, &
          simple_surface_end

!-----------------------------------------------------------------------
character(len=128) :: version = '$Id: simple_surface.f90,v 1.1.2.3 2005/05/21 02:02:04 pjp Exp $'
character(len=128) :: tagname = '$Name:  $'

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=14), parameter :: mod_name = 'simple_surface'

integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,              &
           id_rough_moist, id_rough_heat, id_rough_mom,             &
           id_u_star, id_b_star, id_u_flux, id_v_flux, id_t_surf,   &
           id_t_flux, id_q_flux, id_o_flux, id_r_flux,              &
           id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
           id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref,               &
           id_del_h,  id_del_m,  id_del_q, id_albedo, id_entrop_evap, &
           id_entrop_shflx, id_entrop_lwflx,                        &
           id_heat,id_tdt_horiz !mj

logical :: first_static = .true.
logical :: do_init = .true.
logical :: do_surface_heating = .false.

!-----------------------------------------------------------------------

real ::   z_ref_heat      = 2.,       &
          z_ref_mom       = 10.,      &
          heat_capacity   = 4.e08,    &
          land_capacity   = -1.,      & !mj
          trop_capacity   = -1.,      & !mj
          trop_cap_limit  = 15.,      & !mj
          heat_cap_limit  = 60.,      & !mj
          zsurf_cap_limit = 10.,      & !mj
          np_cap_factor   =  1.,      & !mj
          const_roughness = 3.21e-05, &
          const_albedo    = 0.30,     &
          albedo_exp      = 2.,       & !mj
          albedo_cntr     = 45.,      & !mj
          albedo_wdth     = 10.,      & !mj
	  max_of          = 25.,      &
	  lonmax_of       = 180.,     &
	  latmax_of       = 0.,       &
          latwidth_of     = 15.,      &
	  lonwidth_of     = 90.,      &
	  higher_albedo    = 0.38,    &
	  lat_glacier      = 45.,     &
	  maxofmerid       = .5,      &
	  latmaxofmerid    = 25.,     &
	  Tm               = 265.,    &
	  deltaT           = 40.,     &
          qflux_amp        = 30.,     & !mj
          qflux_width      = 16.        !mj


integer :: surface_choice   = 1
integer :: roughness_choice = 1
integer :: albedo_choice    = 1 ! 1->constant, 2->NH or SH step, 3->N-S symmetric step, 4->profile with albedo_exp, 5->tanh with albedo_cntr,albedo_wdth, 6->sin2\phi between const_albedo and higher_albedo
logical :: do_oflx          = .false.
logical :: do_oflxmerid     = .false.
logical :: do_qflux         = .false. !mj
logical :: do_warmpool      = .false. !mj
logical :: do_read_sst      = .false. !mj
logical :: do_sc_sst        = .false. !mj
character(len=256) :: sst_file
character(len=256) :: land_option = 'none'
character(len=256) :: land_sea_mask_file = 'lmask'
real,dimension(10) :: slandlon=0,slandlat=0,elandlon=-1,elandlat=-1

namelist /simple_surface_nml/ z_ref_heat, z_ref_mom,             &
                              surface_choice,  heat_capacity,    &
                              land_capacity,trop_capacity,       & !mj
                              trop_cap_limit, heat_cap_limit,    & !mj
                              np_cap_factor, zsurf_cap_limit,    & !mj
                              roughness_choice, const_roughness, &
                              albedo_choice, const_albedo, do_oflx, &
			      max_of, lonmax_of, latmax_of, latwidth_of, &
			      lonwidth_of, higher_albedo, lat_glacier, &
			      do_oflxmerid, maxofmerid, latmaxofmerid, Tm, &
			      deltaT,                            &
                              do_qflux,do_warmpool,              &  !mj
                              do_read_sst,do_sc_sst,sst_file,    &  !mj
                              land_option,slandlon,slandlat,     &  !mj
                              elandlon,elandlat,land_sea_mask_file,&!mj
                              albedo_exp,albedo_cntr,albedo_wdth    !mj

!-----------------------------------------------------------------------

!---- allocatable module storage ------

  real, allocatable, dimension(:,:) :: e_t_n, f_t_delt_n, &
                                       e_q_n, f_q_delt_n

  real, allocatable, dimension(:,:) :: dhdt_surf, dedt_surf, dedq_surf, &
                                       drdt_surf, dhdt_atm, dedq_atm, &
                                       flux_t, flux_q, flux_lw

  real, allocatable, dimension(:,:) :: sst, flux_u, flux_v, flux_o
! mj know about topography
  real, allocatable, dimension(:,:) :: zsurf,land_sea_heat_capacity
!mj read sst and land sea mask from input file
  real, allocatable, dimension(:,:) :: land_sea_mask
  logical,allocatable,dimension(:,:):: lmask_navy
  type(interpolate_type),save :: sst_interp, lmask_interp

contains

!#######################################################################

 subroutine compute_flux     ( dt, Time, Atm, land_frac,            &
                               t_surf_atm, albedo, rough_mom,       &
                               flux_u_atm, flux_v_atm, dtaudv_atm,  &
                               u_star, b_star                       )

 real,                   intent(in)  :: dt
 type       (time_type), intent(in)  :: Time
 type (atmos_data_type), intent(in)  :: Atm
 real, dimension(:,:),   intent(out) :: albedo,    rough_mom,       &
                                        land_frac, dtaudv_atm,      &
                                        flux_u_atm, flux_v_atm,     &
                                        u_star, b_star

real, dimension(:,:),   intent(out) :: t_surf_atm

real, dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: &
       u_surf, v_surf, rough_heat, rough_moist,          &
       stomatal, snow, water, max_water,                 &
       q_star, q_surf, cd_q, cd_t, cd_m, wind, dtaudu_atm

logical, dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: &
       mask, glacier, seawater

logical :: used

integer :: j
real :: lat, pi

pi = 4.0*atan(1.)



!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('compute_flux',  &
                 'must call simple_surface_init first', FATAL)

!-----------------------------------------------------------------------
!------ allocate storage also needed in flux_up_to_atmos -----

   allocate (e_t_n      (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             e_q_n      (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             f_t_delt_n (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             f_q_delt_n (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             dhdt_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             dedt_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             dedq_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             drdt_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             dhdt_atm   (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             dedq_atm   (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             flux_t     (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             flux_q     (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
             flux_lw    (size(Atm%t_bot,1), size(Atm%t_bot,2))  )

   u_surf     = 0.0
   v_surf     = 0.0
   stomatal   = 0.0
   snow       = 0.0
   water      = 1.0
   max_water  = 1.0

   mask    = .true.
   glacier = .false.
   seawater= .false.

   if(roughness_choice == 1) then
     rough_mom   = const_roughness
     rough_heat  = const_roughness
     rough_moist = const_roughness
   elseif(roughness_choice == 2) then
!    call compute_ocean_roughness (mask, flux_u, flux_v,   & ! Fez
     call compute_ocean_roughness (mask, u_star, &           ! Inchon and beyond. Changes answers.
                         rough_mom, rough_heat, rough_moist)
   endif

   if(albedo_choice == 1) then
     albedo = const_albedo
   elseif(albedo_choice == 2) then
     do j = 1, size(Atm%t_bot,2)
       lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))*180/pi
       ! mj SH or NH only
       if ( lat_glacier .ge. 0. ) then
          if ( lat > lat_glacier ) then

             albedo(:,j) = higher_albedo

          else

             albedo(:,j) = const_albedo

          endif
       else
          if ( lat < lat_glacier ) then

             albedo(:,j) = higher_albedo

          else

             albedo(:,j) = const_albedo

          endif
       endif
     enddo
   elseif(albedo_choice == 3) then ! NH and SH albedo step
     do j = 1, size(Atm%t_bot,2)
       lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))*180/pi

       if ( abs(lat) > lat_glacier ) then

         albedo(:,j) = higher_albedo

       else

         albedo(:,j) = const_albedo

       endif

     enddo
!mj add symmetric higher_albedo - exponential increase from equator to pole
  elseif(albedo_choice == 4) then
     do j = 1, size(Atm%t_bot,2)
        lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))*180/pi
        lat = abs(lat)
        albedo(:,j) = const_albedo + (higher_albedo-const_albedo)*(lat/90.)**albedo_exp
     enddo
!mj add symmetric higher_albedo - tanh increase around albedo_cntr with width
! albedo_wdth
  elseif(albedo_choice .eq. 5) then
     do j = 1, size(Atm%t_bot,2)
        lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))*180/pi
        lat = abs(lat)
        albedo(:,j) = const_albedo + (higher_albedo-const_albedo)*&
             0.5*(1+tanh((lat-albedo_cntr)/albedo_wdth))
     enddo
!mj add symmetric higher albedo - sin2 increase from equator to pole
   elseif(albedo_choice .eq. 6) then
      do j = 1, size(Atm%t_bot,2)
         lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))
         albedo(:,j) = const_albedo + (higher_albedo-const_albedo)*&
              sin(lat)**2
      enddo
   endif


   cd_t = 0.0
   cd_m = 0.0
   cd_q = 0.0

   t_surf_atm = sst

!  call surface_flux (Atm%t_bot, Atm%q_bot, Atm%u_bot, Atm%v_bot,        & ! Fez
!                     Atm%p_bot, Atm%z_bot,                              &
!                     Atm%p_surf, t_surf_atm, u_surf, v_surf,            &
!                     rough_mom, rough_heat, rough_moist,                &
!                     Atm%gust, stomatal,                                &
!                     snow, water,  max_water,                           &
!                     flux_t, flux_q, flux_lw, flux_u, flux_v,           &
!                     cd_m,   cd_t, cd_q, wind,                          &
!                     u_star, b_star, q_star, q_surf,                    &
!                     dhdt_surf, dedt_surf,  drdt_surf,                  &
!                     dhdt_atm,  dedq_atm,   dtaudv_atm,                 &
!                     dt, mask, glacier)

   call surface_flux (Atm%t_bot, Atm%q_bot, Atm%u_bot, Atm%v_bot,        & ! Lima
                      Atm%p_bot, Atm%z_bot, Atm%p_surf, t_surf_atm,      &
                      t_surf_atm,                                        & ! Required argument, intent(in). t_surf_atm instead of Land%t_ca
                      q_surf, u_surf, v_surf,                            &
                      rough_mom, rough_heat, rough_moist,                &
                      rough_mom,                                         & ! Required argument, intent(in). rough_mom instead of Land%rough_scale
                      Atm%gust, flux_t, flux_q, flux_lw, flux_u, flux_v, &
                      cd_m, cd_t, cd_q, wind, u_star, b_star, q_star,    &
                      dhdt_surf, dedt_surf,                              &
                      dedq_surf,                                         & ! Required argument, intent(out), but not needed by this model.
                      drdt_surf, dhdt_atm, dedq_atm,                     &
                      dtaudu_atm,                                        & ! Required argument, intent(out), but not needed by this model.
                      dtaudv_atm, dt,                                    & ! Required argument, intent(in). Looks like it should be .false. everywhere.
                      .not.mask,                                         &
                      seawater,                                          & ! Required argument, intent(in). Looks like fudgefactor for salt water. Use .false.
                      mask )                                               ! Required argument, intent(in). Looks like it should be .true. everywhere.

! intent(out):: flux_t, flux_q, flux_lw, flux_u, flux_v, cd_m, cd_t, cd_q, wind, u_star, b_star, q_star,
! intent(out):: dhdt_surf, dedt_surf, dedq_surf, drdt_surf, dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm
! intent(inout) :: q_surf
! All others intent(in)

   flux_u_atm = flux_u
   flux_v_atm = flux_v

   land_frac = 0.0

!=======================================================================
!-------------------- diagnostics section ------------------------------

   if ( id_wind       > 0 ) used = send_data ( id_wind,       wind,               Time )
   if ( id_drag_moist > 0 ) used = send_data ( id_drag_moist, cd_q,               Time )
   if ( id_drag_heat  > 0 ) used = send_data ( id_drag_heat,  cd_t,               Time )
   if ( id_drag_mom   > 0 ) used = send_data ( id_drag_mom,   cd_m,               Time )
   if ( id_rough_heat > 0 ) used = send_data ( id_rough_heat, rough_heat,         Time )
   if ( id_rough_mom  > 0 ) used = send_data ( id_rough_mom,  rough_mom,          Time )
   if ( id_u_star     > 0 ) used = send_data ( id_u_star,     u_star,             Time )
   if ( id_b_star     > 0 ) used = send_data ( id_b_star,     b_star,             Time )
   if ( id_t_atm      > 0 ) used = send_data ( id_t_atm,      Atm%t_bot,          Time )
   if ( id_u_atm      > 0 ) used = send_data ( id_u_atm,      Atm%u_bot,          Time )
   if ( id_v_atm      > 0 ) used = send_data ( id_v_atm,      Atm%v_bot,          Time )
   if ( id_albedo     > 0 ) used = send_data ( id_albedo,     albedo,             Time )



!=======================================================================

 end subroutine compute_flux

!#######################################################################

subroutine update_simple_surface (dt, Time, Atm, dt_t_atm, dt_q_atm  )

real, intent(in) :: dt
type       (time_type), intent(in)  :: Time
type (atmos_data_type), intent(in)  :: Atm

real, dimension(:,:),   intent(out) :: dt_t_atm, dt_q_atm

real, dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: &
                     gamma, dtmass, delta_t, delta_q, dflux_t, dflux_q, &
                     flux, deriv, dt_t_surf, &
                     entrop_evap, entrop_shflx, entrop_lwflx
real, dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: &
                     lon2d,lat2d,horiz_heat


 real    :: cp_inv
 logical :: used
! mj input SST
 real,dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: sst_new
! mj shallower ocean in tropics, land-sea contrast
 !real :: lon,lat,pi,loc_cap
 real :: pi
 integer :: i,j

   pi = 4.*atan(1.)

   flux_lw = Atm%flux_lw - flux_lw

   dtmass  = Atm%Surf_Diff%dtmass
   delta_t = Atm%Surf_Diff%delta_t
   delta_q = Atm%Surf_Diff%delta_q
   dflux_t = Atm%Surf_Diff%dflux_t
   dflux_q = Atm%Surf_Diff%dflux_q

 cp_inv = 1.0/cp_air

 ! temperature

   gamma      =  1./ (1.0 - dtmass*(dflux_t + dhdt_atm*cp_inv))
   e_t_n      =  dtmass*dhdt_surf*cp_inv*gamma
   f_t_delt_n = (delta_t + dtmass * flux_t*cp_inv) * gamma

   flux_t     =  flux_t        + dhdt_atm * f_t_delt_n
   dhdt_surf  =  dhdt_surf     + dhdt_atm * e_t_n

! moisture

   gamma      =  1./ (1.0 - dtmass*(dflux_q + dedq_atm))
   e_q_n      =  dtmass*dedt_surf*gamma
   f_q_delt_n = (delta_q  + dtmass * flux_q) * gamma

   flux_q     =  flux_q        + dedq_atm * f_q_delt_n
   dedt_surf  =  dedt_surf     + dedq_atm * e_q_n


   dt_t_surf  = 0.0
   !###############################
   ! heating due to surface fluxes
   !
   if(surface_choice == 1)then
      if(do_sc_sst) then !mj sst read from input file
         call interpolator( sst_interp, Time, sst_new, trim(sst_file) )
         dt_t_surf = sst_new - sst
      else
      !else   !mj ocean depth function of latitude
!!$         land_sea_heat_capacity = heat_capacity
!!$         if ( trop_capacity .ne. heat_capacity .or. np_cap_factor .ne. 1.0 ) then
!!$            do j=1,size(Atm%t_bot,2)
!!$               lat = 0.5*180/pi*( Atm%lat_bnd(j+1) + Atm%lat_bnd(j) )
!!$               if ( lat > 0. ) then
!!$                  loc_cap = heat_capacity*np_cap_factor
!!$               else
!!$                  loc_cap = heat_capacity
!!$               endif
!!$               if ( abs(lat) < trop_cap_limit ) then
!!$                  land_sea_heat_capacity(:,j) = trop_capacity
!!$               elseif ( abs(lat) < heat_cap_limit ) then
!!$                  land_sea_heat_capacity(:,j) = trop_capacity*(1.-(abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)) + (abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)*loc_cap
!!$               elseif ( lat > heat_cap_limit ) then
!!$                  land_sea_heat_capacity(:,j) = loc_cap
!!$               end if
!!$            enddo
!!$         endif
!!$         if (trim(land_option) .eq. 'interpolated')then
!!$            lmask_navy = get_ocean_mask(Atm%lon_bnd,Atm%lat_bnd,lmask_navy)
!!$            where(lmask_navy) land_sea_heat_capacity = land_capacity
!!$! mj land heat capacity function of surface topography
!!$         else if(trim(land_option) .eq. 'zsurf')then
!!$            call get_surf_geopotential(zsurf)
!!$            where ( zsurf > zsurf_cap_limit ) land_sea_heat_capacity = land_capacity
!!$! mj land heat capacity given in inputfile
!!$         else if(trim(land_option) .eq. 'input')then
!!$            where(land_sea_mask .gt. 0) land_sea_heat_capacity = land_capacity
!!$! mj land heat capacity given through ?landlon, ?landlat
!!$         else if(trim(land_option) .eq. 'lonlat')then
!!$            do j=1,size(Atm%t_bot,2)
!!$               lat = 0.5*180/pi*( Atm%lat_bnd(j+1) + Atm%lat_bnd(j) )
!!$               do i=1,size(Atm%t_bot,1)
!!$                  lon = 0.5*180/pi*( Atm%lon_bnd(i+1) + Atm%lon_bnd(i) )
!!$                  do k=1,size(slandlat)
!!$                     if ( lon >= slandlon(k) .and. lon <= elandlon(k) &
!!$                          &.and. lat >= slandlat(k) .and. lat <= elandlat(k) )then
!!$                        land_sea_heat_capacity(i,j) = land_capacity
!!$                     endif
!!$                  enddo
!!$               enddo
!!$            enddo
!!$         endif

         flux    = (flux_lw + Atm%flux_sw - hlf*Atm%fprec &
              - (flux_t + hlv*flux_q) + flux_o)*dt/land_sea_heat_capacity

         deriv   = - (dhdt_surf + hlv*dedt_surf + drdt_surf)*dt/land_sea_heat_capacity
        !flux    = (flux_lw + Atm%flux_sw - hlf*Atm%fprec &
        !        - (flux_t + hlv*flux_q) + flux_o)*dt/heat_capacity

        !   deriv   = - (dhdt_surf + hlv*dedt_surf + drdt_surf)*dt/heat_capacity
! mj end

         dt_t_surf = flux/(1.0 -deriv)
      endif

   elseif(surface_choice == 2) then

      dt_t_surf  = 0.0

   endif
   
   !###############################
   ! additional heating
   !
   if ( do_surface_heating ) then
      do j = 1,size(Atm%t_bot,2)
         do i = 1,size(Atm%t_bot,1)
            lat2d(i,j) = 0.5*( Atm%lat_bnd(j+1) + Atm%lat_bnd(j) )
            lon2d(i,j) = 0.5*( Atm%lon_bnd(i+1) + Atm%lon_bnd(i) )
         enddo
      enddo
      call horizontal_heating(Time,lon2d,lat2d,horiz_heat)
   else
      horiz_heat = 0.0 ! for diagnostics output
   endif
   dt_t_surf = dt_t_surf + horiz_heat*dt
   !###############################
   ! apply all heating
   sst = sst + dt_t_surf


  flux_t     = flux_t      + dt_t_surf*dhdt_surf
  flux_q     = flux_q      + dt_t_surf*dedt_surf
  flux_lw    = flux_lw     - dt_t_surf*drdt_surf
  dt_t_atm   = f_t_delt_n  + dt_t_surf*e_t_n
  dt_q_atm   = f_q_delt_n  + dt_t_surf*e_q_n



!=======================================================================
!-------------------- diagnostics section ------------------------------

   if ( id_t_surf > 0 ) used = send_data ( id_t_surf, sst, Time )
   if ( id_t_flux > 0 ) used = send_data ( id_t_flux, flux_t, Time )
   if ( id_r_flux > 0 ) used = send_data ( id_r_flux, flux_lw, Time )
   if ( id_q_flux > 0 ) used = send_data ( id_q_flux, flux_q, Time )
   if ( id_o_flux > 0 ) used = send_data ( id_o_flux, flux_o, Time )
   if ( id_heat   > 0 ) used = send_data ( id_heat,land_sea_heat_capacity,Time )
   if ( id_tdt_horiz   > 0 ) used = send_data ( id_tdt_horiz,horiz_heat,Time )
   if ( id_entrop_evap > 0 ) then
      entrop_evap = flux_q/sst
      used = send_data ( id_entrop_evap, entrop_evap, Time)
   endif
   if (id_entrop_shflx > 0 ) then
      entrop_shflx = flux_t/sst
      used = send_data ( id_entrop_shflx, entrop_shflx, Time)
   endif
   if (id_entrop_lwflx > 0 ) then
      entrop_lwflx = flux_lw/sst
      used = send_data ( id_entrop_lwflx, entrop_lwflx, Time)
   endif

   call sum_diag_integral_field ('evap', flux_q*86400.)

!=======================================================================
!---- deallocate module storage ----


   deallocate (f_t_delt_n, f_q_delt_n, e_t_n, e_q_n)
   deallocate(dhdt_surf, dedt_surf, dedq_surf, drdt_surf, dhdt_atm, dedq_atm, &
              flux_t, flux_q, flux_lw)



!-----------------------------------------------------------------------

 end subroutine update_simple_surface

!#######################################################################

 subroutine simple_surface_init ( Time, Atm)

 type       (time_type), intent(in)  :: Time
 type (atmos_data_type), intent(in)  :: Atm

 integer :: unit, ierr, io

 integer :: i, j, k, lati
 real :: xx, xx2, lat, lon, pi, y0
 real :: coslat !mj
 real, dimension(100) :: oftabl
 real, dimension(32) :: ssttabl
! mj shallower ocean in tropics, land-sea contrast
 real :: loc_cap
 logical :: ocean_mask_worked

 pi = 4.0*atan(1.)

 !-----------------------------------------------------------------------
!------ read namelist ------

   if ( file_exist('input.nml')) then
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=simple_surface_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'simple_surface_nml')
      enddo
 10   call close_file (unit)
   endif
!mj make choices compatible
   !if(do_read_sst .or. do_sc_sst) call error_mesg ('simple_surface',  &
   !              'THERE IS A BUG WITH DO_READ_SST, SO I AM STOPPING', FATAL)
   if(do_sc_sst) do_read_sst = .true.
   if(trop_capacity .le. 0.) trop_capacity = heat_capacity
   if(land_capacity .le. 0.) land_capacity = heat_capacity

!--------- write version number and namelist ------------------

   call write_version_number ( version, tagname )
   if ( mpp_pe() == mpp_root_pe() ) then
        write (stdlog(), nml=simple_surface_nml)
   endif

   call diag_integral_field_init ('evap', 'f6.3')
   call diag_field_init ( Time, Atm%axes(1:2) )

allocate(sst   (size(Atm%t_bot,1),size(Atm%t_bot,2)))
allocate(flux_u(size(Atm%t_bot,1),size(Atm%t_bot,2)))
allocate(flux_v(size(Atm%t_bot,1),size(Atm%t_bot,2)))
allocate(flux_o(size(Atm%t_bot,1),size(Atm%t_bot,2)))

!mj read fixed SSTs
if( do_read_sst ) then
   call interpolator_init( sst_interp, trim(sst_file)//'.nc',Atm%lon_bnd,Atm%lat_bnd, data_out_of_bounds=(/CONSTANT/) )
endif

!mj land sea mask
if (surface_choice .eq. 1 .and. .not. do_sc_sst)then
    allocate(land_sea_heat_capacity(size(Atm%t_bot,1), size(Atm%t_bot,2)))
    land_sea_heat_capacity = heat_capacity
    if ( trop_capacity .ne. heat_capacity .or. np_cap_factor .ne. 1.0 ) then
        do j=1,size(Atm%t_bot,2)
           lat = 0.5*180/pi*( Atm%lat_bnd(j+1) + Atm%lat_bnd(j) )
           if ( lat > 0. ) then
              loc_cap = heat_capacity*np_cap_factor
           else
              loc_cap = heat_capacity
           endif
           if ( abs(lat) < trop_cap_limit ) then
              land_sea_heat_capacity(:,j) = trop_capacity
           elseif ( abs(lat) < heat_cap_limit ) then
              land_sea_heat_capacity(:,j) = trop_capacity*(1.-(abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)) + (abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)*loc_cap
           elseif ( lat > heat_cap_limit ) then
              land_sea_heat_capacity(:,j) = loc_cap
           end if
        enddo
     endif
     if( trim(land_option) .eq. 'input' ) then
        allocate(land_sea_mask(size(Atm%t_bot,1),size(Atm%t_bot,2)))
        if(mpp_pe() .eq. mpp_root_pe()) write(*,'(a)') 'Reading land-sea mask from file INPUT/'//trim(land_sea_mask_file)//'.nc'
        call read_data('INPUT/'//trim(land_sea_mask_file),trim(land_sea_mask_file),land_sea_mask,domain=Atm%domain)
! mj use navy land-sea mask
     else if (trim(land_option) .eq. 'interpolated')then
        allocate(lmask_navy(size(land_sea_heat_capacity,1),size(land_sea_heat_capacity,2)))
        ocean_mask_worked = get_ocean_mask(Atm%lon_bnd,Atm%lat_bnd,lmask_navy)
        if(.not.ocean_mask_worked) then
           call error_mesg('get_ocean_mask','land_option="'//trim(land_option)//'"'// &
                         ' and ocean_mask is not present but water data file does not exist', FATAL)
        endif
        where(.not. lmask_navy) land_sea_heat_capacity = land_capacity
! mj land heat capacity function of surface topography
     else if(trim(land_option) .eq. 'zsurf')then
        allocate(zsurf(size(Atm%t_bot,1), size(Atm%t_bot,2)))
        call get_surf_geopotential(zsurf)
        where ( zsurf > zsurf_cap_limit ) land_sea_heat_capacity = land_capacity
! mj land heat capacity given in inputfile
     else if(trim(land_option) .eq. 'input')then
        where(land_sea_mask .gt. 0) land_sea_heat_capacity = land_capacity
! mj land heat capacity given through ?landlon, ?landlat
     else if(trim(land_option) .eq. 'lonlat')then
        do j=1,size(Atm%t_bot,2)
           lat = 0.5*180/pi*( Atm%lat_bnd(j+1) + Atm%lat_bnd(j) )
           do i=1,size(Atm%t_bot,1)
              lon = 0.5*180/pi*( Atm%lon_bnd(i+1) + Atm%lon_bnd(i) )
              do k=1,size(slandlat)
                 if ( lon >= slandlon(k) .and. lon <= elandlon(k) &
                      &.and. lat >= slandlat(k) .and. lat <= elandlat(k) )then
                    land_sea_heat_capacity(i,j) = land_capacity
                 endif
              enddo
           enddo
        enddo
     endif
endif


if(file_exist('INPUT/simple_surface.res.nc')) then
  call read_data('INPUT/simple_surface.res.nc', 'sst',    sst,    domain=Atm%domain)
  call read_data('INPUT/simple_surface.res.nc', 'flux_u', flux_u, domain=Atm%domain)
  call read_data('INPUT/simple_surface.res.nc', 'flux_v', flux_v, domain=Atm%domain)
else if(file_exist('INPUT/simple_surface.res')) then
  unit = open_file(file='INPUT/simple_surface.res',form='unformatted',&
                   action='read')
  call set_domain(Atm%domain)

  call read_data(unit, sst)
  call read_data(unit, flux_u)
  call read_data(unit, flux_v)

  call close_file(unit)
!mj read fixed SSTs
else if( do_read_sst ) then
   call interpolator( sst_interp, Time, sst, trim(sst_file) )
else
  do j = 1, size(Atm%t_bot,2)
    lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))
!    xx = 1. - sin(1.5*lat)*sin(1.5*lat)
!    if(abs(lat) .gt. atan(1.0)*60.0/45.0) xx = 0.0
!    sst(:,j) = 273.15 + 27.0*xx
    xx = sin(lat)*sin(lat)
    xx2 = xx*xx
!    sst(:,j) = 305.0 - 10.0*xx2 - 30.0*xx
!    sst(:,j) = 290.0
!    sst(:,j) = 300.0 - 35.0*xx2
! CL functional form
    sst(:,j) = Tm - deltaT*(3.*xx-1.)/3.
! mj equal temperature
   if (Tm.le.0.) then
     sst(:,j) = -Tm
   endif
! APE "control+2" experiment
   if (Tm.ge.399.) then
     sst(:,j) = 273.15+29.0*(1.0 - sin(1.40369*lat)*sin(1.40369*lat))
     if (abs(lat).ge.pi/3.) sst(:,j) = 273.15
   end if
! APE "control" experiment
    if (Tm.ge.499.) then
      sst(:,j) = 273.15+27.0*(1.0 - sin(1.5*lat)*sin(1.5*lat))
      if (abs(lat).ge.pi/3.) sst(:,j) = 273.15
    endif
! APE "qobs" experiment
    if (Tm.ge.599.) then
     sst(:,j)=273.15+.5*27.*(1. - sin(1.5*lat)*sin(1.5*lat))+ .5*27.*(1.-sin(1.5*lat)*sin(1.5*lat)*sin(1.5*lat)*sin(1.5*lat))
       if (abs(lat).ge.pi/3.) sst(:,j) = 273.15
    end if
! APE "qobs+5" experiment
    if (Tm.ge.699.) then
       sst(:,j)=sst(:,j) + 5.
    end if
! Jian's qobs + El Nino forcing
    if (Tm.ge.704.) then
      sst(:,j)=sst(:,j)- 5.+2.*(1.-sin(6.* max(min(abs(lat),pi/12.),0.))**4)
    end if
! Jian's qobs + expanding 1
    if (Tm.ge.709.) then
      sst(:,j) = sst(:,j) - 2.*(1.-sin(6.* max(min(abs(lat),pi/12.),0.))**4)
      y0 = 0.
      sst(:,j) = sst(:,j) + 2.*(1.-sin(3.* max(min((abs(lat)-y0),pi/6.),0.))**4)
    end if
! Jian's qobs + expanding 2
    if (Tm.ge.714.) then
      sst(:,j) = sst(:,j) - 2.*(1.-sin(3.* max(min((abs(lat)-y0),pi/6.),0.))**4)
      y0 = 20.*pi/180.
      sst(:,j) = sst(:,j) + 2.*(1.-sin(3.* max(min((abs(lat)-y0),pi/6.),0.))**4)
    end if
  end do
  if (Tm.ge.800.) then
    data ssttabl / 244.2464, 244.9859, 246.1535, 247.7335, 249.7301, 252.0308, 254.5561,&
       257.2648, 260.1125, 263.0303, 265.9484, 268.7581, 271.4256, 273.9761,&
       276.4128, 278.7544, 281.0249, 283.2028, 285.2909, 287.3179, 289.2411,&
       291.0211, 292.6213, 293.9812, 295.1642, 296.2877, 297.2654, 298.1137,&
       298.8906, 299.5288, 300.2687, 301.1037/
    do j=1, size(Atm%t_bot,2)
      lat = 0.5*180./pi*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))
      sst(:,j) = ssttabl(32)
      if(abs(lat).ge.2.) then
        sst(:,j) = ssttabl(31)
      end if
      if(abs(lat).ge.5.) then
        sst(:,j) = ssttabl(30)
      end if
      if(abs(lat).ge.8.) then
        sst(:,j) = ssttabl(29)
      end if
      if(abs(lat).ge.11.) then
        sst(:,j) = ssttabl(28)
      end if
      if(abs(lat).ge.14.) then
        sst(:,j) = ssttabl(27)
      end if
      if(abs(lat).ge.17.) then
        sst(:,j) = ssttabl(26)
      end if
      if(abs(lat).ge.19.) then
        sst(:,j) = ssttabl(25)
      end if
      if(abs(lat).ge.22.) then
        sst(:,j) = ssttabl(24)
      end if
      if(abs(lat).ge.25.) then
        sst(:,j) = ssttabl(23)
      end if
      if(abs(lat).ge.28.) then
        sst(:,j) = ssttabl(22)
      end if
      if(abs(lat).ge.31.) then
        sst(:,j) = ssttabl(21)
      end if
      if(abs(lat).ge.33.) then
        sst(:,j) = ssttabl(20)
      end if
      if(abs(lat).ge.36.) then
        sst(:,j) = ssttabl(19)
      end if
      if(abs(lat).ge.39.) then
        sst(:,j) = ssttabl(18)
      end if
      if(abs(lat).ge.42.) then
        sst(:,j) = ssttabl(17)
      end if
      if(abs(lat).ge.45.) then
        sst(:,j) = ssttabl(16)
      end if
      if(abs(lat).ge.47.) then
        sst(:,j) = ssttabl(15)
      end if
      if(abs(lat).ge.50.) then
        sst(:,j) = ssttabl(14)
      end if
      if(abs(lat).ge.53.) then
        sst(:,j) = ssttabl(13)
      end if
      if(abs(lat).ge.56.) then
        sst(:,j) = ssttabl(12)
      end if
      if(abs(lat).ge.58.) then
        sst(:,j) = ssttabl(11)
      end if
      if(abs(lat).ge.61.) then
        sst(:,j) = ssttabl(10)
      end if
      if(abs(lat).ge.64.) then
        sst(:,j) = ssttabl(9)
      end if
      if(abs(lat).ge.67.) then
        sst(:,j) = ssttabl(8)
      end if
      if(abs(lat).ge.70.) then
        sst(:,j) = ssttabl(7)
      end if
      if(abs(lat).ge.72.) then
        sst(:,j) = ssttabl(6)
      end if
      if(abs(lat).ge.75.) then
        sst(:,j) = ssttabl(5)
      end if
      if(abs(lat).ge.78.) then
        sst(:,j) = ssttabl(4)
      end if
      if(abs(lat).ge.81.) then
        sst(:,j) = ssttabl(3)
      end if
      if(abs(lat).ge.84.) then
        sst(:,j) = ssttabl(2)
      end if
      if(abs(lat).ge.86.) then
        sst(:,j) = ssttabl(1)
      end if
    enddo
  end if
  flux_u = 0.0
  flux_v = 0.0
endif

if (do_oflx) then
  do i=1, size(Atm%t_bot,1)
    do j=1, size(Atm%t_bot,2)
      lat = 0.5*180./pi*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))
      lon = 0.5*180./pi*(Atm%lon_bnd(i+1) + Atm%lon_bnd(i))
      if (lonmax_of - lonwidth_of < lon .and. lon < lonmax_of + lonwidth_of ) then
        flux_o(i,j) = max_of*sin(pi*(lon - (lonmax_of - lonwidth_of))/lonwidth_of)* &
           exp(-(lat-latmax_of)**2./latwidth_of**2.)
      endif
    enddo
  enddo
else
   flux_o = 0.
endif

if(do_oflxmerid) then
     data oftabl / -36.3215, -36.2461, -36.1313, -35.9763, -35.7805, &
     -35.5430, -35.2629, -34.9392, -34.5707, -34.1562, -33.6940, -33.1824, &
     -32.6194, -32.0029, -31.3302, -30.5986, -29.8049, -28.9454, -28.0163, &
     -27.0132, -25.9313, -24.7655, -23.5104, -22.1599, -20.7081, -19.1486, &
     -17.4754, -15.6824, -13.7641, -11.7158,  -9.5341,  -7.2172,  -4.7656, &
     -2.1827,   0.5248,   3.3459,   6.2647,   9.2603,  12.3059,  15.3690, &
     18.4114,  21.3899,  24.2574,  26.9641,  29.4597,  31.6952,  33.6258, &
     35.2132,  36.4275,  37.2494,  37.6712,  37.6971,  37.3434,  36.6366, &
     35.6125,  34.3136,  32.7867,  31.0806,  29.2438,  27.3228,  25.3601, &
     23.3935,  21.4554,  19.5724,  17.7656,  16.0508,  14.4388,  12.9365, &
     11.5470,  10.2706,   9.1053,   8.0474,   7.0916,   6.2320,   5.4622, &
     4.7753,   4.1646,   3.6232,   3.1449,   2.7234,   2.3529,   2.0280, &
     1.7439,   1.4959,   1.2800,   1.0923,   0.9296,   0.7889,   0.6673, &
     0.5626,   0.4726,   0.3954,   0.3294,   0.2730,   0.2249,   0.1842, &
     0.1497,   0.1205,   0.0960,   0.0425/
   do j=1, size(Atm%t_bot,2)
      lat = 0.5*180./pi*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))
      if(abs(lat).le.latmaxofmerid) then
         lati = floor(4.*abs(lat))
         flux_o(:,j) = flux_o(:,j)+maxofmerid*(oftabl(lati)+(oftabl(lati+1)-oftabl(lati))*(4.*abs(lat)-lati))
      end if
   enddo
endif

if ( do_qflux .or. do_warmpool) then
   call qflux_init
!mj q-flux as in Merlis et al (2013) [Part II]
   if ( do_qflux ) call qflux(Atm%lat_bnd,flux_o)
!mj q-flux to create a tropical temperature perturbation
   if ( do_warmpool) call warmpool(Atm%lon_bnd,Atm%lat_bnd,flux_o)
endif

!mj adding local surface heating
if ( do_local_heating ) then
   do j=1,ngauss
      if ( hamp(j) .ne. 0. .and. pcenter(j) .lt. 0. ) then
         do_surface_heating = .true.
      endif
   enddo
endif


    do_init = .false.

!-----------------------------------------------------------------------

 end subroutine simple_surface_init

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)

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

   id_t_flux     = &
   register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
                        'sensible heat flux',     'w/m2'    )

   id_q_flux     = &
   register_diag_field ( mod_name, 'evap',       atmos_axes, Time, &
                        'evaporation rate',        'kg/m2/s'  )

   id_o_flux     = &
   register_diag_field (mod_name, 'oflx',        atmos_axes, Time, &
                        'prescribed ocean heat divergence', 'w/m2' )

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
   id_albedo      = &
   register_diag_field ( mod_name, 'albedo',      atmos_axes, Time,     &
                        'surface albedo','none' )
   id_heat        = & !mj
   register_diag_field ( mod_name, 'heat_capacity',atmos_axes, Time,     &
                        'mixed layer heat capacity','none' )
   id_tdt_horiz        = & !mj
   register_diag_field ( mod_name, 'tdt_surf',atmos_axes, Time,     &
                        'local surface heating','none' )
   id_entrop_evap      = &
   register_diag_field ( mod_name, 'entrop_evap', atmos_axes, Time,     &
                        'entropy source from evap','kg/m2/s/K' )

   id_entrop_shflx      = &
   register_diag_field ( mod_name, 'entrop_shflx', atmos_axes, Time,     &
                        'entropy source from SH flux','w/m2/K' )

   id_entrop_lwflx      = &
   register_diag_field ( mod_name, 'entrop_lwflx', atmos_axes, Time,     &
                        'entropy source from LW flux','w/m2/K' )




!-----------------------------------------------------------------------

end subroutine diag_field_init


!########################################################################

subroutine simple_surface_end (Atm)

type (atmos_data_type), intent(in)  :: Atm
integer :: unit

call write_data('RESTART/simple_surface.res.nc', 'sst',    sst,    domain=Atm%domain)
call write_data('RESTART/simple_surface.res.nc', 'flux_u', flux_u, domain=Atm%domain)
call write_data('RESTART/simple_surface.res.nc', 'flux_v', flux_v, domain=Atm%domain)

!unit = open_file(file='RESTART/simple_surface.res',form='unformatted', &
!                       action='write')
!call set_domain(Atm%domain)

!call write_data(unit, sst)
!call write_data(unit, flux_u)
!call write_data(unit, flux_v)

!call close_file(unit)

end subroutine simple_surface_end

!#######################################################################

end module simple_surface_mod
