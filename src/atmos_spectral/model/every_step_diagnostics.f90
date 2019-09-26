module every_step_diagnostics_mod

use              fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

use     time_manager_mod, only: time_type

use    field_manager_mod, only: MODEL_ATMOS

use   tracer_manager_mod, only: get_number_tracers, get_tracer_names, get_tracer_index, NO_TRACER

use     diag_manager_mod, only: diag_axis_init, register_diag_field, register_static_field, send_data

use press_and_geopot_mod, only: pressure_variables

use       transforms_mod, only: get_deg_lon, get_deg_lat, get_grid_domain
use          spec_mpp_mod,only: grid_domain

use        constants_mod, only: hlv, cp_air, grav, kappa, rvgas, rdgas, &
                                RADIUS, RADIAN, OMEGA !mj

!===============================================================================================
implicit none
private
!===============================================================================================

public :: every_step_diagnostics_init, every_step_diagnostics, every_step_diagnostics_end

!===============================================================================================

character(len=128), parameter :: version = '$Id: every_step_diagnostics.f90,v 11.0.2.1 2005/05/13 18:16:38 pjp Exp $'
character(len=128), parameter :: tagname = '$Name:  $'

real, parameter :: p00=1.e5, virtual_factor=(rvgas/rdgas) - 1.0

integer :: id_ps, id_u, id_v, id_t, num_levels, num_tracers

integer :: id_wf, id_z, id_drystaten, id_moiststaten, id_uv, id_vq, id_vdse, id_pt,    &
    id_vp, id_wu, id_wv, id_wq, id_wdse, id_wup, id_wvp, id_wqp, id_wdsep, id_wp,      &
    id_wsubt, id_wsubtv, id_vqint, id_udt_damp, id_vdt_damp, id_tdt_damp,              &
    id_entrop_dampuv, id_entrop_dampt, id_tdt_dampuv, id_tdt_tempcor, id_qdt_watercor, &
    id_entrop_tempcor, id_kegen, id_kegenq, id_kegenqtinv
!mj
integer :: id_upvp, id_upwp, id_vpTp, id_psi_star, id_psi_dwc
!mj

integer, allocatable, dimension(:) :: id_tr, id_dt_hadv, id_dt_vadv

logical :: module_is_initialized = .false.

integer :: two_dt_id_ps, two_dt_id_u, two_dt_id_v, two_dt_id_t
integer, allocatable, dimension(:) :: two_dt_id_tr
integer :: iwt, num_time_steps

real, allocatable, dimension(:,:)     :: two_dt_ps
real, allocatable, dimension(:,:,:)   :: two_dt_u, two_dt_v, two_dt_t
real, allocatable, dimension(:,:,:,:) :: two_dt_tr

!mj moved from every_step_diagnostics_init for use in every_step_diagnostics as well
real, allocatable, dimension(:) :: lon, lat
!jm

integer :: is, ie, js, je
integer :: nsphum=NO_TRACER ! tracer index of specific humidity. Initialized to NO_TRACER, then set in subroutine every_step_diagnostics_init.

type(time_type) :: Time_save ! every_step_diagnostics_end needs time because it is a required argument of send_data,
                             ! even though the fields are static. When Time is not available to every_step_diagnostics_end
                             ! it uses Time_save instead.
!===============================================================================================

contains

!===============================================================================================

subroutine every_step_diagnostics_init(Time, lon_max, lat_max, num_levels_in, reference_sea_level_press)

type(time_type), intent(in) :: Time
integer, intent(in) :: lon_max, lat_max, num_levels_in
real, intent(in) :: reference_sea_level_press

integer, dimension(2) :: axes_2d
integer, dimension(2) :: axes_zm !mj
integer, dimension(3) :: axes_3d
integer :: id_lon, id_lat, id_pfull, ntr
real, dimension(2) :: vrange,trange
character(len=128) :: tname, longname, units
character(len=16) :: mod_name = 'dynamics_every'
real, dimension(num_levels_in)   :: p_full, ln_p_full
real, dimension(num_levels_in+1) :: p_half, ln_p_half

call write_version_number(version, tagname)
!mj
allocate(lon(lon_max),lat(lat_max))
call get_deg_lon(lon)
call get_deg_lat(lat)
!jm
id_lon = diag_axis_init('lon_every', lon, 'degrees_E', 'x', 'longitude', set_name=mod_name, Domain2=grid_domain)
id_lat = diag_axis_init('lat_every', lat, 'degrees_N', 'y', 'latitude',  set_name=mod_name, Domain2=grid_domain)
call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, reference_sea_level_press)
p_full = .01*p_full
id_pfull = diag_axis_init('pfull_every', p_full, 'hPa', 'z', 'approx full pressure level', direction=-1, set_name=mod_name)

vrange = (/ -400., 400. /)
trange = (/  100., 400. /)
axes_2d = (/ id_lon, id_lat /)
axes_3d = (/ id_lon, id_lat, id_pfull /)
axes_zm = (/ id_lat, id_pfull /) !mj

id_ps = register_diag_field(mod_name, 'ps_every', axes_2d, Time, 'surface pressure', 'pascals')
id_u  = register_diag_field(mod_name, 'u_every',  axes_3d, Time, 'zonal wind component', 'm/sec', range=vrange)
id_v  = register_diag_field(mod_name, 'v_every',  axes_3d, Time, 'meridional wind component', 'm/sec', range=vrange)
id_t  = register_diag_field(mod_name, 't_every',  axes_3d, Time, 'temperature', 'deg_k', range=trange)

id_vp  = register_diag_field(mod_name, 'vp',      axes_3d, time, 'meridional wind weighter by ps', 'm/sec', range=vrange)
!id_pt   = register_diag_field(mod_name, 'pottemp',axes_3d, time, 'potential temperature',    'deg_k', range=trange)
id_uv   = register_diag_field(mod_name, 'uv',   axes_3d, time, 'uv',        'm**2/s**2')
id_vq   = register_diag_field(mod_name, 'vq',   axes_3d, time, 'vq',        'm/s')
id_vqint   = register_diag_field(mod_name, 'vqint',   (/id_lon,id_lat/), time, 'vqint', 'm^2/s')
id_vdse   = register_diag_field(mod_name, 'vdse',   axes_3d, time, 'v DSE',        'm/s J/kg')
id_wu  = register_diag_field(mod_name, 'wu',   axes_3d, time, 'omega u',        'Pa/s m/s')
id_wv  = register_diag_field(mod_name, 'wv',   axes_3d, time, 'omega v',        'Pa/s m/s')
id_wq  = register_diag_field(mod_name, 'wq',   axes_3d, time, 'omega q',        'Pa/s')
id_wdse  = register_diag_field(mod_name, 'wdse',   axes_3d, time, 'omega DSE',        'Pa/s J/kg')
id_wup  = register_diag_field(mod_name, 'wup',   axes_3d, time, 'omega u pressure weighted',        'Pa/s m/s')
id_wvp  = register_diag_field(mod_name, 'wvp',   axes_3d, time, 'omega v pressure weighted',        'Pa/s m/s')
id_wqp  = register_diag_field(mod_name, 'wqp',   axes_3d, time, 'omega q pressure weighted',        'Pa/s')
id_wdsep  = register_diag_field(mod_name, 'wdsep',   axes_3d, time, 'omega DSE pressure weighted',        'Pa/s J/kg')
id_wp  = register_diag_field(mod_name, 'omegap',   axes_3d, time, 'dp/dt vertical velocity, pressure weighted',  'Pa/sec')
id_wsubt  = register_diag_field(mod_name, 'wsubt',   axes_3d, time, 'kappa omega T/p, pressure weighted',  'K')
id_wsubtv  = register_diag_field(mod_name, 'wsubtv',   axes_3d, time, 'kappa omega Tv/p, pressure weighted',  'K')
id_kegen  = register_diag_field(mod_name, 'kegen',   axes_3d, time, 'kappa omega T/p, pressure weighted',  'K')
id_kegenq  = register_diag_field(mod_name, 'kegenq',   axes_3d, time, 'kappa omega T q/p, pressure weighted',  'K')
id_kegenqtinv  = register_diag_field(mod_name, 'kegenqtinv',   axes_3d, time, 'kappa omega q/p, pressure weighted',  'K')
id_drystaten =  register_diag_field(mod_name, 'dry_stat_en',  axes_3d, time, 'dry static energy',      'J/kg')
id_moiststaten =  register_diag_field(mod_name, 'moist_stat_en',  axes_3d, time, 'moist static energy',      'J/kg')
id_udt_damp =  register_diag_field(mod_name, 'udt_damp',  axes_3d, time, 'Zonal wind tend from horiz diff',      'm/s**2')
id_vdt_damp =  register_diag_field(mod_name, 'vdt_damp',  axes_3d, time, 'Merid wind tend from horiz diff',      'm/s**2')
id_tdt_damp =  register_diag_field(mod_name, 'tdt_damp',  axes_3d, time, 'Temperature tend from horiz diff',      'm/s**2')
id_entrop_dampuv =  register_diag_field(mod_name, 'entrop_dampuv', &
   axes_3d, time, 'Entropy prod from horiz diff of vel',   '1/s')
id_entrop_dampt =  register_diag_field(mod_name, 'entrop_dampt',   &
   axes_3d, time, 'Entropy prod from horiz diff of temp',      '1/s')
id_tdt_dampuv =  register_diag_field(mod_name, 'tdt_dampuv',  &
   axes_3d, time, 'Temp tend (not applied) due to diff of winds',  'K/s')
id_tdt_tempcor =  register_diag_field(mod_name, 'tdt_tempcor',  &
   axes_3d, time, 'Temp tend due to dynamics temp correction',  'K/s')
id_qdt_watercor =  register_diag_field(mod_name, 'qdt_watercor',  &
   axes_3d, time, 'Humidity tend due to dynamics water correction',  '1/s')
id_entrop_tempcor =  register_diag_field(mod_name, 'entrop_tempcor',  &
   axes_3d, time, 'Entropy tend due to hor diff temp corr',  '1/s')

!mj
id_upvp     = register_diag_field(mod_name, 'upvp'    , axes_3d, time, 'meridional eddy momentum flux', 'm^2/s^2')
id_upwp     = register_diag_field(mod_name, 'upwp'    , axes_3d, time, 'vertical eddy momentum flux'  , 'm/s Pa/s')
id_vpTp     = register_diag_field(mod_name, 'vpTp'    , axes_3d, time, 'meridional eddy heat flux'    , 'm/s K')
id_psi_dwc  = register_diag_field(mod_name, 'psi_dwc' , axes_zm, time, 'residual mean streamfunction' , 'm/s Pa')
id_psi_star = register_diag_field(mod_name, 'psi_star', axes_zm, time, 'residual mean streamfunction' , 'm/s Pa')
!jm

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
allocate(id_tr(num_tracers), id_dt_hadv(num_tracers), id_dt_vadv(num_tracers))
do ntr=1,num_tracers
  call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
  id_tr(ntr) = register_diag_field(mod_name, trim(tname)//'_every', axes_3d, Time, longname, units)
  id_dt_hadv(ntr) = register_diag_field(mod_name,trim(tname)//'_hadv', axes_3d, Time, & ! New name=sphum_hadv XXX
                    'Humidity production due to horizontal advection', trim(units)//' * kg_air/(sec*m^2)')
  id_dt_vadv(ntr) = register_diag_field(mod_name,trim(tname)//'_vadv', axes_3d, Time, & ! New name=sphum_vadv XXX
                    'Humidity production due to vertical advection',   trim(units)//' * kg_air/(sec*m^2)')
enddo
nsphum = get_tracer_index(MODEL_ATMOS,'sphum')

call get_grid_domain(is, ie, js, je)
num_levels = num_levels_in
allocate(two_dt_ps(is:ie, js:je))
allocate(two_dt_u (is:ie, js:je, num_levels))
allocate(two_dt_v (is:ie, js:je, num_levels))
allocate(two_dt_t (is:ie, js:je, num_levels))
allocate(two_dt_tr(is:ie, js:je, num_levels, num_tracers))
allocate(two_dt_id_tr(num_tracers))

two_dt_ps = 0.
two_dt_u  = 0.
two_dt_v  = 0.
two_dt_t  = 0.
two_dt_tr = 0.

two_dt_id_ps = register_static_field(mod_name, '2dt_ps', axes_2d, 'Amplitude of 2*dt wave in surface pressure', 'pascals')
two_dt_id_u  = register_static_field(mod_name, '2dt_u',  axes_3d, 'Amplitude of 2*dt wave in zonal wind', 'm/sec')
two_dt_id_v  = register_static_field(mod_name, '2dt_v',  axes_3d, 'Amplitude of 2*dt wave in meridional wind', 'm/sec')
two_dt_id_t  = register_static_field(mod_name, '2dt_t',  axes_3d, 'Amplitude of 2*dt wave in temperature', 'deg_k')
do ntr=1,num_tracers
  call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
  two_dt_id_tr(ntr) = &
       register_static_field(mod_name, '2dt_'//trim(tname), axes_3d, 'Amplitude of 2*dt wave in '//longname, units)
enddo
iwt = 1
num_time_steps = 0

module_is_initialized = .true.

return
end subroutine every_step_diagnostics_init
!===================================================================================
subroutine every_step_diagnostics(Time, p_surf, u_grid, v_grid, t_grid, tr_grid, &
           wg_full, p_full, p_half, z_full, dt_ug_damp, dt_vg_damp, dt_tg_damp,  &
           temperature_correction, water_correction, dt_hadv, dt_vadv, kegen, kegenq, kegenqtinv)

type(time_type), intent(in) :: Time
real, intent(in), dimension(is:ie, js:je)                         :: p_surf
real, intent(in), dimension(is:ie, js:je, num_levels)             :: u_grid, v_grid, t_grid
real, intent(in), dimension(is:ie, js:je, num_levels,num_tracers) :: tr_grid
real, intent(in), dimension(is:ie, js:je, num_levels)             :: wg_full, p_full, z_full
real, intent(in), dimension(is:ie, js:je, num_levels)             :: dt_ug_damp, dt_vg_damp, dt_tg_damp
real, intent(in), dimension(is:ie, js:je, num_levels+1)           :: p_half
real, intent(in)                                                  :: temperature_correction
real, intent(in), dimension(num_tracers)                          :: dt_hadv, dt_vadv
real, intent(in), dimension(is:ie, js:je, num_levels)             :: water_correction, kegen, kegenq, kegenqtinv

real, dimension(is:ie, js:je, num_levels) :: tempdiag_3d, dse, q_grid
real, dimension(is:ie, js:je)             :: tempdiag_2d

logical :: used
integer :: i, j, k, ntr
character(len=8) :: err_msg_1, err_msg_2

!mj
real, dimension(is:ie, js:je) :: dp
real, dimension(is:ie, js:je, num_levels) :: Thta,vThta
real, dimension(js:je, num_levels) :: tempdiag_zm,dp_u,vpthp,upvp,upwp,fhat
real, dimension(num_levels) :: dphi
real, dimension(js:je) :: dpzm
real :: coslat
!jm

if(.not.module_is_initialized) then
  call error_mesg('every_step_diagnostics','module is not initialized', FATAL)
endif

Time_save = Time

if(id_ps > 0) used = send_data(id_ps, p_surf, Time)
if(id_u  > 0) used = send_data(id_u,  u_grid, Time)
if(id_v  > 0) used = send_data(id_v,  v_grid, Time)
if(id_t  > 0) used = send_data(id_t,  t_grid, Time)

if(id_drystaten > 0 .or. id_vdse > 0 .or. id_wdse > 0 .or. id_wdsep > 0) then
   dse = cp_air*t_grid + grav*z_full ! Dry static energy
endif

q_grid = tr_grid(:,:,:,nsphum)

if(id_drystaten > 0) then
   used = send_data(id_drystaten, dse, Time)
endif
if(id_moiststaten > 0) then
   tempdiag_3d = dse + hlv*q_grid ! Moist static energy
   used = send_data(id_moiststaten, tempdiag_3d, Time)
endif
if(id_uv > 0) then
  do k=1,num_levels
     tempdiag_3d(:,:,k) = u_grid(:,:,k)*v_grid(:,:,k)*p_surf/p00
  enddo
  used = send_data(id_uv, tempdiag_3d, Time)
endif
if(id_vq > 0) then
  do k=1,num_levels
    tempdiag_3d(:,:,k) = v_grid(:,:,k)*q_grid(:,:,k)*p_surf/p00
  enddo
  used = send_data(id_vq, tempdiag_3d, Time)
endif
if(id_vqint > 0) then
  tempdiag_2d(:,:) = 0.
  do k=1,num_levels
     tempdiag_2d(:,:) = tempdiag_2d(:,:) + &
       v_grid(:,:,k)*q_grid(:,:,k)*(p_half(:,:,k+1) - p_half(:,:,k))/grav
  enddo
  used = send_data(id_vqint, tempdiag_2d, Time)
endif
if(id_vdse > 0) then
  do k=1,num_levels
    tempdiag_3d(:,:,k) = v_grid(:,:,k)*dse(:,:,k)*p_surf/p00
  enddo
  used = send_data(id_vdse, tempdiag_3d, Time)
endif
if(id_vp > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = v_grid(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_vp, tempdiag_3d, Time)
endif
if(id_wu > 0) then
   tempdiag_3d = wg_full*u_grid
   used = send_data(id_wu, tempdiag_3d, Time)
endif
if(id_wv > 0) then
   tempdiag_3d = wg_full*v_grid
   used = send_data(id_wv, tempdiag_3d, Time)
endif
if(id_wq > 0) then
   tempdiag_3d = wg_full*q_grid
   used = send_data(id_wq, tempdiag_3d, Time)
endif
if(id_wdse > 0) then
   tempdiag_3d = wg_full*dse
   used = send_data(id_wdse, tempdiag_3d, Time)
endif
if(id_wup > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = wg_full(:,:,k)*u_grid(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_wup, tempdiag_3d, Time)
endif
if(id_wvp > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = wg_full(:,:,k)*v_grid(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_wvp, tempdiag_3d, Time)
endif
if(id_wqp > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = wg_full(:,:,k)*q_grid(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_wqp, tempdiag_3d, Time)
endif
if(id_wdsep > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = wg_full(:,:,k)*dse(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_wdsep, tempdiag_3d, Time)
endif
if(id_wp > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = wg_full(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_wp, tempdiag_3d, Time)
endif
if(id_wsubt > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = kappa*wg_full(:,:,k)*t_grid(:,:,k)/p_full(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_wsubt, tempdiag_3d, Time)
endif
if(id_wsubtv > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = kappa*wg_full(:,:,k)*(t_grid(:,:,k)+virtual_factor*q_grid(:,:,k))/&
         p_full(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_wsubtv, tempdiag_3d, Time)
endif
if(id_kegen > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = kegen(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_kegen, tempdiag_3d, Time)
endif
if(id_kegenq > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = kegenq(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_kegenq, tempdiag_3d, Time)
endif
if(id_kegenqtinv > 0) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = kegenqtinv(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_kegenqtinv, tempdiag_3d, Time)
endif
if(id_udt_damp > 0) then
  used = send_data(id_udt_damp, dt_ug_damp, Time)
endif
if(id_vdt_damp > 0) then
  used = send_data(id_vdt_damp, dt_vg_damp, Time)
endif
if(id_tdt_damp > 0) then
  used = send_data(id_tdt_damp, dt_tg_damp, Time)
endif
if(id_entrop_dampuv > 0) then
  do k=1,num_levels
     tempdiag_3d(:,:,k) = -1./cp_air*(u_grid(:,:,k)*dt_ug_damp(:,:,k) + &
                   v_grid(:,:,k)*dt_vg_damp(:,:,k))*p_surf/p00/t_grid(:,:,k)
  enddo
  used = send_data(id_entrop_dampuv, tempdiag_3d, Time)
endif
if(id_entrop_dampt > 0) then
  do k=1,num_levels
     tempdiag_3d(:,:,k) = dt_tg_damp(:,:,k)/t_grid(:,:,k)*p_surf/p00
  enddo
  used = send_data(id_entrop_dampt, tempdiag_3d, Time)
endif
if(id_tdt_dampuv > 0) then
  tempdiag_3d = -1./cp_air*(u_grid*dt_ug_damp + v_grid*dt_vg_damp)
  used = send_data(id_tdt_dampuv, tempdiag_3d, Time)
endif
if(id_tdt_tempcor > 0 ) then
   tempdiag_3d = temperature_correction*t_grid/t_grid
   used = send_data(id_tdt_tempcor, tempdiag_3d, Time)
endif
if(id_qdt_watercor > 0) then
   used = send_data(id_qdt_watercor, water_correction, Time)
endif
if(id_entrop_tempcor > 0 ) then
   do k=1,num_levels
      tempdiag_3d(:,:,k) = temperature_correction/t_grid(:,:,k)*p_surf/p00
   enddo
   used = send_data(id_entrop_tempcor, tempdiag_3d, Time)
endif
!mj
if(id_psi_dwc > 0) then
  !! F1 = -u'v' + d_puv'Th'/d_pTh
  !! F2 = fhat*v'Th'/d_pTh - u'w'
  ! Thta = Theta_bar
  Thta = t_grid*(p00/p_full)**kappa
  do i=is,ie
    ! vThta = v'Theta' [m/s K]
    vThta(i,:,:) = (v_grid(i,:,:) - sum(v_grid,1)/max(1,size(v_grid,1)))* &
                    (Thta(i,:,:) - sum(Thta,1)/max(1,size(Thta,1)))

  enddo
  ! vThta = v'Theta'/d_pTheta
  ! tempdiag_3d = d_pu
  vThta(:,:,1) = vThta(:,:,1) * &
                (p_full(:,:,1) - p_full(:,:,2))/(Thta(:,:,1) - Thta(:,:,2))
  tempdiag_3d(:,:,1) = ( u_grid(:,:,1)-u_grid(:,:,2) )/( p_full(:,:,1)-p_full(:,:,2))
  do k=2,num_levels-1
    dp = p_full(:,:,k+1) - p_full(:,:,k-1)
    vThta(:,:,k) = vThta(:,:,k) * &
                    dp/(Thta(:,:,k+1) - Thta(:,:,k-1))
    tempdiag_3d(:,:,k) = ( u_grid(:,:,k+1)-u_grid(:,:,k-1) )/dp
  enddo
  dp = p_full(:,:,num_levels) - p_full(:,:,num_levels-1)
  vThta(:,:,num_levels) = vThta(:,:,num_levels) * &
                dp/(Thta(:,:,num_levels) - Thta(:,:,num_levels-1))
  tempdiag_3d(:,:,num_levels) = ( u_grid(:,:,num_levels)-u_grid(:,:,num_levels-1) )/dp
  ! zonal means of v'Theta'/d_pTheta and d_pu
  vpthp = sum(vThta,1)/max(1,size(vThta,1))
  dp_u  = sum(tempdiag_3d,1)/max(1,size(tempdiag_3d,1))
  ! Thta becomes now fhat
  do j=js,je
    coslat = cos(lat(j)/RADIAN)
    if ( j .eq. 1 ) then
      Thta(:,j,:) = (u_grid(:,j+1,:)*cos(lat(j+1)/RADIAN) - u_grid(:,j,:)*coslat) / &
                      ((lat(j+1)-lat(j))/RADIAN)
    elseif ( j .eq. je ) then
      Thta(:,j,:) = (u_grid(:,j,:)*coslat - u_grid(:,j-1,:)*cos(lat(j-1)/RADIAN)) / &
                      ((lat(j)-lat(j-1))/RADIAN)
    else
      Thta(:,j,:) = (u_grid(:,j+1,:)*cos(lat(j+1)/RADIAN) - u_grid(:,j-1,:)*cos(lat(j-1)/RADIAN)) / &
                      ((lat(j+1) - lat(j-1))/RADIAN)
    endif
    Thta(:,j,:) = -Thta(:,j,:)/RADIUS/coslat + 2*OMEGA*sin(lat(j)/RADIAN)
  enddo
  ! zonal mean
  fhat = sum(Thta,1)/max(1,size(Thta,1))
  ! Now put things together for tempdiag_zm = v_star_downward_control
  ! u'v'
  do i=is,ie
    tempdiag_3d(i,:,:) = (u_grid(i,:,:) - sum(u_grid,1)/max(1,size(u_grid,1)))*&
                            (v_grid(i,:,:) - sum(v_grid,1)/max(1,size(v_grid,1)))
  enddo
  upvp = sum(tempdiag_3d,1)/max(1,size(tempdiag_3d,1))
  ! d_xF1
  dphi = (lat(1) - lat(2))/RADIAN
  tempdiag_zm(1,:) = (-(upvp(1,:)-upvp(2,:)) + &
                      (dp_u(1,:)*vpthp(1,:)-dp_u(2,:)*vpthp(2,:))) / dphi
  do j=js+1,je-1
    dphi = (lat(j+1)-lat(j-1))/RADIAN
    tempdiag_zm(j,:) = (-(upvp(j+1,:)-upvp(j-1,:)) + &
                      (dp_u(j+1,:)*vpthp(j+1,:)-dp_u(j-1,:)*vpthp(j-1,:))) / dphi
  enddo
  dphi = (lat(je) - lat(je-1))/RADIAN
  tempdiag_zm(je,:) = (-(upvp(je,:)-upvp(je-1,:)) + &
                      (dp_u(je,:)*vpthp(je,:)-dp_u(je-1,:)*vpthp(je-1,:))) / dphi
  tempdiag_zm = tempdiag_zm / RADIUS
  ! d_pF2
  ! u'w'
  do i=is,ie
    tempdiag_3d(i,:,:) = (u_grid(i,:,:) - sum(u_grid,1)/max(1,size(u_grid,1)))*&
                            (wg_full(i,:,:) - sum(wg_full,1)/max(1,size(wg_full,1)))
  enddo
  upwp = sum(tempdiag_3d,1)/max(1,size(tempdiag_3d,1))
  dpzm = sum(p_full(:,:,1) - p_full(:,:,2),1)/max(1,size(p_full,1))
  tempdiag_zm(:,1) = tempdiag_zm(:,1) + &
                        (fhat(:,1)*vpthp(:,1)-fhat(:,2)*vpthp(:,2) - &
                        (upwp(:,1)-upwp(:,2))) / dpzm
  do k=2,num_levels-1
    dpzm = sum(p_full(:,:,k+1)-p_full(:,:,k-1),1)/max(1,size(p_half,1))
    tempdiag_zm(:,k) = tempdiag_zm(:,k) + &
                          (fhat(:,k+1)*vpthp(:,k+1)-fhat(:,k-1)*vpthp(:,k-1) - &
                          (upwp(:,k+1)-upwp(:,k-1))) / dpzm
  enddo
  dpzm = sum(p_full(:,:,num_levels)-p_full(:,:,num_levels-1),1)/max(1,size(p_full,1))
  tempdiag_zm(:,num_levels) = tempdiag_zm(:,num_levels) + &
                                (fhat(:,num_levels)*vpthp(:,num_levels)-fhat(:,num_levels-1)*vpthp(:,num_levels-1) - &
                                (upwp(:,num_levels)-upwp(:,num_levels-1))) / dpzm
  ! v_star_downward_control = 1/fhat*(-diff(F))
  tempdiag_zm = -tempdiag_zm / fhat
  !now compute psi_star_dwc using the trapezoidal rule
  !temporarily, upvp = v_star_downward_control
  upvp = tempdiag_zm
  tempdiag_zm = 0.
  do k=2,num_levels
    dpzm = sum(p_full(:,:,k) - p_full(:,:,k-1),1)/max(1,size(p_full,1))
    tempdiag_zm(:,k) = tempdiag_zm(:,k-1) + 0.5*dpzm*(upvp(:,k)+upvp(:,k-1))
  enddo
  used = send_data(id_psi_dwc, tempdiag_zm, Time)
endif
if(id_upvp > 0) then
  do k=1,num_levels
    do i=is,ie
      tempdiag_3d(i,:,k) = (u_grid(i,:,k) - sum(u_grid(:,:,k),1)/max(1,size(u_grid,1)))* &
                           (v_grid(i,:,k) - sum(v_grid(:,:,k),1)/max(1,size(v_grid,1)))
    enddo
  enddo
  used = send_data(id_upvp, tempdiag_3d, Time)
endif
if(id_upwp > 0) then
  do k=1,num_levels
    do i=is,ie
      tempdiag_3d(i,:,k) = (u_grid(i,:,k) - sum(u_grid(:,:,k),1)/max(1,size(u_grid,1)))* &
                           (wg_full(i,:,k)- sum(wg_full(:,:,k),1)/max(1,size(wg_full,1)))
    enddo
  enddo
  used = send_data(id_upwp, tempdiag_3d, Time)
endif
if(id_vpTp > 0) then
  do k=1,num_levels
    do i=is,ie
      tempdiag_3d(i,:,k) = (v_grid(i,:,k) - sum(v_grid(:,:,k),1)/max(1,size(v_grid,1)))* &
                           (t_grid(i,:,k) - sum(t_grid(:,:,k),1)/max(1,size(t_grid,1)))
    enddo
  enddo
  used = send_data(id_vpTp, tempdiag_3d, Time)
endif
if(id_psi_star > 0) then
  do i=is,ie
    !tempdiag_3d = v'Theta' [m/s K]
    Thta(i,:,:) = (t_grid(i,:,:) - sum(t_grid,1)/max(1,size(t_grid,1)))* &
                  (p00/p_full(i,:,:))**kappa
    tempdiag_3d(i,:,:) = (v_grid(i,:,:) - sum(v_grid,1)/max(1,size(v_grid,1)))* &
                  Thta(i,:,:)
  enddo
  !now, tempdiag_3d = v'Theta'/dpTheta [m/s Pa]
  tempdiag_3d(:,:,1) = (Thta(:,:,2) - Thta(:,:,1))/(p_full(:,:,2)-p_full(:,:,1))
  do k=2,num_levels-1
    dp = p_full(:,:,k+1) - p_full(:,:,k-1)
    tempdiag_3d(:,:,k) = tempdiag_3d(:,:,k)*dp/(Thta(:,:,k+1)-Thta(:,:,k-1))
  enddo
  tempdiag_3d(:,:,num_levels) = (Thta(:,:,num_levels) - Thta(:,:,num_levels-1))/&
                                (p_full(:,:,num_levels) - p_full(:,:,num_levels-1))
  !store this in the 2D variable
  tempdiag_zm = sum(tempdiag_3d,3)/max(1,size(tempdiag_3d,3))
  !now compute psi using the trapezoidal rule
  tempdiag_3d(:,:,1) = 0.
  do k=2,num_levels
    dp = p_full(:,:,k) - p_full(:,:,k-1)
    tempdiag_3d(:,:,k) = tempdiag_3d(:,:,k-1) + 0.5*dp*(v_grid(:,:,k)+v_grid(:,:,k-1))
  enddo
  !now put it all together
  tempdiag_zm = sum(tempdiag_3d,1)/max(1,size(tempdiag_3d,1)) - tempdiag_zm
  used = send_data(id_psi_star, tempdiag_zm, Time)
endif
!jm

!--------------------
!if(id_pt > 0) then
!   tempdiag_3d = t_grid*(p00/p_full)**kappa
!   used = send_data(id_pt, tempdiag_3d, Time)
!endif
!--------------------

if(size(tr_grid,4) /= num_tracers) then
  write(err_msg_1,'(i8)') size(tr_grid,4)
  write(err_msg_2,'(i8)') num_tracers
  call error_mesg('every_step_diagnostics','size(tracers)='//err_msg_1//' Should be='//err_msg_2, FATAL)
endif
do ntr=1,num_tracers
  if(id_tr(ntr) > 0) used = send_data(id_tr(ntr), tr_grid(:,:,:,ntr), Time)
  if(id_dt_hadv(ntr) > 0) used = send_data(id_dt_hadv(ntr), (dt_hadv(ntr)*u_grid/u_grid), Time) ! Why output a 3d array? XXX
  if(id_dt_vadv(ntr) > 0) used = send_data(id_dt_vadv(ntr), (dt_vadv(ntr)*u_grid/u_grid), Time) ! Why output a 3d array? XXX
enddo

if(two_dt_id_ps > 0) then
  do j=js,je
    do i=is,ie
      two_dt_ps(i,j) = two_dt_ps(i,j) + iwt*p_surf(i,j)
    enddo
  enddo
endif
if(two_dt_id_u  > 0) then
  do k=1,num_levels
    do j=js,je
      do i=is,ie
        two_dt_u(i,j,k) = two_dt_u(i,j,k) + iwt*u_grid(i,j,k)
      enddo
    enddo
  enddo
endif
if(two_dt_id_v  > 0) then
  do k=1,num_levels
    do j=js,je
      do i=is,ie
        two_dt_v(i,j,k) = two_dt_v(i,j,k) + iwt*v_grid(i,j,k)
      enddo
    enddo
  enddo
endif
if(two_dt_id_t  > 0) then
  do k=1,num_levels
    do j=js,je
      do i=is,ie
        two_dt_t(i,j,k) = two_dt_t(i,j,k) + iwt*t_grid(i,j,k)
      enddo
    enddo
  enddo
endif
do ntr=1,num_tracers
  if(two_dt_id_tr(ntr) > 0) then
    do k=1,num_levels
      do j=js,je
        do i=is,ie
          two_dt_tr(i,j,k,ntr) = two_dt_tr(i,j,k,ntr) + iwt*tr_grid(i,j,k,ntr)
        enddo
      enddo
    enddo
  endif
enddo

iwt = -1*iwt
num_time_steps = num_time_steps + 1

return
end subroutine every_step_diagnostics
!===================================================================================
subroutine every_step_diagnostics_end(Time_in)
type(time_type), intent(in), optional :: Time_in
logical :: used
integer :: ntr
type(time_type) :: Time

if(.not.module_is_initialized) return

if(present(Time_in)) then
  Time = Time_in
else
  Time = Time_save
endif

two_dt_ps = two_dt_ps / num_time_steps
two_dt_u  = two_dt_u  / num_time_steps
two_dt_v  = two_dt_v  / num_time_steps
two_dt_t  = two_dt_t  / num_time_steps
two_dt_tr = two_dt_tr / num_time_steps

if(two_dt_id_ps > 0) used = send_data(two_dt_id_ps, two_dt_ps, Time)
if(two_dt_id_u  > 0) used = send_data(two_dt_id_u,  two_dt_u,  Time)
if(two_dt_id_v  > 0) used = send_data(two_dt_id_v,  two_dt_v,  Time)
if(two_dt_id_t  > 0) used = send_data(two_dt_id_t,  two_dt_t,  Time)
do ntr=1,num_tracers
  if(two_dt_id_tr(ntr) > 0) used = send_data(two_dt_id_tr(ntr), two_dt_tr(:,:,:,ntr), Time)
enddo

deallocate(id_tr, two_dt_id_tr)
deallocate(two_dt_ps, two_dt_u, two_dt_v, two_dt_t, two_dt_tr)
module_is_initialized = .false.

return
end subroutine every_step_diagnostics_end
!===================================================================================

end module every_step_diagnostics_mod
