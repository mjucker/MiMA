
module vert_diff_driver_mod

!-----------------------------------------------------------------------
!   module performs vertical diffusion of atmospheric variables
!-----------------------------------------------------------------------

use    vert_diff_mod, only:  surf_diff_type,     &
                             vert_diff_init, &
                             gcm_vert_diff_down, &
                             gcm_vert_diff_up,   &
                             vert_diff_end


use diag_manager_mod, only:  register_diag_field, send_data

use time_manager_mod, only:  time_type

use          fms_mod, only:  file_exist, open_namelist_file, error_mesg,  &
                             check_nml_error, FATAL, mpp_pe, mpp_root_pe, &
                             close_file, write_version_number, stdlog

use    constants_mod, only:  CP_AIR, GRAV

!-----------------------------------------------------------------------

implicit none
private

public :: vert_diff_driver_down, vert_diff_driver_up,  &
          vert_diff_driver_init, vert_diff_driver_end
public :: surf_diff_type


!-----------------------------------------------------------------------
!---- namelist ----

logical :: do_conserve_energy         = .false.
logical :: do_mcm_no_neg_q            = .false.
logical :: use_virtual_temp_vert_diff = .true.
logical :: do_mcm_plev                = .false.
logical :: do_mcm_vert_diff_tq        = .false.

namelist /vert_diff_driver_nml/ do_conserve_energy,         &
                                do_mcm_no_neg_q,            &
                                use_virtual_temp_vert_diff, &
                                do_mcm_plev, do_mcm_vert_diff_tq

!-----------------------------------------------------------------------

real, allocatable, dimension(:,:,:) :: dt_t_save, dt_q_save

!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_vdif, id_qdt_vdif, id_udt_vdif, id_vdt_vdif,  &
           id_sens_vdif, id_evap_vdif,                          &
           id_tdt_diss_vdif, id_diss_heat_vdif, id_entrop_vdif_sens, id_entrop_vdif_kediss

real :: missing_value = -999.

character(len=9), parameter :: mod_name = 'vert_diff'

!-----------------------------------------------------------------------
!---- version number ----

character(len=128) :: version = '$Id: vert_diff_driver.f90,v 12.0.4.1 2005/05/13 18:16:38 pjp Exp $'
character(len=128) :: tagname = '$Name:  $'

logical :: module_is_initialized = .false.


contains

!#######################################################################

 subroutine vert_diff_driver_down (is, js, Time, delt, p_half, p_full, &
                                   z_full, diff_mom, diff_heat,        &
                                   u, v, t, q, trs,                    &
                                   dtau_du, dtau_dv, tau_x, tau_y,     &
                                   dt_u, dt_v, dt_t, dt_q, dt_trs,     &
                                   Surf_diff,  mask, kbot              )

integer, intent(in)                     :: is, js
type(time_type),   intent(in)           :: Time
real, intent(in)                        :: delt
real, intent(in)   , dimension(:,:,:)   :: p_half, p_full, z_full,  &
                                           diff_mom, diff_heat
real, intent(in),    dimension(:,:,:)   :: u, v, t, q
real, intent(in),    dimension(:,:,:,:) :: trs
real, intent(in),    dimension(:,:)     :: dtau_du, dtau_dv

real, intent(inout), dimension(:,:)     :: tau_x, tau_y
real, intent(inout), dimension(:,:,:)   :: dt_u, dt_v, dt_t, dt_q
real, intent(inout), dimension(:,:,:,:) :: dt_trs

type(surf_diff_type), intent(inout)     :: Surf_diff

real   , intent(in), dimension(:,:,:), optional :: mask
integer, intent(in), dimension(:,:),   optional :: kbot

real, dimension(size(t,1),size(t,2),size(t,3)) :: tt, dpg, q_2, entrop_vdif_sens, entrop_vdif_kediss
real, dimension(size(t,1),size(t,2),size(t,3)) :: dissipative_heat
integer :: k, ntp, kx
logical :: used
real, dimension(size(t,1),size(t,2)) :: diag2
integer :: ie, je

!-----------------------------------------------------------------------

  if (.not. module_is_initialized) call error_mesg       &
                  ('vert_diff_driver_mod',  &
                   'vert_diff_driver_init must be called first', FATAL)

!-----------------------------------------------------------------------

  ntp = size(dt_trs,4) ! number of prognostic tracers
  if (size(trs,4) < ntp) call error_mesg ('vert_diff_driver', &
             'Number of tracers .lt. number of tracer tendencies',FATAL)

  ie = is + size(t,1) -1
  je = js + size(t,2) -1


    if(do_mcm_vert_diff_tq) then
      dt_t_save(is:ie,js:je,:) = dt_t
      dt_q_save(is:ie,js:je,:) = dt_q
      dt_t = 0.0
      dt_q = 0.0
    endif

!-----------------------------------------------------------------------
!---- to do diagnostics on dt_t, dt_q, dt_u, and dt_v at this point add 
!-----in the negative value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------




!------- diagnostics for uwnd_diff -------
    if ( id_udt_vdif > 0 ) then
       used = send_data ( id_udt_vdif, -2.*dt_u, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for vwnd_diff -------
    if ( id_vdt_vdif > 0 ) then
       used = send_data ( id_vdt_vdif, -2.*dt_v, Time, is, js, 1, &
                          rmask=mask )
    endif



!-----------------------------------------------------------------------
!---- local temperature ----
!     (avoids divid by zero when computing virtual temp)

   tt = t
   if (present(mask)) where (mask < 0.5) tt = 200.

!-----------------------------------------------------------------------
!---- momentum diffusion ----
!---- heat/moisture diffusion (down only) ----
!---- tracer diffusion (no surface flux) ----

 q_2 = q
 if (do_mcm_no_neg_q) then
   where (q_2 < 0.0)  q_2 = 0.0
 endif

 call gcm_vert_diff_down (is, js, delt, u, v, tt, q_2, trs(:,:,:,1:ntp), &
                          diff_mom, diff_heat, p_half, p_full, z_full,   &
                          tau_x, tau_y, dtau_du, dtau_dv,                &
                          dt_u, dt_v, dt_t, dt_q, dt_trs(:,:,:,1:ntp), &
                          dissipative_heat, Surf_diff,  kbot           )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---- to do diagnostics on dt_u, and dt_v at this point add 
!-----in the value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------

!------ preliminary calculations for vert integrals -----
    if ( id_sens_vdif > 0 .or. id_evap_vdif > 0 .or. id_diss_heat_vdif > 0 ) then
            do k = 1, size(p_half,3)-1
               dpg(:,:,k) = (p_half(:,:,k+1)-p_half(:,:,k))/GRAV
            enddo
            if (present(mask)) dpg = dpg*mask
    endif
    
!------- diagnostics for dt/dt_diff -------
    if ( id_tdt_vdif > 0 ) then
       used = send_data ( id_tdt_vdif, -2.*dt_t, Time, is, js, 1, &
                           rmask=mask )
    endif
    if ( id_entrop_vdif_sens > 0 ) then
       kx = size(p_half,3)-1
       do k=1,kx
          entrop_vdif_sens(:,:,k) = -2.*dt_t(:,:,k)/t(:,:,k)*p_half(:,:,kx)/1.e5
       enddo
       used = send_data ( id_entrop_vdif_sens, entrop_vdif_sens, Time, is, js, 1, &
                           rmask=mask )
    endif
    
!------- diagnostics for dq/dt_diff -------
    if ( id_qdt_vdif > 0 ) then
       used = send_data ( id_qdt_vdif, -2.*dt_q, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for sens_diff -------
    if ( id_sens_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_t*dpg, 3 )
          used = send_data ( id_sens_vdif, -2.*CP_AIR*diag2, Time, is, js )
    endif

!------- diagnostics for evap_diff -------
    if ( id_evap_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_q*dpg, 3 )
          used = send_data ( id_evap_vdif, -2.*diag2, Time, is, js )
    endif

!------- diagnostics for uwnd_diff -------
    if ( id_udt_vdif > 0 ) then
       used = send_data ( id_udt_vdif, 2.*dt_u, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for vwnd_diff -------
    if ( id_vdt_vdif > 0 ) then
       used = send_data ( id_vdt_vdif, 2.*dt_v, Time, is, js, 1, &
                          rmask=mask )
    endif
    
!------- diagnostics for dissipative heating -------
    if ( id_tdt_diss_vdif > 0 ) then
       used = send_data ( id_tdt_diss_vdif, dissipative_heat, Time, &
                          is, js, 1, &
                          rmask=mask)
    endif
    if ( id_entrop_vdif_kediss > 0 ) then
       kx = size(p_half,3)-1
       do k=1,kx
          entrop_vdif_kediss(:,:,k) = dissipative_heat(:,:,k)/t(:,:,k)*p_half(:,:,kx)/1.e5
       enddo

       used = send_data ( id_entrop_vdif_kediss, entrop_vdif_kediss, Time, &
                          is, js, 1, &
                          rmask=mask)
    endif

!------- diagnostics for vertically integrated dissipative heating -------
    if ( id_diss_heat_vdif > 0 ) then
          diag2 = sum( CP_AIR*dissipative_heat*dpg, 3 )
          used = send_data ( id_diss_heat_vdif, diag2, Time, is, js )
    endif

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_down

!#######################################################################

 subroutine vert_diff_driver_up (is, js, Time, delt, p_half, &
                                 Surf_diff, dt_t, dt_q, mask, kbot, t)

 integer,           intent(in)            :: is, js
 type(time_type),   intent(in)            :: Time
 real,    intent(in)                      :: delt
 real,    intent(in),    dimension(:,:,:) :: p_half
 type(surf_diff_type),   intent(in)       :: Surf_diff
 real,    intent(inout), dimension(:,:,:) :: dt_t, dt_q
 real   , intent(in), dimension(:,:,:), optional :: mask
 integer, intent(in),   dimension(:,:), optional :: kbot
 real,    intent(in), dimension(:,:,:), optional :: t

 integer :: k, kx
 logical :: used
 real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: dpg, entrop_vdif_sens
 real, dimension(size(p_half,1),size(p_half,2)) :: diag2
 integer :: ie, je

!-----------------------------------------------------------------------
    ie = is + size(p_half,1) -1
    je = js + size(p_half,2) -1
!-----------------------------------------------------------------------

    call gcm_vert_diff_up (is, js, delt, Surf_diff, dt_t, dt_q, kbot)

!-----------------------------------------------------------------------
!---- to do diagnostics on dt_t and dt_q at this point add in the
!---- the postive value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------
!------- diagnostics for dt/dt_diff -------

    if ( id_tdt_vdif > 0 ) then
       used = send_data ( id_tdt_vdif, 2.*dt_t, Time, is, js, 1, &
                          rmask=mask )
    endif
    if ( present(t)) then
      if ( id_entrop_vdif_sens > 0 ) then
         kx = size(p_half,3)-1
         do k=1,kx
            entrop_vdif_sens(:,:,k) = 2.*dt_t(:,:,k)/t(:,:,k)*p_half(:,:,kx)/1.e5
         enddo
         used = send_data ( id_entrop_vdif_sens, entrop_vdif_sens, Time, is, js, 1, &
                             rmask=mask )
      endif
    endif

!------- diagnostics for dq/dt_diff -------
    if ( id_qdt_vdif > 0 ) then
       used = send_data ( id_qdt_vdif, 2.*dt_q, Time, is, js, 1, &
                          rmask=mask )
    endif

!------ preliminary calculations for vert integrals -----
    if ( id_sens_vdif > 0 .or. id_evap_vdif > 0 ) then
            do k = 1, size(p_half,3)-1
               dpg(:,:,k) = (p_half(:,:,k+1)-p_half(:,:,k))/GRAV
            enddo
            if (present(mask)) dpg = dpg*mask
    endif

!------- diagnostics for sens_diff -------
    if ( id_sens_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_t*dpg, 3 )
          used = send_data ( id_sens_vdif, 2.*CP_AIR*diag2, Time, is, js )
    endif

!------- diagnostics for evap_diff -------
    if ( id_evap_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_q*dpg, 3 )
          used = send_data ( id_evap_vdif, 2.*diag2, Time, is, js )
    endif

    if(do_mcm_vert_diff_tq) then
      dt_t = dt_t + dt_t_save(is:ie,js:je,:)
      dt_q = dt_q + dt_q_save(is:ie,js:je,:)
    endif

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_up

!#######################################################################

 subroutine vert_diff_driver_init ( Surf_diff, idim, jdim, kdim,  &
                                    axes, Time )

 type(surf_diff_type), intent(inout) :: Surf_diff
 integer             , intent(in)    :: idim, jdim, kdim, axes(4)
 type     (time_type), intent(in)    :: Time

 integer :: unit, io, ierr

!-----------------------------------------------------------------------
!------ read namelist ------

   if ( file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=vert_diff_driver_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'vert_diff_driver_nml')
      enddo
 10   call close_file (unit)
   endif

!--------- write version number and namelist ------------------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) write(stdlog(),nml=vert_diff_driver_nml)

!-------- initialize gcm vertical diffusion ------

   call vert_diff_init (Surf_diff, idim, jdim, kdim, do_conserve_energy, &
                        use_virtual_temp_vert_diff, do_mcm_plev)

!-----------------------------------------------------------------------

   if(do_mcm_vert_diff_tq) then
     allocate(dt_t_save(idim,jdim,kdim)) ; dt_t_save = 0.0
     allocate(dt_q_save(idim,jdim,kdim)) ; dt_q_save = 0.0
   endif

!--------------- initialize diagnostic fields --------------------

   id_entrop_vdif_sens = &
   register_diag_field ( mod_name, 'entrop_vdif_sens', axes(1:3), Time, &
                         'Entropy tendency from vert diff sensible heating',  &
                         '1/s', missing_value=missing_value  )

   id_entrop_vdif_kediss = &
   register_diag_field ( mod_name, 'entrop_vdif_kediss', axes(1:3), Time, &
                         'Entropy tendency from vert diff kinetic energy dissipation',  &
                         '1/s', missing_value=missing_value  )

   id_tdt_vdif = &
   register_diag_field ( mod_name, 'tdt_vdif', axes(1:3), Time, &
                        'Temperature tendency from vert diff',  &
                        'deg_K/s', missing_value=missing_value  )

   id_qdt_vdif = &
   register_diag_field ( mod_name, 'qdt_vdif', axes(1:3), Time, &
                        'Spec humidity tendency from vert diff',&
                        'kg/kg/s', missing_value=missing_value  )

   id_udt_vdif = &
   register_diag_field ( mod_name, 'udt_vdif', axes(1:3), Time, &
                        'Zonal wind tendency from vert diff',   &
                        'm/s2', missing_value=missing_value     )

   id_vdt_vdif = &
   register_diag_field ( mod_name, 'vdt_vdif', axes(1:3), Time,    &
                        'Meridional wind tendency from vert diff', &
                        'm/s2', missing_value=missing_value        )

   id_sens_vdif = &
   register_diag_field ( mod_name, 'sens_vdif', axes(1:2), Time,  &
                        'Integrated heat flux from vert diff',    &
                        'W/m2' )

   id_evap_vdif = &
   register_diag_field ( mod_name, 'evap_vdif', axes(1:2), Time,    &
                        'Integrated moisture flux from vert diff',  &
                        'kg/m2/s' )

   id_tdt_diss_vdif = &
   register_diag_field ( mod_name, 'tdt_diss_vdif', axes(1:3), Time,  &
                        'Dissipative heating from vert_diff', 'deg_K/s', &
                         missing_value=missing_value  ) 

   id_diss_heat_vdif = &
   register_diag_field ( mod_name, 'diss_heat_vdif', axes(1:2), Time,  &
                        'Integrated dissipative heating from vert diff',  &
                        'W/m2' )


!-----------------------------------------------------------------------

   module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_init

!#######################################################################

 subroutine vert_diff_driver_end

   call vert_diff_end
   if(do_mcm_vert_diff_tq) deallocate(dt_t_save, dt_q_save)

!-----------------------------------------------------------------------

   module_is_initialized = .false.

!-----------------------------------------------------------------------
 end subroutine vert_diff_driver_end

!#######################################################################

end module vert_diff_driver_mod

