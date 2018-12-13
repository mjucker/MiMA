module spectral_initialize_fields_mod

! epg: we added "error_mesg, FATAL, and file_exist" here so that we can report an error if 
!      the initial_conditions.nc file is missing.
use              fms_mod, only: mpp_pe, mpp_root_pe, write_version_number, file_exist, FATAL, error_mesg

use        constants_mod, only: rdgas

use       transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, vor_div_from_uv_grid, &
                                uv_grid_from_vor_div, get_grid_domain, get_spec_domain, area_weighted_global_mean

implicit none
private

public :: spectral_initialize_fields

character(len=128), parameter :: version = &
'$Id: spectral_initialize_fields.f90,v 10.0 2003/10/24 22:00:59 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: lima $'

logical :: entry_to_logfile_done = .false.

! epg: this netcdf include file is needed to load in specified initial
! conditions.  Only used if choice_of_init == 3
include 'netcdf.inc'

contains

!-------------------------------------------------------------------------------------------------
subroutine spectral_initialize_fields(reference_sea_level_press, triang_trunc, choice_of_init, initial_temperature, &
                        surf_geopotential, ln_ps, vors, divs, ts, psg, ug, vg, tg, vorg, divg)

real,    intent(in) :: reference_sea_level_press
logical, intent(in) :: triang_trunc
integer, intent(in) :: choice_of_init
real,    intent(in) :: initial_temperature

real,    intent(in),  dimension(:,:    ) :: surf_geopotential
complex, intent(out), dimension(:,:    ) :: ln_ps
complex, intent(out), dimension(:,:,:  ) :: vors, divs, ts
real,    intent(out), dimension(:,:    ) :: psg
real,    intent(out), dimension(:,:,:  ) :: ug, vg, tg
real,    intent(out), dimension(:,:,:  ) :: vorg, divg

real, allocatable, dimension(:,:) :: ln_psg

real :: initial_sea_level_press, global_mean_psg
real :: initial_perturbation   = 1.e-7

integer :: ms, me, ns, ne, is, ie, js, je, num_levels

! epg: needed to load in initial conditions from a netcdf file
!      code was initially developed by Lorenzo Polvani; hence lmp
real, allocatable,dimension(:,:,:) :: lmptmp
integer :: ncid,vid,err,counts(3)
! --------

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

num_levels = size(ug,3)
call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)
allocate(ln_psg(is:ie, js:je))

initial_sea_level_press = reference_sea_level_press  

ug      = 0.
vg      = 0.
tg      = 0.
psg     = 0.
vorg    = 0.
divg    = 0.

vors  = (0.,0.)
divs  = (0.,0.)
ts    = (0.,0.)
ln_ps = (0.,0.)

tg     = initial_temperature
ln_psg = log(initial_sea_level_press) - surf_geopotential/(rdgas*initial_temperature)
psg    = exp(ln_psg)

if(choice_of_init == 1) then  ! perturb temperature field
  if(is <= 1 .and. ie >= 1 .and. js <= 1 .and. je >= 1) then
    tg(1,1,:) = tg(1,1,:) + 1.0
  endif
endif

if(choice_of_init == 2) then   ! initial vorticity perturbation used in benchmark code
  if(ms <= 1 .and. me >= 1 .and. ns <= 3 .and. ne >= 3) then
    vors(2-ms,4-ns,num_levels  ) = initial_perturbation
    vors(2-ms,4-ns,num_levels-1) = initial_perturbation
    vors(2-ms,4-ns,num_levels-2) = initial_perturbation
  endif
  if(ms <= 5 .and. me >= 5 .and. ns <= 3 .and. ne >= 3) then
    vors(6-ms,4-ns,num_levels  ) = initial_perturbation
    vors(6-ms,4-ns,num_levels-1) = initial_perturbation
    vors(6-ms,4-ns,num_levels-2) = initial_perturbation
  endif
  if(ms <= 1 .and. me >= 1 .and. ns <= 2 .and. ne >= 2) then
    vors(2-ms,3-ns,num_levels  ) = initial_perturbation
    vors(2-ms,3-ns,num_levels-1) = initial_perturbation
    vors(2-ms,3-ns,num_levels-2) = initial_perturbation
  endif
  if(ms <= 5 .and. me >= 5 .and. ns <= 2 .and. ne >= 2) then
    vors(6-ms,3-ns,num_levels  ) = initial_perturbation
    vors(6-ms,3-ns,num_levels-1) = initial_perturbation
    vors(6-ms,3-ns,num_levels-2) = initial_perturbation
  endif
  call uv_grid_from_vor_div(vors, divs, ug, vg)
endif

! epg: This was written by Lorenzo Polvani to load in the initial conditions
!      from a netcdf file, which must be called "initial_conditions.nc" and must be placed in the
!      INPUT directory from where the code is being run.
If (choice_of_init == 3) then !initialize with prescribed input
   if (.not.file_exist('INPUT/initial_conditions.nc')) then
      call error_mesg('spectral_initialize_fields','Could not find INPUT/initial_conditions.nc!',FATAL)
   end if

   ! first, open up the netcdf file
   ncid = ncopn('INPUT/initial_conditions.nc',NCNOWRIT,err)
   ! This array tells us the size of input variables.
   counts(1) = size(ug,1)
   counts(2) = size(ug,2)
   counts(3) = size(ug,3)
   ! Allocate space to put the initial condition information, temporarily.
   allocate(lmptmp(counts(1),counts(2),counts(3)))
   
   ! load the zonal wind initial conditions
   vid = ncvid(ncid,'ucomp',err)
   call ncvgt(ncid,vid,(/is,js,1/),counts,lmptmp,err)
   ug(:,:,:) = lmptmp
 
   ! load the meridional wind
   vid = ncvid(ncid,'vcomp',err)
   call ncvgt(ncid,vid,(/is,js,1/),counts,lmptmp,err)
   vg(:,:,:) = lmptmp
 
   ! load temperature
   vid = ncvid(ncid,'temp',err)
   call ncvgt(ncid,vid,(/is,js,1/),counts,lmptmp,err)
   tg(:,:,:) = lmptmp

   ! load surface pressure
   vid = ncvid(ncid,'ps',err)
   call ncvgt(ncid,vid,(/is,js/),counts(1:2),lmptmp(:,:,1),err)
   psg(:,:) = lmptmp(:,:,1)
   ln_psg = log(psg(:,:))
 
   ! close up the input netcdf file.
   call ncclos(ncid,err)
 
   ! and lastly, let us know that it worked!
   if(mpp_pe() == mpp_root_pe()) then
      print *, 'initial dynamical fields read in from initial_conditions.nc'
   endif
 
endif
!epg: end of lorenzo polvani's script --------


!  initial spectral fields (and spectrally-filtered) grid fields

call trans_grid_to_spherical(tg, ts)
call trans_spherical_to_grid(ts, tg)

call trans_grid_to_spherical(ln_psg, ln_ps)
call trans_spherical_to_grid(ln_ps,  ln_psg)
psg = exp(ln_psg)

call vor_div_from_uv_grid(ug, vg, vors, divs, triang=triang_trunc)
call uv_grid_from_vor_div(vors, divs, ug, vg)
call trans_spherical_to_grid(vors, vorg)
call trans_spherical_to_grid(divs, divg)

!  compute and print mean surface pressure
global_mean_psg = area_weighted_global_mean(psg)
if(mpp_pe() == mpp_root_pe()) then
  print '("mean surface pressure=",f9.4," mb")',.01*global_mean_psg
endif

return
end subroutine spectral_initialize_fields
!================================================================================

end module spectral_initialize_fields_mod
