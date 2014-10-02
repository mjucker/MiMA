!FDOC_TAG_GFDL

module cu_mo_trans_mod
! <CONTACT EMAIL="Isaac.Held@noaa.gov">
!  Isaac Held
! </CONTACT>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    A simple module that computes a diffusivity proportional to the 
!    convective mass flux, for use with diffusive 
!    convective momentum transport closure
! </OVERVIEW>
! <DESCRIPTION>
!   A diffusive approximation to convective momentum transport is crude but
!    has been found to be useful in improving the simulation of tropical
!     precipitation in some models.  The diffusivity computed here is
!     simply 
!<PRE>
! diffusivity = c*W*L 
! W = M/rho  (m/sec) 
! M = convective mass flux (kg/(m2 sec)) 
! rho - density of air <p>
! L = depth of convecting layer (m)
! c = normalization constant = diff_norm/g 
!   (diff_norm is a namelist parameter;
!      the factor of g = acceleration of gravity here is an historical artifact) <p>
! for further discussion see 
!     <LINK SRC="cu_mo_trans.pdf">cu_mo_trans.pdf</LINK>
!</PRE>
! </DESCRIPTION>


!=======================================================================
!
!                 DIFFUSIVE CONVECTIVE MOMENTUM TRANSPORT MODULE
!
!=======================================================================

  use   constants_mod, only:  GRAV, RDGAS, RVGAS
 

  use         fms_mod, only: file_exist, check_nml_error,    &
                             open_namelist_file, close_file, &
                             write_version_number,           &
                             mpp_pe, mpp_root_pe, stdlog,    &
                             error_mesg, FATAL, NOTE

  use  Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use  Time_Manager_Mod, ONLY: time_type

implicit none
private


! public interfaces
!=======================================================================
public :: cu_mo_trans_init, &
          cu_mo_trans,      &
          cu_mo_trans_end

!=======================================================================

! form of interfaces
!=======================================================================

      
logical :: module_is_initialized = .false.


!---------------diagnostics fields------------------------------------- 

integer :: id_diff_cmt

character(len=11) :: mod_name = 'cu_mo_trans'

real :: missing_value = -999.

!--------------------- namelist variables with defaults -------------

real    :: diff_norm =   1.0

namelist/cu_mo_trans_nml/ diff_norm


!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: cu_mo_trans.f90,v 10.0 2003/10/24 22:00:25 fms Exp $'
character(len=128) :: tagname = '$Name: lima $'

contains

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans_init">
!  <OVERVIEW>
!   initializes module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Reads namelist and registers one diagnostic field
!     (diff_cmt:  the kinematic diffusion coefficient)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans_init( axes, Time )
!		
!  </TEMPLATE>
!  <IN NAME=" axes" TYPE="integer">
!    axes identifier needed by diag manager
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!    time at initialization needed by diag manager
!  </IN>
! </SUBROUTINE>
!
subroutine cu_mo_trans_init( axes, Time )

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time

integer :: unit, ierr, io

!------ read namelist ------

   if ( file_exist('input.nml')) then
         unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=cu_mo_trans_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'cu_mo_trans_nml')
      enddo
 10   call close_file (unit)
   endif

!--------- write version number and namelist ------------------

      call write_version_number ( version, tagname )
      if ( mpp_pe() == mpp_root_pe() ) &
      write ( stdlog(), nml=cu_mo_trans_nml )

! --- initialize quantities for diagnostics output -------------

   id_diff_cmt = &
   register_diag_field ( mod_name, 'diff_cmt', axes(1:3), Time,    &
                        'cu_mo_trans coeff for momentum',  'm2/s', &
                         missing_value=missing_value               )

!--------------------------------------------------------------

  module_is_initialized = .true.


end subroutine cu_mo_trans_init

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans_end">
!  <OVERVIEW>
!   terminates module
!  </OVERVIEW>
!  <DESCRIPTION>
!   This is the destructor for cu_mo_trans
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cu_mo_trans_end()

  module_is_initialized = .false.

end subroutine cu_mo_trans_end

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans">
!  <OVERVIEW>
!   returns a diffusivity proportional to the 
!    convective mass flux, for use with diffusive 
!    convective momentum transport closure
!  </OVERVIEW>
!  <DESCRIPTION>
!   A diffusive approximation to convective momentum transport is crude but
!    has been found to be useful in inproving the simulation of tropical
!    precipitation in some models.  The diffusivity computed here is
!    simply 
!<PRE>
! diffusivity = c*W*L
! W = M/rho  (m/sec)
! M = convective mass flux (kg/(m2 sec)) 
! rho - density of air
! L = depth of convecting layer (m)
! c = normalization constant = diff_norm/g
!   (diff_norm is a namelist parameter;
!      the factor of g here is an historical artifact)
! for further discussion see cu_mo_trans.ps
!</PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans (is, js, Time, mass_flux, t,           &
!		p_half, p_full, z_half, z_full, diff)
!		
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! ! horizontal domain on which computation to be performed is
!    (is:is+size(t,1)-1,ie+size(t,2)-1) in global coordinates
!   (used by diag_manager only)
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! current time, used by diag_manager
!  </IN>
!  <IN NAME="mass_flux" TYPE="real">
! convective mass flux (Kg/(m**2 s)), dimension(:,:,:), 3rd dimension is
!    vertical level (top down) -- defined at interfaces, so that
!    size(mass_flux,3) = size(p_half,3); entire field processed;
!   all remaining fields are 3 dimensional;
!  size of first two dimensions must confrom for all variables
!  </IN>
!  <IN NAME="t" TYPE="real">
! temperature (K) at full levels, size(t,3) = size(p_full,3)
!  </IN>
!  <IN NAME="p_half" TYPE="real">
! pressure at interfaces (Pascals) 
!  size(p_half,3) = size(p_full,3) + 1
!  </IN>
!  <IN NAME="p_full" TYPE="real">
! pressure at full levels (levels at which temperature is defined)
!  </IN>
!  <IN NAME="z_half" TYPE="real">
! height at half levels (meters); size(z_half,3) = size(p_half,3)
!  </IN>
!  <IN NAME="z_full" TYPE="real">
! height at full levels (meters); size(z_full,3) = size(p_full,3)
!  </IN>
!  <OUT NAME="diff" TYPE="real">
! kinematic diffusivity (m*2/s); defined at half levels 
!   size(diff,3) = size(p_half,3)
!  </OUT>
! </SUBROUTINE>
!
subroutine cu_mo_trans (is, js, Time, mass_flux, t,           &
                        p_half, p_full, z_half, z_full, diff)

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: is, js

real, intent(in)   , dimension(:,:,:) :: mass_flux, t, &
                                         p_half, p_full, z_half, z_full
real, intent(out), dimension(:,:,:) :: diff                      

real, dimension(size(t,1),size(t,2),size(t,3)) :: rho
real, dimension(size(t,1),size(t,2))           :: zbot, ztop

integer :: k, nlev
logical :: used

!-----------------------------------------------------------------------

 if (.not.module_is_initialized) call error_mesg ('cu_mo_trans',  &
                      'cu_mo_trans_init has not been called.', FATAL)

!-----------------------------------------------------------------------

nlev = size(t,3)

zbot = z_half(:,:,nlev+1)
ztop = z_half(:,:,nlev+1)
  
do k = 2, nlev
  where(mass_flux(:,:,k) .ne. 0.0 .and. mass_flux(:,:,k+1) == 0.0) 
    zbot = z_half(:,:,k)
  endwhere
  where(mass_flux(:,:,k-1) == 0.0 .and. mass_flux(:,:,k) .ne. 0.0) 
    ztop = z_half(:,:,k)
  endwhere
end do

rho  = p_full/(RDGAS*t)  ! density 
   ! (including the virtual temperature effect here might give the 
   ! impression that this theory is accurate to 2%!)

! diffusivity = c*W*L
! W = M/rho  (m/sec)
! M = convective mass flux (kg/(m2 sec)) 
! L = ztop - zbot = depth of convecting layer (m)
! c = normalization constant = diff_norm/g
!   (the factor of g here is an historical artifact)

diff(:,:,1) = 0.0
do k = 2, nlev
  diff(:,:,k) = diff_norm*mass_flux(:,:,k)*(ztop-zbot)/(rho(:,:,k)*GRAV)
end do


! --- diagnostics
     if ( id_diff_cmt > 0 ) then
        used = send_data ( id_diff_cmt, diff, Time, is, js, 1 )
     endif

end subroutine cu_mo_trans


!#######################################################################


end module cu_mo_trans_mod

! <INFO>

! </INFO>

