
                 module donner_deep_clouds_W_mod

use time_manager_mod,       only: time_type
use       fms_mod,          only: error_mesg, FATAL, & 
                                  mpp_pe, mpp_root_pe, &
                                  write_version_number
use rad_utilities_mod,      only: microphysics_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!          donner deep cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: donner_deep_clouds_W.f90,v 12.0 2005/04/14 15:49:35 fms Exp $'
   character(len=128)  :: tagname =  '$Name: lima $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          donner_deep_clouds_W_init, donner_deep_clouds_calc,  &
          donner_deep_clouds_W_end , donner_deep_clouds_amt


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


  logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------

contains 

subroutine donner_deep_clouds_W_init  (pref, lonb, latb, axes, Time)

real, dimension(:,:), intent(in) :: pref
real, dimension(:), intent(in) :: lonb, latb
integer, dimension(4), intent(in)      :: axes
type(time_type),       intent(in)      :: Time

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

!---------------------------------------------------------------------

       module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('donner_deep_clouds_W_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_clouds_W_init

subroutine donner_deep_clouds_W_end
       
!----------------------------------------------------------------------
!    diag_clouds_end is the destructor for diag_clouds_W_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!---------------------------------------------------------------------

      call error_mesg('donner_deep_clouds_W_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_clouds_W_end


!#################################################################

subroutine donner_deep_clouds_amt (is, ie, js, je, Cell_microphys,  &
                                   Meso_microphys)

!---------------------------------------------------------------------
!    donner_deep_clouds_amt defines the distribution of cloud water and
!    cloud ice concentration and particle size and total cloud fraction
!    in both the mesoscale and convective cell-scale components of the
!    clouds associated with donner_deep convection. these values will
!    be combined with the large-scale cloud fields to produce the dist-
!    ribution of cloud radiative properties that will be seen by the
!    radiation package.
!----------------------------------------------------------------------

integer,                 intent(in)    :: is,ie,js,je
type(microphysics_type), intent(inout) :: Cell_microphys, Meso_microphys

!---------------------------------------------------------------------

      call error_mesg('donner_deep_clouds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_clouds_amt  


!#####################################################################


subroutine donner_deep_clouds_calc (             &
                  is,ie,js,je,deltaz,press,temp,                 &
                                    cld_cell,               &
 cldext_cell, cldsct_cell, cldasymm_cell,  &
    abscoeff_cell,          &
                                    cld_meso,               &
 cldext_meso, cldsct_meso, cldasymm_meso,  &
    abscoeff_meso)


integer, intent(in) :: is,ie,js,je
real, dimension(:,:,:), intent(in) :: deltaz, press, temp
real, dimension(:,:,:), intent(inout) :: cld_cell, cld_meso
real, dimension(:,:,:,:), intent(out) :: cldext_cell, cldsct_cell,  &
                                         cldasymm_cell, abscoeff_cell
real, dimension(:,:,:,:), intent(out) :: cldext_meso, cldsct_meso,  &
                                         cldasymm_meso, abscoeff_meso


!---------------------------------------------------------------------

      call error_mesg('donner_deep_clouds_calc', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_clouds_calc


!####################################################################


       end module donner_deep_clouds_W_mod



