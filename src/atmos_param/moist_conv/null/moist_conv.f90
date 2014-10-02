
module moist_conv_mod

!-----------------------------------------------------------------------

 use           mpp_mod, only : mpp_pe,             &
                               mpp_root_pe,        &
                               stdlog
use   time_manager_mod, only : time_type
 use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
use  sat_vapor_pres_mod, ONLY: EsComp, DEsComp
use             fms_mod, ONLY:  error_mesg, file_exist, open_namelist_file,  &
                                check_nml_error, close_file,        &
                                FATAL, WARNING, NOTE, mpp_pe, mpp_root_pe, &
                                write_version_number, stdlog
use       constants_mod, ONLY: HLv, HLs, cp_air, grav, rdgas, rvgas

use           fms_mod, only : write_version_number, ERROR_MESG, FATAL
use field_manager_mod, only : MODEL_ATMOS
use tracer_manager_mod, only : get_tracer_index,   &
                               get_number_tracers, &
                               get_tracer_names,   &
                               get_tracer_indices, &
                               query_method,       &
                               NO_TRACER

implicit none
private

!------- interfaces in this module ------------

public :: moist_conv, moist_conv_init, moist_conv_end


!-----------------------------------------------------------------------
!---- VERSION NUMBER -----

 character(len=128) :: version = '$Id: moist_conv.f90,v 11.0 2004/09/28 19:19:48 fms Exp $'
 character(len=128) :: tagname = '$Name: lima $'
 logical            :: module_is_initialized = .false.

!-----------------------------------------------------------------------

CONTAINS

!#######################################################################

 subroutine moist_conv ( Tin, Qin, Pfull, Phalf, coldT,    &
                         Tdel, Qdel, Rain, Snow, Lbot, &
                         do_strat, ql, qi, cf,  &
                         qldel, qidel, cfdel,   &
                         dtinv, Time, mask, is, js, Conv, &
                         tracers, qtrmca )

!-----------------------------------------------------------------------
!
!                       MOIST CONVECTIVE ADJUSTMENT
!
!-----------------------------------------------------------------------
!
!   INPUT:   Tin     temperature at full model levels
!            Qin     specific humidity of water vapor at full
!                      model levels
!            Pfull   pressure at full model levels
!            Phalf   pressure at half model levels
!            coldT   Should MCA produce snow in this column?
!
!   OUTPUT:  Tdel    temperature adjustment at full model levels (deg k)
!            Qdel    specific humidity adjustment of water vapor at
!                       full model levels
!            Rain    liquid precipitiation (in Kg m-2)
!            Snow    ice phase precipitation (kg m-2)
!  OPTIONAL
!
!   INPUT:   Lbot    integer index of the lowest model level,
!                      Lbot is always <= size(Tin,3)
!              cf    stratiform cloud fraction (used only when
!                    operating with stratiform cloud scheme) (fraction)
!
!  OUTPUT:   Conv    logical flag; TRUE then moist convective
!                       adjustment was performed at that model level.
!            cfdel   change in stratiform cloud fraction (fraction)
!            qldel   change in liquid water condensate due to
!                    convective detrainment (kg condensate /kg air)
!            qidel   change in ice condensate due to
!                    convective detrainment (kg condensate /kg air)
!
!-----------------------------------------------------------------------
!----------------------PUBLIC INTERFACE ARRAYS--------------------------
    real, intent(INOUT) , dimension(:,:,:)          :: Tin, Qin
    real, intent(IN) , dimension(:,:,:)             :: Pfull, Phalf 
 logical, intent(IN) , dimension(:,:)               :: coldT
    real, intent(OUT), dimension(:,:,:)             :: Tdel, Qdel
    real, intent(OUT), dimension(:,:)               :: Rain, Snow
 integer, intent(IN) , dimension(:,:)  , optional   :: Lbot
 logical, intent(OUT), dimension(:,:,:), optional   :: Conv
logical, intent(in)                                 :: do_strat
    real, intent(IN)                                :: dtinv
 integer, intent(IN)                                :: is, js     
    real, intent(INOUT), dimension(:,:,:)           :: ql, qi, cf
    real, intent(IN), dimension(:,:,:,:), optional  :: tracers
    real, dimension(:,:,:,:), intent(out), optional :: qtrmca
    real, intent(INOUT), dimension(:,:,:)           :: qldel, qidel, cfdel
type(time_type), intent(in)                         :: Time
   real, intent(in) , dimension(:,:,:), optional    :: mask
         

!---------------------------------------------------------------------

      call error_mesg('moist_conv', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine moist_conv

!#######################################################################

 subroutine moist_conv_init (axes, Time, tracers_in_mca)

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 logical, dimension(:), intent(in), optional :: tracers_in_mca

!-----------------------------------------------------------------------
      

!---------- output namelist --------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number (version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('moist_conv_init', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine moist_conv_init


!#######################################################################
 subroutine moist_conv_end

     module_is_initialized = .FALSE.

!---------------------------------------------------------------------

      call error_mesg('moist_conv_end', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine moist_conv_end

!#######################################################################

end module moist_conv_mod

