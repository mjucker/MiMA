
                    module longwave_params_mod

! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code that contains parameters for the longwave code
! </OVERVIEW>  
! <DESCRIPTION>
!  This code has the number of bands defined for the longwave gases.
! </DESCRIPTION>
! 

!   shared modules:

 use fms_mod,               only: open_namelist_file, fms_init, &
                                  mpp_pe, mpp_root_pe, stdlog, &
                                  file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, NOTE, WARNING, close_file

!--------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!    longwave_params_mod defines basic parameters used by the
!    longwave radiation code.
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: longwave_params.f90,v 10.0 2003/10/24 22:00:42 fms Exp $'
character(len=128)  :: tagname =  '$Name: lima $'


!--------------------------------------------------------------------
!----- interfaces ------

public     &
         longwave_params_init, &
         longwave_params_end

!private   &


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=8)   :: dummy  = '     '


namelist /longwave_params_nml/    &
                               dummy                 

!-------------------------------------------------------------------
!----- public data --------

!--------------------------------------------------------------------
!       NBCO215
!       NBLY_RSB
!       NBLY_CKD
!       NBLW
!       NBLX
!---------------------------------------------------------------------
integer, parameter, public   :: NBCO215     = 3
integer, parameter, public   :: NBLY_RSB    = 16
integer, parameter, public   :: NBLY_CKD    = 48
integer, parameter, public   :: NBLW        = 300
integer, parameter, public   :: NBLX        = 48




!-------------------------------------------------------------------
!----- private data --------

logical :: module_is_initialized = .false.  ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!####################################################################

! <SUBROUTINE NAME="longwave_params_init">
!  <OVERVIEW>
!   Subroutine to initialize longwave parameter module
!  </OVERVIEW>
!  <DESCRIPTION>
!   This is the longwave_params constructor method
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_params_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_params_init

!------------------------------------------------------------------
!    longwave_params_init is the constructor for longwave_params_mod.
!------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      integer    ::  unit, ierr, io

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
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

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=longwave_params_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'longwave_params_nml')
        end do
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                          write (stdlog(), nml=longwave_params_nml)
      if (mpp_pe() == mpp_root_pe() ) &
           write (stdlog(),9000) NBCO215, NBLY_RSB, NBLY_CKD,   &
                                 NBLW, NBLX 

!----------------------------------------------------------------------
!    mark the module as initialized.
!----------------------------------------------------------------------
     module_is_initialized = .true.

!------------------------------------------------------------------
9000 format ( 'NBCO215=', i3,'  NBLY_RSB=', i4,   &
              '  NBLY_CKD=', i4, '  NBLW= ', i4, '  NBLX=', i4 )

!-------------------------------------------------------------------


end  subroutine longwave_params_init



!####################################################################

! <SUBROUTINE NAME="longwave_params_end">
!  <OVERVIEW>
!   Subroutine to close out longwave parameter module
!  </OVERVIEW>
!  <DESCRIPTION>
!   This is the longwave_params destructor subroutine
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_params_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_params_end

!-------------------------------------------------------------------
!    longwave_params_end is the destructor for longwave_params_mod.
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_params_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------



end subroutine longwave_params_end


!####################################################################


                   end module longwave_params_mod
