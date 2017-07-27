module local_heating_mod
  #ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
  #else
  use fms_mod, only: open_namelist_file
  #endif

  use           fms_mod, only: error_mesg, FATAL, file_exist,       &
                               open_namelist_file, set_domain,      &
  			                       read_data, check_nml_error,          &
                               mpp_pe, mpp_root_pe, close_file,     &
                               write_version_number, stdlog,        &
                               uppercase,&  !pjk
                               mpp_clock_id,mpp_clock_begin,mpp_clock_end,CLOCK_COMPONENT!,mpp_chksum

  use  time_manager_mod, only: time_type, get_time

  use  diag_manager_mod, only: register_diag_field, send_data

  use  field_manager_mod, only: MODEL_ATMOS, parse
  implicit none
  
  !-----------------------------------------------------------------------
  !---------- interfaces ------------
  public :: local_heating
  
  !-----------------------------------------------------------------------
  !-------------------- namelist -----------------------------------------
  !-----------------------------------------------------------------------
  real :: hamp = -10.       ! heating amplitude [K/s]
  real :: latwidth = 15.       ! latitudinal width of Gaussian heating [deg]
  real :: latcenter = 60.      ! latitudinal center of Gaussian heating [deg]
  real :: pwidth = 1.         ! height of Gaussian heating in log-pressure [log(hPa)]
  real :: pcenter = 10.        ! center of Gaussian heating in pressure [hPa]
  
  namelist /local_heating_nml/ hamp,latwidth,latcenter,pwdith,pcenter
  
  
  ! local variables
  real :: logpc
  character(len=14) :: mod_name = 'local_heating'
  
contains
  
  subroutine local_heating_init()
    implicit none
    !-----------------------------------------------------------------------
    integer  unit, io, ierr

    !     ----- read namelist -----

    #ifdef INTERNAL_FILE_NML
         read (input_nml_file, nml=local_heating_nml, iostat=io)
         ierr = check_nml_error(io, 'local_heating_nml')
    #else
          if (file_exist('input.nml')) then
             unit = open_namelist_file ( )
             ierr=1; do while (ierr /= 0)
                read  (unit, nml=local_heating_nml, iostat=io, end=10)
                ierr = check_nml_error (io, 'local_heating_nml')
             enddo
      10     call close_file (unit)
          endif
    #endif
    
    logpc = log(pcenter)
  end subroutine local_heating_init
    
  !-------------------- computing localized heating ----------------------
  !-----------------------------------------------------------------------
  
  subroutine local_heating(lon,lat,p_full,tdt_tot)
    implicit none
    real, dimension(:,:)  ,intent(in)    :: lon,lat
    real, dimension(:,:,:),intent(in)    :: p_full
    real, dimension(:,:,:),intent(inout) :: tdt_tot
    ! local variables
    integer :: i,j,k
    real, dimension(size(lon,1),size(lon,2)) :: lat_factor
    real, dimension(size(tdt_tot,1),size(tdt_tot,2),size(tdt_tot,3)) :: tdt
    real    :: logp,logpc,p_factor,lat_factor
    
    do j = 1,size(lon,2)
      do i = 1,size(lon,1)
        lat_factor(i,j) = exp( -(lat(i,j)-latcenter)**2/(2*(latwidth)**2) )
      enddo
    enddo
  
    do k=1,size(p_full,3)
      do j = 1,size(lon,2)
        do i = 1,size(lon,1)
          logp = log(p_full(i,j,k))
          p_factor = exp(-(logp-logpc)**2/(2*(pwidth)**2))
          tdt(i,j,k) = hamp*lat_factor(i,j)*p_factor
        enddo
      enddo
    enddo
     
     tdt_tot = tdt_tot + tdt
     
   end subroutine local_heating