module local_heating_mod


  use           fms_mod, only: error_mesg, FATAL, file_exist,       &
                               open_namelist_file, set_domain,      &
  			                       read_data, check_nml_error,          &
                               mpp_pe, mpp_root_pe, close_file,     &
                               write_version_number, stdlog,        &
                               uppercase,&  !pjk
                               mpp_clock_id,mpp_clock_begin,mpp_clock_end,CLOCK_COMPONENT!,mpp_chksum

  use diag_manager_mod, only: register_diag_field, send_data
  use time_manager_mod, only: time_type

  use  time_manager_mod, only: time_type, get_time, length_of_year

  use  diag_manager_mod, only: register_diag_field, send_data

  use  field_manager_mod, only: MODEL_ATMOS, parse

  use constants_mod, only   : RADIAN,PI
  
  use rrtm_astro,only : equinox_day

  implicit none
  
  !-----------------------------------------------------------------------
  !---------- interfaces ------------
  public :: local_heating,local_heating_init
  
  !---------------------------------------------------------------------------------------------------------------
  !
  !-------------------- diagnostics fields -------------------------------
  integer :: id_tdt_lheat
  character(len=14) :: mod_name = 'local_heating'
  real :: missing_value = -999.
  !-----------------------------------------------------------------------
  !-------------------- namelist -----------------------------------------
  !-----------------------------------------------------------------------
  integer,parameter :: ngauss = 10
  real,dimension(ngauss)   :: hamp      = 0.        ! heating amplitude [K/d]
  real,dimension(ngauss)   :: lonwidth  = -1.       ! zonal width of Gaussian heating [deg], zonally symmetric if <0
  real,dimension(ngauss)   :: loncenter = 90.       ! zonal center of Gaussian heating [deg]
  real,dimension(ngauss)   :: lonmove   = 0.        ! zonal center motion [deg/day]
  real,dimension(ngauss)   :: latwidth  = 15.       ! meridional width of Gaussian heating [deg]
  real,dimension(ngauss)   :: latcenter = 0.        ! meridional center of Gaussian heating [deg]
  real,dimension(ngauss)   :: latmove   = 0.        ! metidional center motion [deg/day]
  real,dimension(ngauss)   :: pwidth    = 2.        ! height of Gaussian heating in log-pressure [log10(hPa)]
  real,dimension(ngauss)   :: pcenter   = 1.        ! center of Gaussian heating in pressure [hPa]
  real,dimension(ngauss)   :: pmove     = 0.        ! vertical center motion [hPa/day]
  logical,dimension(ngauss):: is_periodic= .false.  ! if .true., reset location in accordance with temporal evolution
                                                    !  in this case, tphase and tperiod also apply to spatial position
  real,dimension(ngauss)   :: twidth    =-1.        ! temporal width of Gaussian heating [days],
                                                    !   constant in time if <0
  real,dimension(ngauss)   :: tphase    = 0.        ! temporal phase of Gaussian heating [days]
  real,dimension(ngauss)   :: tperiod   =-1.        ! temporal period of Gaussian heating
                                                    ! if < 0, period is in fraction of year
                                                    ! if > 0, period is in days
  
  namelist /local_heating_nml/ hamp \
                               ,lonwidth,loncenter,lonmove \
                               ,latwidth,latcenter,latmove \
                               ,pwidth,pcenter,pmove       \
                               ,is_periodic                \
                               ,twidth,tphase,tperiod
  
  
  ! local variables
  real,dimension(ngauss) :: logpc
  integer                :: daysperyear
  
contains
  
  subroutine local_heating_init(axes, Time)
    implicit none
    integer, intent(in), dimension(4) :: axes
    type(time_type), intent(in)       :: Time
    !-----------------------------------------------------------------------
    integer :: seconds
    integer :: unit, io, ierr, n

    !     ----- read namelist -----

    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read  (unit, nml=local_heating_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'local_heating_nml')
       enddo
10     call close_file (unit)
    endif
    
   ! ---- convert input units to code units  -----
    call get_time(length_of_year(),seconds,daysperyear)
    do n = 1,ngauss
       pcenter(n)   = pcenter(n)*100      ! convert hPa to Pa
       hamp(n)      = hamp(n)/86400.      ! convert K/d to K/s
       loncenter(n) = loncenter(n)/RADIAN ! convert degrees to radians
       lonwidth(n)  = lonwidth(n)/RADIAN  ! convert degrees to radians
       lonmove(n)   = lonmove(n)/RADIAN/86400. ! convert degrees/day to radians/s
       latcenter(n) = latcenter(n)/RADIAN ! convert degrees to radians
       latwidth(n)  = latwidth(n)/RADIAN  ! convert degrees to radians
       latmove(n)   = latmove(n)/RADIAN/86400. ! convert degrees/day to radians/s
       pmove(n)     = pmove(n)*100./86400.! convert hPa/day to Pa/s
       if ( tperiod(n) .lt. 0.0 ) tperiod(n) = -tperiod(n)*daysperyear ! convert year fraction to day of year
       twidth(n)    = twidth(n)*86400     ! convert to seconds
       tphase(n)    = tphase(n)*86400     ! convert to seconds
       tperiod(n)   = tperiod(n)*86400    ! convert to seconds
    enddo
  !----
  !------------ initialize diagnostic fields ---------------
    id_tdt_lheat = &
         register_diag_field ( mod_name, 'tdt_lheat', axes(1:3), Time, &
         'Temperature tendency due to local heating', &
         'K/s', missing_value=missing_value               )
    
  end subroutine local_heating_init
    
  !-----------------------------------------------------------------------
  !-------------------- computing localized heating ----------------------
  !-----------------------------------------------------------------------
  
  subroutine local_heating(is,js,Time,lon,lat,p_full,tdt_tot)
    implicit none
    integer, intent(in)                  :: is, js
    type(time_type),intent(in)           :: Time
    real, dimension(:,:)  ,intent(in)    :: lon,lat
    real, dimension(:,:,:),intent(in)    :: p_full
    real, dimension(:,:,:),intent(inout) :: tdt_tot
    ! local variables
    integer :: i,j,k,n
    real, dimension(size(lon,1),size(lon,2)) :: lon_factor
    real, dimension(size(lon,1),size(lon,2)) :: lat_factor
    real, dimension(size(tdt_tot,1),size(tdt_tot,2),size(tdt_tot,3)) :: tdt
    real    :: logp,p_factor,t_factor,targ,tcenter(3),halfper
    logical :: used
    integer :: seconds,days,fullseconds,startseconds,deltasecs
    
    call get_time(Time,seconds,days)
    fullseconds = days*86400+seconds
    

    tdt = 0.
    do n = 1,ngauss
       if ( hamp(n) .ne. 0. ) then
          ! local heating position is determined at peak heating time
          halfper = 0.5*tperiod(n)
          if ( is_periodic(n) ) then
             deltasecs = mod(fullseconds-tphase(n)+halfper,tperiod(n))-halfper
          else
             deltasecs = fullseconds
          endif
          ! compute the center position
          tcenter(1) =   mod(loncenter(n) + lonmove(n)*deltasecs,2*PI)
          tcenter(2) =       latcenter(n) + latmove(n)*deltasecs
          tcenter(3) = log10(pcenter  (n) + pmove  (n)*deltasecs)
          ! temporal component
          if (twidth(n) .lt. 0.0 ) then
             t_factor = 1.0
          else
             targ = mod(fullseconds-tphase(n)+halfper,tperiod(n))-halfper
             t_factor = exp( -(targ)**2/(2*(twidth(n))**2) )
          endif
          ! meridional and zonal components
          do j = 1,size(lon,2)
             do i = 1,size(lon,1)
                if ( loncenter(n) .ge. 0.0 ) then
                   lon_factor(i,j) = exp( -(lon(i,j)-tcenter(1))**2/(2*(lonwidth(n))**2) )
                   ! there is a problem when the heating is close to 360/0
                   lon_factor(i,j) = max(lon_factor(i,j), \
                                     exp( -(lon(i,j)+2*PI-tcenter(1))**2/(2*(lonwidth(n))**2) ) )
                   lon_factor(i,j) = max(lon_factor(i,j), \
                                     exp( -(lon(i,j)-2*PI-tcenter(1))**2/(2*(lonwidth(n))**2) ) )
                else
                   lon_factor(i,j) = 1.0
                endif
                lat_factor(i,j) = exp( -(lat(i,j)-tcenter(2))**2/(2*(latwidth(n))**2) )
             enddo
          enddo
    
          do k=1,size(p_full,3)
             do j = 1,size(lon,2)
                do i = 1,size(lon,1)
                   logp = log10(p_full(i,j,k))
                   ! vertical component
                   p_factor = exp(-(logp-tcenter(3))**2/(2*(pwidth(n))**2))
                   ! everything together
                   tdt(i,j,k) = tdt(i,j,k) + hamp(n)*t_factor*lon_factor(i,j)*lat_factor(i,j)*p_factor
                enddo
             enddo
          enddo
       endif
    enddo
     
    tdt_tot = tdt_tot + tdt
     
   !------- diagnostics ------------
    if ( id_tdt_lheat > 0 ) then
       used = send_data ( id_tdt_lheat, tdt, Time, is, js, 1 )
    endif

   end subroutine local_heating
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

end module local_heating_mod
