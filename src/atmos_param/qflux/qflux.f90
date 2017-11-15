
module qflux_mod

  use constants_mod, only : pi
  use       fms_mod, only : file_exist, open_namelist_file, check_nml_error, &
                            error_mesg, FATAL, close_file

implicit none

real ::    qflux_amp      = 30.,  & ! amplitude of meridional Q-flux [W/m2]
           qflux_width    = 16.,  & ! half-width of Q-flux [deg lat]
           warmpool_amp   =  5.,  & ! amplitude of warmpool [W/m2]
           warmpool_width = 20.,  & ! width of warmpool (square profile) [deg lat]
           warmpool_centr =  0.,  & ! center of warmpool [deg lat]
           warmpool_phase = 0.0     ! phase of warmpool [deg lon]
integer :: warmpool_k     = 1       ! wave number of warmpool []

logical :: qflux_initialized = .false.


namelist /qflux_nml/ qflux_amp,qflux_width,&
                     warmpool_amp,warmpool_width,warmpool_centr,&
                     warmpool_k,warmpool_phase

private 

public :: qflux_init,qflux,warmpool

contains

!########################################################
  subroutine qflux_init
    implicit none
    integer :: unit, ierr, io

    if ( file_exist('input.nml') )then
       unit = open_namelist_file()
       ierr=1; 
       do while (ierr /= 0)
          read( unit, nml=qflux_nml, iostat=io, end=10 )
          ierr = check_nml_error(io,'qflux_nml')
       enddo
10     call close_file(unit)
    endif

    qflux_initialized = .true.
    
  end subroutine qflux_init
!########################################################

  subroutine qflux(latb, flux)
! compute Q-flux as in Merlis et al (2013) [Part II]
    implicit none
    real,dimension(:)  ,intent(in)    :: latb   !latitude boundary
    real,dimension(:,:),intent(inout) :: flux   !total ocean heat flux
!
    integer j
    real lat,coslat

    if( .not. qflux_initialized ) then
       call error_mesg('qflux','qflux module not initialized',FATAL)
    endif
    
    do j=1, size(latb)-1
       lat = 0.5*(latb(j+1) + latb(j))
       coslat = cos(lat)
       lat = lat*180./pi
       flux(:,j) = flux(:,j) - qflux_amp*(1-2.*lat**2/qflux_width**2) * &
            exp(- ((lat)**2/(qflux_width)**2))/coslat
    enddo

  end subroutine qflux

!########################################################

  subroutine warmpool(lonb, latb, flux)
    implicit none
    real,dimension(:)  ,intent(in)   :: lonb,latb  !lon and lat boundaries
    real,dimension(:,:),intent(inout):: flux       !total ocean heat flux
!
    integer i,j
    real lon,lat,piphase

    piphase = warmpool_phase*pi/180.
    do j=1,size(latb)-1
       lat = 0.5*(latb(j+1) + latb(j))*180./pi
       lat = (lat-warmpool_centr)/warmpool_width
       if( abs(lat) .le. 1.0 ) then
          do i=1,size(lonb)-1
             lon = 0.5*(lonb(i+1) + lonb(i))
             flux(i,j) = flux(i,j) &
                  &+ (1.-lat**2.)*warmpool_amp*cos(warmpool_k*lon+piphase)
          enddo
       endif
    enddo

  end subroutine warmpool

!########################################################

end module qflux_mod
