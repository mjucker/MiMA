
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

integer :: gulf_k         = 4       ! wave number of gulfstream perturbation []

real ::    warmpool_k     = 1,    & ! wave number of warmpool []
	   gulf_phase     = 300., & ! phase of warmpool [deg lon]
           gulf_amp       = 0.  , & ! amplitude of gulf stream perturbation [W/m2]
	   kuroshio_amp   = 0. , &  ! amplitude of kuroshio perturbation [W/m2]
	   trop_atlantic_amp  = 0.  ! amplitude of tropical atlantic perturbation [W/m2]

integer :: warmpool_localization_choice    = 1 ! 1->cos, 2->cos but restricted to Indo-Pacific
logical :: qflux_initialized = .false.


namelist /qflux_nml/ qflux_amp,qflux_width,&
                     warmpool_amp,warmpool_width,warmpool_centr,&
                     warmpool_k,warmpool_phase,warmpool_localization_choice,&
		     gulf_k,gulf_phase,gulf_amp,kuroshio_amp,trop_atlantic_amp

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
    real lon,lat,piphase, pigulfphase, latgulf

!    piphase = warmpool_phase/pi
    piphase = warmpool_phase*pi/180  !modified by cig, nov 15 2017
    pigulfphase = gulf_phase*pi/180  !modified by cig, nov 15 2017
    do j=1,size(latb)-1
       lat = 0.5*(latb(j+1) + latb(j))*180./pi
       latgulf = (lat-37)/10
       lat = (lat-warmpool_centr)/warmpool_width
      
       if( abs(lat) .le. 1.0 ) then
          do i=1,size(lonb)-1
             lon = 0.5*(lonb(i+1) + lonb(i))
 	     if (warmpool_localization_choice == 1 ) then
               flux(i,j) = flux(i,j) &
                    &+ (1.-lat**4.)*warmpool_amp*cos(warmpool_k*(lon-piphase))  !modified by cig, nov 15 2017, note that I use 4th power not 2nd as in MJ to better match Pac warm pool
	     elseif (warmpool_localization_choice == 2 ) then  !assumes k=5/3 for warmpool, 
		if( lon .ge. (warmpool_phase - 54)*pi/180. .AND. lon .lt. (warmpool_phase + 162)*pi/180.) then !modified by cig, nov 15 2017
                  flux(i,j) = flux(i,j) &
                     &+ (1.-lat**4.)*warmpool_amp*cos(warmpool_k*(lon-piphase))
	        endif !assumes k=4 for gulfstream-Africa perturbation
		if( lon .ge. (gulf_phase - 22)*pi/180. .AND. lon .lt. (gulf_phase + 68)*pi/180.) then !modified by cig, april 30 2018
                  flux(i,j) = flux(i,j) &
                     &+ (1.-lat**4.)*trop_atlantic_amp*cos(gulf_k*(lon-pigulfphase))
	        endif
             endif
	   !  write (*,*) "lon",lon," phiphase",piphase
          enddo
       endif
       if( abs(latgulf) .le. 1.0 .and. warmpool_localization_choice == 2) then
          do i=1,size(lonb)-1 !add Kuroshio and Gulf streams
             lon = 0.5*(lonb(i+1) + lonb(i))
 	         if( lon .ge. (gulf_phase - 42)*pi/180. .AND. lon .lt. (gulf_phase + 48)*pi/180.) then !modified by cig, june 3 2018
                  flux(i,j) = flux(i,j) &
                     &+ (1.-latgulf**4.)*gulf_amp*cos(gulf_k*(lon-(pigulfphase-20*pi/180.)))
	        endif
                if( lon .ge. (warmpool_phase - 17)*pi/180. .AND. lon .lt. (warmpool_phase + 73)*pi/180.) then !modified by cig, june 3 2018, use k=4 for kuroshio
                  flux(i,j) = flux(i,j) &
                     &+ (1.-latgulf**4.)*kuroshio_amp*cos(gulf_k*(lon-(piphase+5*pi/180)))
	        endif
              
	   !  write (*,*) "lon",lon," phiphase",piphase
          enddo
       endif
    enddo

  end subroutine warmpool

!########################################################

end module qflux_mod

