
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
	   gulf_phase     = 140., & ! phase of warmpool [deg lon]
           gulf_amp       = 0.  , & ! amplitude of gulf stream perturbation [W/m2]
	   kuroshio_amp   = 0. , &  ! amplitude of kuroshio perturbation [W/m2]
	   trop_atlantic_amp  = 0. , &  ! amplitude of tropical atlantic perturbation [W/m2]
           north_sea_heat = 0., & !add extra perturbation to move heat from Canada to North Sea
	   Pac_ITCZextra = 0., & !extra q flux in tropical South Pacific to strengthen local ITCZ	 
	   Pac_SPCZextra = 0., & !extra q flux in subtropical pacific to modulate SPCZ
	   Africaextra = 0.     !extra q flux by Agulhaus
	
integer :: warmpool_localization_choice    = 1 ! 1->Jucker and Gerber 2017. 2->Garfinkel et al 2020, NH stationary waves. 3->Garfinkel et al 2021, SH biases
logical :: qflux_initialized = .false.


namelist /qflux_nml/ qflux_amp,qflux_width,&
                     warmpool_amp,warmpool_width,warmpool_centr,&
                     warmpool_k,warmpool_phase,warmpool_localization_choice,&
		     gulf_k,gulf_phase,gulf_amp,kuroshio_amp,trop_atlantic_amp,&
		     north_sea_heat,  Pac_ITCZextra, &
		        Pac_SPCZextra,  Africaextra 

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
    real lon,lat,piphase, pigulfphase, latgulf,latgreen,  latorig, africaamp

	africaamp=trop_atlantic_amp*2/4.
!    piphase = warmpool_phase/pi
    piphase = warmpool_phase*pi/180  !modified by cig, nov 15 2017
    pigulfphase = gulf_phase*pi/180  !modified by cig, nov 15 2017
    do j=1,size(latb)-1
       lat = 0.5*(latb(j+1) + latb(j))*180./pi
       latgulf = (lat-37.)/10.
       latgreen = (lat - 67.)/10.

	
        latorig=lat   
	lat = (lat-warmpool_centr)/warmpool_width
       
	 
       
          do i=1,size(lonb)-1
             lon = 0.5*(lonb(i+1) + lonb(i))
	     if (abs(lat) .le. 1.0 ) then
	 	     if (warmpool_localization_choice == 1 ) then
        	       flux(i,j) = flux(i,j) &
        	            &+ (1.-lat**4.)*warmpool_amp*cos(warmpool_k*(lon-piphase))  !modified by cig, nov 15 2017, note that I use 4th power not 2nd as in MJ to better match Pac warm pool
		     elseif (warmpool_localization_choice == 2 .or. warmpool_localization_choice == 3  ) then  !assumes k=5/3 for warmpool, 
			if( lon .ge. (warmpool_phase - 54)*pi/180. .AND. lon .le. (warmpool_phase + 162)*pi/180.) then !modified by cig, nov 15 2017
        	          flux(i,j) = flux(i,j) &
        	             &+ (1.-lat**4.)*warmpool_amp*cos(warmpool_k*(lon-piphase))
		        endif 
			if( lon .ge. (warmpool_phase + 117)*pi/180. .AND. lon .le. (warmpool_phase + 162)*pi/180. .AND. warmpool_localization_choice == 3  ) then !modified by cig, mar 28 2019
        	          flux(i,j) = flux(i,j) &
        	             &+ (1.-lat**4.)*warmpool_amp*sin(8*(lon-piphase-139.5*pi/180.))
		        endif 
		     endif
		     if( lon .ge. (warmpool_phase -130)*pi/180. .AND. lon .le. (warmpool_phase - 58)*pi/180. .and. warmpool_localization_choice == 3) then !modified by cig, may 13 2019
        	          flux(i,j) = flux(i,j) &
        	             &+ (1.-lat**2.)*africaamp*cos(5*(lon-(piphase-112*pi/180)))
		     endif
		     if( (lon .ge. (gulf_phase - 22)*pi/180. .OR. lon .le. (gulf_phase + 68 - 360)*pi/180.) .and. warmpool_localization_choice == 2) then !modified by cig, april 30 2018
        	          flux(i,j) = flux(i,j) &
        	             &+ (1.-lat**4.)*trop_atlantic_amp*cos(gulf_k*(lon-pigulfphase))
		     endif
              endif     
	   
        
      	    if( abs(latgulf) .le. 1.0 .and. warmpool_localization_choice == 2) then !add Kuroshio and Gulf streams
          
 	         if( lon .ge. (gulf_phase - 42)*pi/180. .AND. lon .lt. (gulf_phase + 48)*pi/180.) then !modified by cig, june 3 2018
                  flux(i,j) = flux(i,j) &
                     &+ (1.-latgulf**4.)*gulf_amp*cos(gulf_k*(lon-(pigulfphase-19.5*pi/180.)))
		 endif	
		 if( lon .ge. (gulf_phase - 42)*pi/180. .AND. lon .lt. (gulf_phase +3)*pi/180.) then !modified by cig, Nov 18 2018 to localize gulfstream more over ocean
                  flux(i,j) = flux(i,j) &	                    
		     &+ 0.535*(1.-latgulf**4.)*gulf_amp*sin(gulf_k*2*(lon-(pigulfphase-19.5*pi/180.)))
	        endif
		
                if( lon .ge. (warmpool_phase - 17.5)*pi/180. .AND. lon .lt. (warmpool_phase + 72.5)*pi/180.) then !modified by cig, june 3 2018, use k=4 for kuroshio
                  flux(i,j) = flux(i,j) &
                     &+ (1.-latgulf**2.)*kuroshio_amp*cos(gulf_k*(lon-(piphase+5*pi/180))) !modified by cig, mar 29 2019,
	        endif
  		 if( lon .ge. (warmpool_phase +30)*pi/180. .AND. lon .lt. (warmpool_phase + 90)*pi/180.) then !modified by cig, jan 2 2019, shift cooling kuroshio east k=6
                  flux(i,j) = flux(i,j) &
                     &- (1.-latgulf**2.)*0.65*kuroshio_amp*cos((2*gulf_k-2)*(lon-(piphase+15*pi/180)))!modified by cig, mar 29 2019
	        endif           
       endif
 	
       if( abs(latgreen) .le. 1.0 .and. warmpool_localization_choice == 2 .and. north_sea_heat .gt. 0.001) then
	    
 	       if( lon .ge. (gulf_phase - 52)*pi/180. .OR. lon .lt. (gulf_phase + 68 - 360)*pi/180.) then !modified by cig, Nov 22 2018 to add opposite of Gulfstream further polewrd
                  flux(i,j) = flux(i,j) &	                    
		     & +north_sea_heat*(1.-latgreen**4.)*gulf_amp*cos((gulf_k-1)*(lon-pigulfphase-38*pi/180.))
		 
	      endif
	       if( lon .ge. (gulf_phase - 67)*pi/180. .AND. lon .lt. (gulf_phase - 7)*pi/180.) then !modified by cig, Nov 22 2018 to smear out the cooling over Northern Canada over broad region
                  flux(i,j) = flux(i,j) &	                    
		     & +north_sea_heat*0.25*(1.-latgreen**4.)*gulf_amp*cos(2*(gulf_k-1)*(lon-pigulfphase+22*pi/180.))
		 
	      endif
	endif


   

       if(  warmpool_localization_choice == 3 .and. kuroshio_amp .gt. 0.001 .AND. latorig .le. 47. .AND. latorig .ge. 5. ) then !add Kuroshio  
               if( lon .ge. (warmpool_phase - 30.)*pi/180. .AND. lon .lt. (warmpool_phase + 130.)*pi/180.) then !modified by cig, may 13 2019, Pacific sector is exp 
                  flux(i,j) = flux(i,j) & 
                     &- kuroshio_amp*59.4/100.*exp(-(lon*180./pi+latorig-268.)**2./(2*49.))*exp(-(lon*180./pi-latorig-215.)**2./(2*625.))  &	       
                     &+ kuroshio_amp*exp(-(lon*180./pi-3*latorig-45.)**2./(2*100.))*exp(-(lon*180./pi+latorig-170.)**2./(2*400.))
	        endif

	endif

   

       if(  warmpool_localization_choice == 3 .and. warmpool_amp .gt. 0.001 .and.   lon .ge. 70.*pi/180. .AND. latorig .le. 60. .AND. latorig .ge. -10. .AND. lon .lt. 240.*pi/180.) then !add more Kurishio at expense of Southeast Asia  modified by cig, may 14 2019
	     if (kuroshio_amp .gt. 0.001) then
	     	flux(i,j) = flux(i,j) & 
                     &- 27.60*exp(-(lon*180./pi-140)**2./(2*1521))*(exp(-(latorig-19.7)**2./(2*49)))   &
		     &- (5.2)*exp(-(lon*180./pi-140)**2./(2*64))*(exp(-(latorig-20.)**2./(2*16)))   &
                     &+ 35.4*exp(-(lon*180./pi-160)**2./(2*400))*exp(-(latorig-35)**2./(2*36)) &
	   	     &+ 49.5*exp(-(lon*180./pi-3*latorig-45.)**2./(2*100.))*exp(-(lon*180./pi+latorig-160.)**2./(2*400.)) &
		     &+ 22.9*exp(-(lon*180./pi-90)**2./(2*144))*exp(-(latorig-0)**2./(2*25)) 
	     else 
		flux(i,j) = flux(i,j) & 
                     &- 5.3*exp(-(lon*180./pi-140)**2./(2*1521))*(exp(-(latorig-19.7)**2./(2*49)))   &
                     &+ 22.9*exp(-(lon*180./pi-90)**2./(2*144))*exp(-(latorig-0)**2./(2*25)) 
	     
	     endif
	endif 

  if(  warmpool_localization_choice == 3 .and. warmpool_amp .gt. 0.001 .AND. latorig .le. 24. .AND. latorig .ge. -78. .AND. lon .ge. 129.*pi/180. .AND. lon .lt. 290.*pi/180.) then !add south pacific and SPCZ at expense of further cold tongue  modified by cig, may 28 2019, also includes Australia bit from Kuroshio above
		 flux(i,j) = flux(i,j) & 
                     &-(50. -Pac_SPCZextra*.28)*exp(-(lon*180./pi-270)**2./(2*81))*(exp(-(latorig+0.)**2./(2*9)))   &
		     &-(50. -Pac_SPCZextra*.28)*exp(-(lon*180./pi-250)**2./(2*81))*(exp(-(latorig+1.)**2./(2*9)))   &
		     &-(50.-Pac_SPCZextra*.28)*exp(-(lon*180./pi-230)**2./(2*81))*(exp(-(latorig+2.)**2./(2*9)))   &
		     &-39.*exp(-(lon*180./pi-210)**2./(2*81))*(exp(-(latorig+2.)**2./(2*9)))   &
		     & -36.*exp(-(lon*180./pi-190)**2./(2*81))*(exp(-(latorig+0.)**2./(2*9))) &		  
		     & -16.*exp(-(lon*180./pi-170)**2./(2*81))*(exp(-(latorig+0.)**2./(2*9))) &		  
		     &-40.*exp(-(lon*180./pi-287)**2./(2*4))*(exp(-(latorig+25.)**2./(2*81)))   &
		     &-15.*exp(-(lon*180./pi-282)**2./(2*25))*(exp(-(latorig+15.)**2./(2*81)))   &
		     &-(25.+ Pac_ITCZextra+Pac_SPCZextra)*exp(-(lon*180./pi-240)**2./(2*1600))*(exp(-(latorig+21.)**2./(2*121)))   &	 	     
		     &-38.0*exp(-(lon*180./pi-195)**2./(2*169))*(exp(-(latorig-16.)**2./(2*49)))   &
		     &-51.4*exp(-(lon*180./pi-225)**2./(2*169))*(exp(-(latorig-16.)**2./(2*49)))   &
                     &+ (28.2+ Pac_ITCZextra*.8623)*exp(-(lon*180./pi-220)**2./(2*1600))*exp(-(latorig+57)**2./(2*225)) &
		     &+ (14.+Pac_SPCZextra*1.1195)*exp(-(lon*180./pi-165)**2./(2*400))*exp(-(latorig+20)**2./(2*25))&
		     &+ (16.+Pac_SPCZextra*1.1195)*exp(-(lon*180./pi-195)**2./(2*400))*exp(-(latorig+33)**2./(2*49)) &		   
		     &+ (50.+Pac_SPCZextra*1.1195)*exp(-(lon*180./pi-155)**2./(2*9))*exp(-(latorig+30)**2./(2*49)) &
		     &+ (40.+Pac_SPCZextra*1.1195)*exp(-(lon*180./pi-180)**2./(2*25))*exp(-(latorig+40)**2./(2*25))&
		     &+ (41.+ Pac_ITCZextra)*exp(-(lon*180./pi-240)**2./(2*900))*exp(-(latorig+62)**2./(2*64)) &
  		     &+ 60.*exp(-(lon*180./pi-180)**2./(2*169))*(exp(-(latorig-6.97)**2./(2*4))) &   
		     &+ 47.*exp(-(lon*180./pi-210)**2./(2*169))*(exp(-(latorig-6.97)**2./(2*4))) &   
		     &+ 45.*exp(-(lon*180./pi-240)**2./(2*169))*(exp(-(latorig-6.97)**2./(2*4))) &
		     &+ (19.5+Pac_SPCZextra)*exp(-(lon*180./pi-145)**2./(2*196))*(exp(-(latorig-3.)**2./(2*16)))  &
		     &+ (40.+Pac_SPCZextra*.435)*exp(-(lon*180./pi-150)**2./(2*169))*(exp(-(latorig-7.)**2./(2*9)))    
	endif



 if(  warmpool_localization_choice == 3 .and. warmpool_amp .gt. 0.001 .AND. latorig .le. 10. .AND. latorig .ge. -36. .AND. lon .ge. 50.*pi/180. .AND. lon .lt. 220.*pi/180.) then !add south pacific and south indian at expense of Australia  modified by cig, may 13 2019
		 flux(i,j) = flux(i,j) & 
                     &-(qflux_amp+warmpool_amp)*1.02*exp(-(lon*180./pi-135)**2./(2*225))*(exp(-(latorig+20.)**2./(2*36)))   &
                     &-(qflux_amp+warmpool_amp)*1.02*10/44*exp(-(lon*180./pi-147)**2./(2*64))*(exp(-(latorig+27.)**2./(2*49)))   &
                     &+ (qflux_amp+warmpool_amp)*1.02*16.6/44*exp(-(lon*180./pi-120)**2./(2*900))*exp(-(latorig+20)**2./(2*36)) &
		     &+ (qflux_amp+warmpool_amp)*1.02*(27.89)/44*exp(-(lon*180./pi-100)**2./(2*100))*exp(-(latorig+10)**2./(2*16)) &
		      &+ (qflux_amp+warmpool_amp)*1.02*4.9/44*exp(-(lon*180./pi-135)**2./(2*225))*exp(-(latorig+0)**2./(2*16)) 
		   
	endif


        if(  warmpool_localization_choice == 3 .and. gulf_amp .gt. 0.001 .AND. latorig .le. 52. .AND. latorig .ge. 10. ) then !add  Gulf streams
               if( lon .ge. (warmpool_phase + 135.)*pi/180. .AND. lon .lt. (warmpool_phase + 195.)*pi/180.) then !modified by cig, may 13 2019, Atlantic sector is exp 
                  flux(i,j) = flux(i,j) & 
                     &+ gulf_amp*exp(-(lon*180./pi-2.*latorig-220.)**2./(2*9.))*exp(-(lon*180./pi+latorig-335.)**2./(2*625.))
	        endif
  		if( lon .ge. (warmpool_phase +158.)*pi/180. .AND. lon .lt. (warmpool_phase + 218.)*pi/180.) then !modified by cig, may 13 2019
                  flux(i,j) = flux(i,j) & 
                     &- gulf_amp*63.9/70.*exp(-(lon*180./pi-.5*latorig-325.)**2./(2*9.))*exp(-(latorig-25.)**2./(2*49.))
	        endif

	endif

 	if(  warmpool_localization_choice == 3 .and. trop_atlantic_amp .gt. 0.001 .AND. latorig .le. 77. .AND. latorig .ge. -35. .AND. lon .ge. 275.*pi/180. .AND. gulf_amp .gt. 0.001) then !add more Gulf at expense of tropical Atlantic  modified by cig, may 14 2019
		 flux(i,j) = flux(i,j) & 
                     &- trop_atlantic_amp*(exp(-(lon*180./pi-342)**2./(2*81))+exp(-(lon*180./pi-360.)**2./(2*64)))*(exp(-(latorig+5.)**2./(2*25)))   &
	             &- 12.6*(exp(-(lon*180./pi-345)**2./(2*256)))*(exp(-(latorig+16.)**2./(2*64)))   &
                     &+ trop_atlantic_amp*30.65/28.*exp(-(lon*180./pi-2*latorig-220)**2./(2*100))*exp(-(latorig+lon*180./pi-375)**2./(2*900))
	endif  
       if( warmpool_localization_choice == 3 .and. trop_atlantic_amp .gt. 0.001  .AND. gulf_amp .gt. 0.001) then
	      if( lon .ge. (318.)*pi/180. .OR. lon .lt. (18.)*pi/180.) then !second part of Gulf at expense of South America, replaced greenland perturbation
		 if (abs(latgreen) .le. 1.0) then
                   flux(i,j) = flux(i,j) &	                    
		      & +trop_atlantic_amp*36./28*(1.-latgreen**4.)*cos((gulf_k-1)*(lon-pigulfphase-38.*pi/180.))  
		 elseif (latorig .le. 15. .AND. latorig .ge. -35.) then
		    flux(i,j) = flux(i,j) &
		      &- trop_atlantic_amp*exp(-(lon*180./pi - 0.)**2./(2*64))*(exp(-(latorig+5.)**2./(2*25)))  &
 		      &- 12.6*exp(-(lon*180./pi + 15.)**2./(2*256))*(exp(-(latorig+16.)**2./(2*64)))      
		 endif
	      endif
	endif


	


	if(  warmpool_localization_choice == 3 .and. trop_atlantic_amp  .gt. 0.001 .AND.  lon .ge. 250.*pi/180.  .AND. lon .lt. 344.*pi/180. .AND. latorig .le. 40. .AND. latorig .ge. -35.  .AND. gulf_amp .gt. 0.001) then !add Caribean and South America at expense of tropical North/South Atlantic  modified by cig, may 14 2019
		    flux(i,j) = flux(i,j) & 
                     &- qflux_amp*.92*exp(-(lon*180./pi-290)**2./(2*400))*(exp(-(latorig+20)**2./(2*49)))   &
		     &- 16.8*exp(-(lon*180./pi-325)**2./(2*484))*(exp(-(latorig-19.5)**2./(2*64)))   &
                     &+ qflux_amp*1.2*exp(-(lon*180./pi-270)**2./(2*49))*exp(-(latorig-22)**2./(2*25)) &
		     &+ qflux_amp*1.58*exp(-(lon*180./pi-283)**2./(2*25))*(exp(-(latorig+0)**2./(2*36)))   &
		     &+ qflux_amp*1.06415*exp(-(lon*180./pi-304)**2./(2*36))*(exp(-(latorig+2)**2./(2*49)))  &
     		     &+ qflux_amp*.85*exp(-(lon*180./pi-284)**2./(2*25))*(exp(-(latorig+10)**2./(2*36)))   &                    
		     &+ qflux_amp*.63*exp(-(lon*180./pi-317)**2./(2*25))*(exp(-(latorig+6)**2./(2*16)))  &
		     &+ (42.54)*exp(-(lon*180./pi-325)**2./(2*121))*(exp(-(latorig-4.2)**2./(2*4)))    
	endif  

        if(  warmpool_localization_choice == 3 .and. africaamp .gt. 0.001 .AND. latorig .le. 35. .AND. latorig .ge. -60. .AND. lon .ge. 2.*pi/180. .AND. lon .lt. 100.*pi/180.) then !add Agulhas current add expense of Africa  modified by cig, may 13 2019
		 flux(i,j) = flux(i,j) & 
                     &- 30.*exp(-(lon*180./pi-28)**2./(2*100))*(exp(-(latorig-18.)**2./(2*50))+exp(-(latorig+18)**2./(2*60))) &
		     &- (38.5+Africaextra*.7709)*exp(-(lon*180./pi-11)**2./(2*4))*exp(-(latorig+15)**2./(2*100)) &    
                     &+ (83.+Africaextra)*exp(-(lon*180./pi-50)**2./(2*625))*exp(-(latorig+40)**2./(2*16)) &
		     &- (64.22+Africaextra*1.3)*exp(-(lon*180./pi-50)**2./(2*400))*exp(-(latorig+48)**2./(2*16)) &
		     &+  (38.0+Africaextra/3)*exp(-(lon*180./pi-2./3.*latorig-57.)**2./(2*16))*exp(-(lon*180./pi+latorig-10.)**2./(2*225.)) &
		     &+ 20.*exp(-(lon*180./pi-14)**2./(2*30))*(exp(-(latorig-0.)**2./(2*50))) &
		     &+ 11.*exp(-(lon*180./pi-36)**2./(2*30))*(exp(-(latorig-0.)**2./(2*50))) 
		    
	endif

      if(  warmpool_localization_choice == 3  .AND. latorig .ge. -61. .AND. latorig .le. -30. .AND. lon .ge. 290.*pi/180. ) then !add dipole in south atlantic  modified by cig, august 4 2019
		 flux(i,j) = flux(i,j) &  
                     &+ (37.4+Africaextra*0*.9349/2)*exp(-(lon*180./pi-323)**2./(2*121))*exp(-(latorig+36)**2./(2*16)) &
		     &- (40.+0*Africaextra/2)*exp(-(lon*180./pi-311)**2./(2*121))*exp(-(latorig+45)**2./(2*16)) 
		
	endif


	if(  warmpool_localization_choice == 3 .and. trop_atlantic_amp  .gt. 0.001  .AND. gulf_amp .gt. 0.001) then !add Norwegian Sea at expense of Africa modified by cig, may 28 2019
		if (lon .ge. 310.*pi/180. .AND. latorig .ge. 10. .AND. latorig .le. 35.) then  
		    flux(i,j) = flux(i,j) & 
                       &- qflux_amp*14.5/26.*exp(-(lon*180./pi-357.)**2./(2*400))*(exp(-(latorig-20)**2./(2*49)))   
		endif
		if (lon .le. 30.*pi/180. .AND. latorig .ge. 10. .AND. latorig .le. 35.) then  
		    flux(i,j) = flux(i,j) & 
                     &- qflux_amp*14.5/26.*exp(-(lon*180./pi+3)**2./(2*400))*(exp(-(latorig-20)**2./(2*49)))  
		endif
		if ((lon .le. 75.*pi/180. .OR. lon .ge. (345.)*pi/180.) .AND. latorig .ge. 71. .AND. latorig .le. 83.) then  
		   flux(i,j) = flux(i,j) &	                    
		     & + qflux_amp*68./26.*(1.-((latorig-76.)/6.5)**4.)*cos(2*(lon-30.*pi/180.)) 
		endif 
	endif  

	if(  warmpool_localization_choice == 3 .and. trop_atlantic_amp  .gt. 0.001 ) then !wave1 over Arctic and Hudson Bay modified by cig, may 30 2019
		if (latorig .ge. 69. .AND. latorig .le. 83.) then  
		    flux(i,j) = flux(i,j) & 
                       & + 25.*(1.-((latorig-76.)/7.)**4.)*cos(1.*(lon-10.*pi/180.))
		endif
		if (latorig .ge. 60. .AND. latorig .le. 76.) then  
		    if ((lon .le. 17.*pi/180. .OR. lon .ge. (347.)*pi/180.) ) then  
		   	flux(i,j) = flux(i,j) &	                    
		     	   & + 68.2*(1.-((latorig-68.)/8.)**4.)*cos(6*(lon-2.*pi/180.)) 
		    endif 
		endif
		if (lon .ge. 260.*pi/180.  .AND. lon .le. 310.*pi/180. .AND. latorig .ge. 55. .AND. latorig .le. 85.) then  
		    flux(i,j) = flux(i,j) & 
                      &-38*exp(-(lon*180./pi-2.*latorig-152.)**2./(2*100.))*exp(-(lon*180./pi+latorig-342.)**2./(2*400.)) &
                      &-100*exp(-(lon*180./pi-275.)**2./(2*25.))*exp(-(latorig-58.)**2./(2*16.)) 
		endif
		if (lon .ge. 275.*pi/180.  .AND. lon .le. 335.*pi/180. .AND. latorig .ge. 10. .AND. latorig .le. 52.) then  
                   flux(i,j) = flux(i,j) & 
                       &+ 10.8*exp(-(lon*180./pi-2.*latorig-220.)**2./(2*100.))*exp(-(lon*180./pi+latorig-335.)**2./(2*625.))
		endif		
	endif  

  
      
       enddo
  enddo
  end subroutine warmpool

!########################################################

end module qflux_mod

