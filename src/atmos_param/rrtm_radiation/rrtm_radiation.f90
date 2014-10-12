
      module rrtm_vars
!
!   Martin Jucker, 2014.
!
!   Contains all variables needed to 
!   run the RRTM code, version for GCMs (hence the 'G'),
!   i.e. all variables needed for radiation that are not
!   within the spectral_physics module otherwise
! 
!   external modules
        use parkind, only         : im => kind_im, rb => kind_rb
        use interpolator_mod, only: interpolate_type
!
!  rrtm_radiation variables
!
        implicit none

        logical                                    :: rrtm_init=.false.    ! has radiation been initialized?
        type(interpolate_type),save                :: o3_interp            ! use external file for ozone
        integer(kind=im)                           :: ncols_rrt,nlay_rrt   ! RRTM field sizes
                                                                           ! ncols_rrt = (size(lon)/lonstep*
                                                                           !             size(lat)
                                                                           ! nlay_rrt  = size(pfull)
        ! gas volume mixing ratios, dimensions (ncols_rrt,nlay=nlevels)
        ! vmr = mass mixing ratio [mmr,kg/kg] scaled with molecular weight [g/mol]
        ! any modification of the units must be accounted for in rrtm_xw_rad.nomica.f90, where x=(s,l)
        real(kind=rb),allocatable,dimension(:,:)   :: h2o                  ! specific humidity [kg/kg]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: o3                   ! ozone [mmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: co2                  ! CO2 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: zeros                ! place holder for any species set
                                                                           ! to zero
        ! the following species are only set if use_secondary_gases=.true.
        real(kind=rb),allocatable,dimension(:,:)   :: ch4                  ! CH4 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: n2o                  ! N2O [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: o2                   ! O2  [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc11                ! CFC11 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc12                ! CFC12 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc22                ! CFC22 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cc14                 ! CC14 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: emis                 ! surface LW emissivity per band
                                                                           ! dimension (ncols_rrt x nbndlw)
                                                                           ! =1 for black body
        ! clouds stuff = set to zero
        !mjdb CHECK IF I CAN ACTUALLY REPLACE THOSE WITH ZEROS
        real(kind=rb),allocatable,dimension(:,:)   :: cldfr, cicewp, cliqwp, &
             reice, reliq
        ! cloud & aerosol optical depths, cloud and aerosol specific parameters
        !mjdb AGAIN, CHECK IF CAN REPLACE WITH ZEROS
        real(kind=rb),allocatable,dimension(:,:,:) :: taucld,tauaer, &
             ssacld,asmcld,fsfcld,ssaaer,asmaer,ecaer
        ! heating rates and fluxes, zenith angle when in-between radiation time steps
        real(kind=rb),allocatable,dimension(:,:)   :: sw_flux,lw_flux,zencos! surface fluxes, cos(zenith angle) 
                                                                            ! dimension (lon x lat)
        real(kind=rb),allocatable,dimension(:,:,:) :: tdt_rad               ! heating rate [K/s]
                                                                            ! dimension (lon x lat x pfull)
        real(kind=rb),allocatable,dimension(:,:,:) :: tdt_sw_rad,tdt_lw_rad ! SW, LW radiation heating rates,
                                                                            ! diagnostics only [K/s]
                                                                            ! dimension (lon x lat x pfull)
        real(kind=rb),allocatable,dimension(:,:,:) :: t_half                ! temperature at half levels [K]
                                                                            ! dimension (lon x lat x phalf)
        real(kind=rb),allocatable,dimension(:,:)   :: rrtm_precip           ! total time of precipitation
                                                                            ! between radiation steps to
                                                                            ! determine precip_albedo
                                                                            ! dimension (lon x lat)
        integer(kind=im)                           :: num_precip            ! number of times precipitation
                                                                            ! has been summed in rrtm_precip
        integer(kind=im)                           :: dt_last               ! time of last radiation calculation
                                                                            ! used for alarm
!---------------------------------------------------------------------------------------------------------------
! some constants
        real(kind=rb)      :: daypersec=1./86400,deg2rad
! no clouds in the radiative scheme
        integer(kind=im) :: icld=0,idrv=0, &
             inflglw=0,iceflglw=0,liqflglw=0, &
             iaer=0
!---------------------------------------------------------------------------------------------------------------
!                                namelist values
!---------------------------------------------------------------------------------------------------------------
        logical            :: include_secondary_gases=.false. ! non-zero values for above listed secondary gases?
        logical            :: do_read_ozone=.false.           ! read ozone from an external file?
        character(len=256) :: ozone_file='ozone_1990'         !  file name of ozone file to read
        real(kind=rb)      :: h2o_lower_limit = 2.e-7         ! never use smaller than this in radiative scheme
        real(kind=rb)      :: temp_lower_limit = 100.         ! never go below this in radiative scheme
        real(kind=rb)      :: temp_upper_limit = 370.         ! never go above this in radiative scheme
        real(kind=rb)      :: co2ppmv=300.                    ! CO2 ppmv concentration
        real(kind=rb)      :: solr_cnst= 1368.22              ! solar constant [W/m2]
        logical            :: use_dyofyr=.true.               ! use day of year for solrad calculation?
        real(kind=rb)      :: solrad=1.0                      ! distance Earth-Sun [AU] if use_dyofyr=.false.
        real(kind=rb)      :: solday=0.                       ! if >0, do perpetual run corresponding to 
                                                              !  day of the year = solday
        logical            :: store_intermediate_rad =.true.  ! Keep rad constant over entire dt_rad?
                                                              ! Else only heat radiatively at every dt_rad
        logical            :: do_rad_time_avg =.true.         ! Average radiation over dt_rad?
        integer(kind=im)   :: dt_rad=900                      ! Radiation time step
        integer(kind=im)   :: lonstep=1                       ! Subsample fields along longitude
                                                              !  for faster radiation calculation  
        logical            :: do_precip_albedo=.false.        ! Modify albedo depending on large scale
                                                              !  precipitation (crude cloud parameterization)
        real(kind=rb)      :: precip_albedo=0.8               ! If so, what's the cloud albedo?
!---------------------------------------------------------------------------------------------------------------
!
!-------------------- diagnostics fields -------------------------------

        integer :: id_tdt_rad,id_tdt_sw,id_tdt_lw,id_coszen,id_flux_sw,id_flux_lw,id_albedo
        character(len=14), parameter :: mod_name = 'rrtm_radiation'
        real :: missing_value = -999.

!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------

        namelist/rrtm_radiation_nml/ include_secondary_gases, do_read_ozone, ozone_file, &
             &h2o_lower_limit,temp_lower_limit,temp_upper_limit,co2ppmv, &
             &solr_cnst, solrad, use_dyofyr, solday, &
             &store_intermediate_rad, do_rad_time_avg, dt_rad, lonstep, &
             &do_precip_albedo,precip_albedo

      end module rrtm_vars
!*****************************************************************************************
!*****************************************************************************************
      module rrtm_radiation
        use parkind, only : im => kind_im, rb => kind_rb
        implicit none
    
      contains

!*****************************************************************************************
        subroutine rrtm_radiation_init(axes,Time,ncols,nlay,lonb,latb)
! Martin Jucker 2014
!
! Initialize diagnostics, allocate variables, set constants
!
! Modules
          use rrtm_vars
          use parrrtm, only:          nbndlw
          use parrrsw, only:          nbndsw
          use diag_manager_mod, only: register_diag_field, send_data
          use interpolator_mod, only: interpolate_type, interpolator_init, &
                                      &CONSTANT, INTERP_WEIGHTED_P
          use fms_mod, only:          open_namelist_file, check_nml_error,  &
                                      &mpp_pe, mpp_root_pe, close_file, &
                                      &write_version_number, stdlog, &
                                      &error_mesg, NOTE, WARNING
          use time_manager_mod, only: time_type
! Local variables
          implicit none
          
          integer, intent(in), dimension(4) :: axes
          type(time_type), intent(in)       :: Time
          integer(kind=im),intent(in)       :: ncols,nlay
          real(kind=rb),dimension(:),intent(in),optional :: lonb,latb

          integer :: i,k,seconds

          integer :: ierr, io, unit


! read namelist and copy to logfile
          unit = open_namelist_file ( )
          ierr=1
          do while (ierr /= 0)
             read  (unit, nml=rrtm_radiation_nml, iostat=io, end=10)
             ierr = check_nml_error (io, 'rrtm_radiation_nml')
          enddo
10        call close_file (unit)
          
          !call write_version_number ( version, tagname )
          if ( mpp_pe() == mpp_root_pe() ) then
             write (stdlog(), nml=rrtm_radiation_nml)
          endif
          call close_file (unit)
!----
!------------ initialize diagnostic fields ---------------

          id_tdt_rad = &
               register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
                 'Temperature tendency due to radiation', &
                 'K/s', missing_value=missing_value               )
          id_tdt_sw = &
               register_diag_field ( mod_name, 'tdt_sw', axes(1:3), Time, &
                 'Temperature tendency due to SW radiation', &
                 'K/s', missing_value=missing_value               )
          id_tdt_lw = &
               register_diag_field ( mod_name, 'tdt_lw', axes(1:3), Time, &
                 'Temperature tendency due to LW radiation', &
                 'K/s', missing_value=missing_value               )
          id_coszen  = &
               register_diag_field ( mod_name, 'coszen', axes(1:2), Time, &
                 'cosine of zenith angle', &
                 'none', missing_value=missing_value               )
          id_flux_sw = &
               register_diag_field ( mod_name, 'flux_sw', axes(1:2), Time, &
                 'Net SW surface flux', &
                 'W/m2', missing_value=missing_value               )
          id_flux_lw = &
               register_diag_field ( mod_name, 'flux_lw', axes(1:2), Time, &
                 'LW surface flux', &
                 'W/m2', missing_value=missing_value               )
          id_albedo  = &
               register_diag_field ( mod_name, 'rrtm_albedo', axes(1:2), Time, &
                 'Interactive albedo', &
                 'none', missing_value=missing_value               )
! 
!------------ make sure namelist choices are consistent -------
! this does not work at the moment, as dt_atmos from coupler_mod induces a circular dependency at compilation
!          if(dt_rad .le. dt_atmos .and. store_intermediate_rad)then
!             call error_mesg ( 'rrtm_gases_init', &
!                  ' dt_rad <= dt_atmos, for conserving memory, I am setting store_intermediate_rad=.false.', &
!                  WARNING)
!             store_intermediate_rad = .false.
!          endif
!          if(dt_rad .gt. dt_atmos .and. .not.store_intermediate_rad)then
!             call error_mesg( 'rrtm_gases_init', &
!                  ' dt_rad > dt_atmos, but store_intermediate_rad=.false. might cause time steps with zero radiative forcing!', &
!                  WARNING)
!          endif
          if(solday .gt. 0)then
             call error_mesg( 'rrtm_gases_init', &
                  ' running perpetual simulation', NOTE)
          endif
!------------ set some constants and parameters -------

          deg2rad = acos(0.)/90.

          dt_last = -dt_rad

          ncols_rrt = ncols/lonstep
          nlay_rrt  = nlay

!------------ allocate arrays to be used later  -------
          allocate(t_half(size(lonb,1)-1,size(latb)-1,nlay+1))

          allocate(h2o(ncols_rrt,nlay_rrt),o3(ncols_rrt,nlay_rrt), &
               co2(ncols_rrt,nlay_rrt))
          if(include_secondary_gases)then
             allocate(ch4(ncols_rrt,nlay_rrt), &
               n2o(ncols_rrt,nlay_rrt),o2(ncols_rrt,nlay_rrt), &
               cfc11(ncols_rrt,nlay_rrt),cfc12(ncols_rrt,nlay_rrt), &
               cfc22(ncols_rrt,nlay_rrt),cc14(ncols_rrt,nlay_rrt), &
               reice(ncols_rrt,nlay_rrt),reliq(ncols_rrt,nlay_rrt))
          else
             allocate(zeros(ncols_rrt,nlay_rrt))
          endif
          allocate(emis(ncols_rrt,nbndlw))
          allocate(cldfr(ncols_rrt,nlay_rrt),cicewp(ncols_rrt,nlay_rrt), &
               cliqwp(ncols_rrt,nlay_rrt))
          allocate(taucld(nbndlw,ncols_rrt,nlay_rrt),tauaer(ncols_rrt, &
               nlay_rrt,nbndlw))
          allocate(ssacld(nbndsw,ncols_rrt,nlay_rrt), &
               asmcld(nbndsw,ncols_rrt,nlay_rrt), &
               fsfcld(nbndsw,ncols_rrt,nlay_rrt), &
               ssaaer(ncols_rrt,nlay_rrt,nbndsw), &
               asmaer(ncols_rrt,nlay_rrt,nbndsw), &
               ecaer(ncols_rrt,nlay_rrt,nbndsw))
          if(store_intermediate_rad .or. id_flux_sw > 0) &
               allocate(sw_flux(size(lonb,1)-1,size(latb,1)-1))
          if(store_intermediate_rad .or. id_flux_lw > 0) &
               allocate(lw_flux(size(lonb,1)-1,size(latb,1)-1))
          if(id_coszen > 0)allocate(zencos (size(lonb,1)-1,size(latb,1)-1))
          if(do_precip_albedo)allocate(rrtm_precip(size(lonb,1)-1,size(latb,1)-1))
          if(store_intermediate_rad .or. id_tdt_rad > 0)&
               allocate(tdt_rad(size(lonb,1)-1,size(latb,1)-1,nlay))
          if(id_tdt_sw > 0)allocate(tdt_sw_rad(size(lonb,1)-1,size(latb,1)-1,nlay)) 
          if(id_tdt_lw > 0)allocate(tdt_lw_rad(size(lonb,1)-1,size(latb,1)-1,nlay)) 

          h2o   = 0. !this will be set by the water vapor tracer
          o3    = 0. !this will be set by an input file if do_read_ozone=.true.
          co2   = co2ppmv*1.e-6
          if(include_secondary_gases)then !change values here
             ch4   = 0.
             n2o   = 0.
             o2    = 0.
             cfc11 = 0.
             cfc12 = 0.
             cfc22 = 0.
             cc14  = 0.
          else
             zeros = 0.
          endif

          emis  = 1. !black body: 1.0

          cldfr  = 0.
          cicewp = 0.
          cliqwp = 0.
          reice = 10.
          reliq = 10.
          
          taucld = 0.
          tauaer = 0.

          ssacld = 0.
          asmcld = 0.
          fsfcld = 0.
          ssaaer = 0.
          asmaer = 0.
          ecaer  = 0.

          if(do_read_ozone)then
             call interpolator_init (o3_interp, trim(ozone_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
          endif

          if(do_precip_albedo)then
             rrtm_precip = 0.
             num_precip  = 0
          endif

          rrtm_init=.true.

        end subroutine rrtm_radiation_init
!*****************************************************************************************
        subroutine interp_temp(z_full,z_half,t)
          use rrtm_vars
          implicit none

          real(kind=rb),dimension(:,:,:),intent(in)  :: z_full,z_half,t

          integer i,j,k,kend
          real dzk,dzk1,dzk2

! note: z_full(kend) = z_half(kend), so there's something fishy
! also, for some reason, z_half(k=1)=0. so we need to deal with k=1 separately
          kend=size(z_full,3)
          do k=2,kend
             do j=1,size(t,2)
                do i=1,size(t,1)
                   dzk2 = 1./( z_full(i,j,k-1)   - z_full(i,j,k) )
                   dzk  = ( z_half(i,j,k  )   - z_full(i,j,k) )*dzk2
                   dzk1 = ( z_full(i,j,k-1)   - z_half(i,j,k) )*dzk2 
                   t_half(i,j,k) = t(i,j,k)*dzk1 + t(i,j,k-1)*dzk
                enddo
             enddo
          enddo
! top of the atmosphere: need to extrapolate. z_half(1)=0, so need to use values on full grid
          do j=1,size(t,2)
             do i=1,size(t,1)
                !standard linear extrapolation
                !top: use full points, and distance is 1.5 from k=2
                t_half(i,j,1) = 0.5*(3*t(i,j,1)-t(i,j,2))
                !bottom: z=0 => distance is -z_full(kend-1)/(z_full(kend)-z_full(kend-1))
                t_half(i,j,kend+1) = t(i,j,kend-1) &
                     + (z_half(i,j,kend+1) - z_full(i,j,kend-1))&
                     * (t     (i,j,kend  ) - t     (i,j,kend-1))&
                     / (z_full(i,j,kend  ) - z_full(i,j,kend-1))
             enddo
          enddo
          

        end subroutine interp_temp
!*****************************************************************************************
!*****************************************************************************************
        subroutine run_rrtmg(is,js,Time,lat,lon,p_full,p_half,albedo,q,t,t_surf_rad,tdt,coszen,flux_sw,flux_lw)
!
! Martin Jucker 2014
!
! Driver for RRTMG radiation scheme.
! Prepares all inputs, calls SW and LW radiation schemes, 
!  transforms outputs back into FMS form
!
! Modules
          use fms_mod, only:         error_mesg, FATAL
          use mpp_mod, only:         mpp_pe,mpp_root_pe
          use rrtmg_lw_rad, only:    rrtmg_lw
          use rrtmg_sw_rad, only:    rrtmg_sw
          use rrtm_astro, only:      compute_zenith
          use rrtm_vars
          use time_manager_mod,only: time_type,get_time,set_time
          use interpolator_mod,only: interpolator
!---------------------------------------------------------------------------------------------------------------
! In/Out variables
          implicit none

          integer, intent(in)                               :: is, js          ! index range for each CPU
          type(time_type),intent(in)                        :: Time            ! global time in calendar
          real(kind=rb),dimension(:,:,:),intent(in)         :: p_full,p_half   ! pressure, full and half levels
                                                                               ! dimension (lat x lon x p*)
          real(kind=rb),dimension(:,:,:),intent(in)         :: q               ! water vapor mixing ratio [g/g]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:,:),intent(in)         :: t               ! temperature [K]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:),intent(in)           :: lat,lon         ! latitude, longitude
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(in)           :: albedo          ! surface albedo
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(in)           :: t_surf_rad      ! surface temperature [K]
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:,:),intent(inout)      :: tdt             ! heating rate [K/s]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:),intent(out)          :: coszen          ! cosine of zenith angle
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(out),optional :: flux_sw,flux_lw ! surface fluxes [W/m2]
                                                                               ! dimension (lat x lon)
                                                                               ! need to have both or none!
!---------------------------------------------------------------------------------------------------------------
! Local variables
          integer k,j,i,ij,j1,i1,ij1,kend,dyofyr,seconds
          integer si,sj,sk,locmin(3)
          real(kind=rb),dimension(size(q,1),size(q,2),size(q,3)) :: o3f
          real(kind=rb),dimension(ncols_rrt,nlay_rrt) :: pfull,tfull,fracday&
               , hr,hrc, swhr, swhrc
          real(kind=rb),dimension(size(tdt,1),size(tdt,2),size(tdt,3)) :: tdt_rrtm
          real(kind=rb),dimension(ncols_rrt,nlay_rrt+1) :: uflx, dflx, uflxc, dflxc&
               ,swuflx, swdflx, swuflxc, swdflxc
          real(kind=rb),dimension(size(q,1)/lonstep,size(q,2),size(q,3)  ) :: swijk,lwijk
          real(kind=rb),dimension(size(q,1)/lonstep,size(q,2)) :: swflxijk,lwflxijk
          real(kind=rb),dimension(ncols_rrt,nlay_rrt+1):: phalf,thalf
          real(kind=rb),dimension(ncols_rrt)   :: tsrf,cosz_rr,albedo_rr
          real(kind=rb) :: dlon,dlat,dj,di 
          type(time_type) :: Time_loc
          real(kind=rb),dimension(size(q,1),size(q,2)) :: albedo_loc
!---------------------------------------------------------------------------------------------------------------

          if(.not. rrtm_init)&
               call error_mesg('run_rrtm','module not initialized', FATAL)

!check if we really want to recompute radiation (alarm)
          call get_time(Time,seconds)
          if(seconds < dt_last) dt_last=dt_last-86400 !it's a new day
          if(seconds - dt_last .ge. dt_rad) then
             dt_last = seconds
          else
             if(store_intermediate_rad)then
                tdt_rrtm = tdt_rad
                flux_sw = sw_flux
                flux_lw = lw_flux
             else
                tdt_rrtm = 0.
                flux_sw  = 0.
                flux_lw  = 0.
             endif
             tdt = tdt + tdt_rrtm
             call write_diag_rrtm(Time,is,js)
             return !not time yet
          endif
!
! we know now that we want to run radiation
!
!make sure we run perpetual when solday > 0)
          if(solday > 0.)then
             Time_loc = set_time(seconds,floor(solday))
          else
             Time_loc = Time
          endif
!compute zenith angle
          if(do_rad_time_avg) then
             call compute_zenith(Time_loc,dt_rad,lat,lon,coszen,dyofyr)
          else
             call compute_zenith(Time_loc,0     ,lat,lon,coszen,dyofyr)
          end if

          si=size(tdt,1)
          sj=size(tdt,2)
          sk=size(tdt,3)

          if(.not. use_dyofyr) dyofyr=0 !use solrad instead of day of year

          !get ozone 
          if(do_read_ozone)then
             call interpolator( o3_interp, Time_loc, p_half, o3f, trim(ozone_file))
          endif

          !interactive albedo
          if(do_precip_albedo .and. num_precip>0)then
             albedo_loc = albedo + (precip_albedo - albedo)*rrtm_precip/num_precip
             rrtm_precip = 0.
             num_precip  = 0
          else
             albedo_loc = albedo
          endif
!---------------------------------------------------------------------------------------------------------------
          !RRTM's first pressure level is at the surface - need to inverse order
          !also, RRTM's pressures are in hPa
          !reshape arrays
          pfull = reshape(p_full(1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))*0.01
          phalf = reshape(p_half(1:si:lonstep,:,sk+1:1:-1),(/ si*sj/lonstep,sk+1 /))*0.01
          tfull = reshape(t     (1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))
          thalf = reshape(t_half(1:si:lonstep,:,sk+1:1:-1),(/ si*sj/lonstep,sk+1 /))
          h2o   = reshape(q     (1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))
          if(do_read_ozone)o3 = reshape(o3f(1:si:lonstep,:,sk :1:-1),(/ si*sj/lonstep,sk  /))
          
          cosz_rr   = reshape(coszen    (1:si:lonstep,:),(/ si*sj/lonstep /))
          albedo_rr = reshape(albedo_loc(1:si:lonstep,:),(/ si*sj/lonstep /))
          tsrf      = reshape(t_surf_rad(1:si:lonstep,:),(/ si*sj/lonstep /))
          
!---------------------------------------------------------------------------------------------------------------
! now actually run RRTM
!
          swhr = 0.
          swdflx = 0.
          swuflx = 0.
          ! make sure we don't go beyond 'dangerous' values for radiation. Same as Merlis spectral_am2rad
          h2o   = max(h2o  , h2o_lower_limit)
          tfull = max(tfull, temp_lower_limit)
          tfull = min(tfull, temp_upper_limit)
          thalf = max(thalf, temp_lower_limit)
          thalf = min(thalf, temp_upper_limit)

          
          if(include_secondary_gases)then
             call rrtmg_sw &
                  (ncols_rrt, nlay_rrt , icld     , iaer     , &
                  pfull     , phalf    , tfull    , thalf    , tsrf , &
                  h2o       , o3       , co2      , ch4      , n2o  , o2, &
                  albedo_rr , albedo_rr, albedo_rr, albedo_rr, &
                  cosz_rr   , solrad   , dyofyr   , solr_cnst, &
                  inflglw   , iceflglw , liqflglw , cldfr    , &
                  taucld    , ssacld   , asmcld   , fsfcld   , &
                  cicewp    , cliqwp   , reice    , reliq    , &
                  tauaer    , ssaaer   , asmaer   , ecaer    , &
                  swuflx    , swdflx   , swhr     , swuflxc  , swdflxc, swhrc)
          else
             call rrtmg_sw &
                  (ncols_rrt, nlay_rrt , icld     , iaer     , &
                  pfull     , phalf    , tfull    , thalf    , tsrf , &
                  h2o       , o3       , co2      , zeros    , zeros, zeros, &
                  albedo_rr , albedo_rr, albedo_rr, albedo_rr, &
                  cosz_rr   , solrad   , dyofyr   , solr_cnst, &
                  inflglw   , iceflglw , liqflglw , cldfr    , &
                  taucld    , ssacld   , asmcld   , fsfcld   , &
                  cicewp    , cliqwp   , reice    , reliq    , &
                  tauaer    , ssaaer   , asmaer   , ecaer    , &
                  swuflx    , swdflx   , swhr     , swuflxc  , swdflxc, swhrc)
          endif
          
          ! make sure we don't have SW radiation at night
          ! there is some optimization possible here: only feed grid points to rrtm_sw where cosz_rr>0
          do i=1,size(swhr,2)
             where( cosz_rr == 0. )
                swuflx(:,i) = 0.
                swdflx(:,i) = 0.
                swhr  (:,i) = 0.
             endwhere
          enddo
             

          swijk   = reshape(swhr(:,sk:1:-1),(/ si/lonstep,sj,sk /))*daypersec

          hr = 0.
          dflx = 0.
          uflx = 0.
          if(include_secondary_gases)then
             call rrtmg_lw &
                  (ncols_rrt, nlay_rrt, icld    , idrv , &
                  pfull     , phalf   , tfull   , thalf, tsrf , &
                  h2o       , o3      , co2     , ch4  , n2o  , o2, &
                  cfc11     , cfc12   , cfc22   , cc14 , emis , &
                  inflglw   , iceflglw, liqflglw, cldfr,  &
                  taucld    , cicewp  , cliqwp  , reice, reliq, &
                  tauaer    , &
                  uflx      , dflx    , hr      , uflxc, dflxc, hrc)
          else
             call rrtmg_lw &
                  (ncols_rrt, nlay_rrt, icld    , idrv , &
                  pfull     , phalf   , tfull   , thalf, tsrf , &
                  h2o       , o3      , co2     , zeros, zeros, zeros, &
                  zeros     , zeros   , zeros   , zeros, emis , &
                  inflglw   , iceflglw, liqflglw, cldfr,  &
                  taucld    , cicewp  , cliqwp  , reice, reliq, &
                  tauaer    , &
                  uflx      , dflx    , hr      , uflxc, dflxc, hrc)
          endif

          lwijk   = reshape(hr(:,sk:1:-1),(/ si/lonstep,sj,sk /))*daypersec

!---------------------------------------------------------------------------------------------------------------
! interpolate back onto GCM grid (latitude is kept the same due to parallelisation)
          dlon=1./lonstep
          do i=1,size(swijk,1)
             i1 = i+1
             ! close toroidally
             if(i1 > size(swijk,1)) i1=1
             do ij=1,lonstep
                di = (ij-1)*dlon
                ij1 = (i-1)*lonstep + ij
                tdt_rrtm(ij1,:,:) =  &
                     +       di *(swijk(i1,:,:) + lwijk(i1,:,:)) &
                     +   (1.-di)*(swijk(i ,:,:) + lwijk(i ,:,:))
                if(id_tdt_sw>0)tdt_sw_rad(ij1,:,:)=di*swijk(i1,:,:)+(1.-di)*swijk(i,:,:)
                if(id_tdt_lw>0)tdt_lw_rad(ij1,:,:)=di*lwijk(i1,:,:)+(1.-di)*lwijk(i,:,:)
             enddo
          enddo
          tdt = tdt + tdt_rrtm
          ! store radiation between radiation time steps
          if(store_intermediate_rad .or. id_tdt_rad > 0) tdt_rad = tdt_rrtm
         
          ! get the surface fluxes
          if(present(flux_sw).and.present(flux_lw))then
             !only surface fluxes are needed
             swflxijk = reshape(swdflx(:,1)-swuflx(:,1),(/ si/lonstep,sj /)) ! net down SW flux
             lwflxijk = reshape(  dflx(:,1)            ,(/ si/lonstep,sj /)) ! down LW flux
             dlon=1./lonstep
             do i=1,size(swijk,1)
                i1 = i+1
                ! close toroidally
                if(i1 > size(swijk,1)) i1=1
                do ij=1,lonstep
                   di = (ij-1)*dlon
                   ij1 = (i-1)*lonstep + ij
                   flux_sw(ij1,:) = di*swflxijk(i1,:) + (1.-di)*swflxijk(i ,:)
                   flux_lw(ij1,:) = di*lwflxijk(i1,:) + (1.-di)*lwflxijk(i ,:)
                enddo
             enddo
             ! store between radiation steps
             if(store_intermediate_rad)then
                sw_flux = flux_sw
                lw_flux = flux_lw
             else
                if(id_flux_sw > 0)sw_flux = flux_sw
                if(id_flux_lw > 0)lw_flux = flux_lw
             endif
             if(id_coszen  > 0)zencos  = coszen
          endif
          
          if(do_precip_albedo)then
             call write_diag_rrtm(Time,is,js,albedo_loc)
          else
             call write_diag_rrtm(Time,is,js)
          endif
        end subroutine run_rrtmg

!*****************************************************************************************
!*****************************************************************************************
        subroutine write_diag_rrtm(Time,is,js,albedo_loc)
!
! Martin Jucker
! 
! write out diagnostics fields
!
! Modules
          use rrtm_vars,only:         sw_flux,lw_flux,zencos,tdt_rad,tdt_sw_rad,tdt_lw_rad,&
                                      &id_tdt_rad,id_tdt_sw,id_tdt_lw,id_coszen,&
                                      &id_flux_sw,id_flux_lw,id_albedo
          use diag_manager_mod, only: register_diag_field, send_data
          use time_manager_mod,only:  time_type
! Input variables
          implicit none
          type(time_type),intent(in)              :: Time
          integer, intent(in)                     :: is, js
          real,dimension(:,:),intent(in),optional :: albedo_loc
! Local variables
          logical :: used

!------- temperature tendency due to radiation ------------
          if ( id_tdt_rad > 0 ) then
             used = send_data ( id_tdt_rad, tdt_rad, Time, is, js, 1 )
          endif
!------- temperature tendency due to SW radiation ---------
          if ( id_tdt_sw > 0 ) then
             used = send_data ( id_tdt_sw, tdt_sw_rad, Time, is, js, 1 )
          endif
!------- temperature tendency due to LW radiation ---------
          if ( id_tdt_lw > 0 ) then
             used = send_data ( id_tdt_lw, tdt_lw_rad, Time, is, js, 1 )
          endif
!------- cosine of zenith angle                ------------
          if ( id_coszen > 0 ) then
             used = send_data ( id_coszen, zencos, Time, is, js )
          endif
!------- Net SW surface flux                   ------------
          if ( id_flux_sw > 0 ) then
             used = send_data ( id_flux_sw, sw_flux, Time, is, js )
          endif
!------- Net LW surface flux                   ------------
          if ( id_flux_lw > 0 ) then
             used = send_data ( id_flux_lw, lw_flux, Time, is, js )
          endif
!------- Interactive albedo                    ------------
          if ( present(albedo_loc)) then
             used = send_data ( id_albedo, albedo_loc, Time, is, js )
          endif
        end subroutine write_diag_rrtm
!*****************************************************************************************

        subroutine rrtm_radiation_end
          use rrtm_vars, only: do_read_ozone,o3_interp
          use interpolator_mod, only: interpolator_end
          implicit none

          if(do_read_ozone)call interpolator_end(o3_interp)
        end subroutine rrtm_radiation_end


!*****************************************************************************************
      end module rrtm_radiation
       


        
