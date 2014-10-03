      module rrtm_vars
        use parkind, only : im => kind_im, rb => kind_rb
        use interpolator_mod, only: interpolate_type
        implicit none

        logical :: rrtm_init=.false.
        type(interpolate_type),save :: o3_interp
        integer(kind=im) :: ncols_rrt,nlay_rrt
        ! gas volume mixing ratios, dimensions (ncols=nlon*nlat,nlay=nlevels)
        ! vmr = mass mixing ratio [g/g] scaled with molecular weight [g/mol]
        real(kind=rb),allocatable,dimension(:,:) :: h2o, o3, co2, zeros, &
             &ch4, n2o, o2, cfc11, cfc12, cfc22, cc14
        ! clouds
        real(kind=rb),allocatable,dimension(:,:) :: cldfr, cicewp, cliqwp, &
             reice, reliq
        ! half-level temperature, emssivity
        real(kind=rb),allocatable,dimension(:,:)   :: emis
        ! cloud & aerosol optical depths, cloud and aerosol specific parameters
        real(kind=rb),allocatable,dimension(:,:,:) :: taucld,tauaer, &
             ssacld,asmcld,fsfcld,ssaaer,asmaer,ecaer
        ! heating rates and fluxes when not re-computing
        real(kind=rb),allocatable,dimension(:,:)   :: sw_flux,lw_flux,zencos
        real(kind=rb),allocatable,dimension(:,:,:) :: tdt_rad,tdt_sw_rad,tdt_lw_rad,t_half
        real(kind=rb)      :: daypersec=1./86400,deg2rad
        integer(kind=im)   :: dt_last
! namelist values
        logical            :: include_secondary_gases=.false.
        logical            :: do_read_ozone=.false.
        character(len=256) :: ozone_file='ozone_1990'
        real(kind=rb)      :: h2o_lower_limit = 2.e-7  !never use small than this in radiation
        real(kind=rb)      :: temp_lower_limit = 100. !never go below this in radiation
        real(kind=rb)      :: temp_upper_limit = 370. !never go above this in radiation
        real(kind=rb)      :: co2ppmv=300.
        real(kind=rb)      :: solr_cnst= 1368.22 ! solar constant [1368.22 W/m2]
        real(kind=rb)      :: solrad=1.0 ! distance Earth-Sun [AU]
        logical            :: use_dyofyr=.true. ! if true, use day of year for solrad calculation
        real(kind=im)      :: solday=0.   ! if >0, do perpetual run corresponding to day of the year = solday
        logical            :: store_intermediate_rad =.true. !if true, keep rad constant over entire dt_rad, else only head radiatively at every dt_rad
        logical            :: do_rad_time_avg =.true. !if true, average radiation over dt_rad
        integer(kind=im)   :: dt_rad=900,lonstep=1
!
        integer(kind=im) :: icld=0,idrv=0, &
             inflglw=0,iceflglw=0,liqflglw=0, &
             iaer=0
!-------------------- diagnostics fields -------------------------------

        integer :: id_tdt_rad,id_tdt_sw,id_tdt_lw,id_coszen,id_flux_sw,id_flux_lw
        character(len=14), parameter :: mod_name = 'rrtm_radiation'
        real :: missing_value = -999.

        namelist/rrtm_radiation_nml/ include_secondary_gases, do_read_ozone, ozone_file, &
             &h2o_lower_limit,temp_lower_limit,temp_upper_limit,co2ppmv, &
             &solr_cnst, solrad, use_dyofyr, solday, &
             &store_intermediate_rad, do_rad_time_avg, dt_rad, lonstep, &
             &icld, idrv, inflglw, iceflglw, liqflglw, iaer

      end module rrtm_vars
!*****************************************************************************************
!*****************************************************************************************
      module rrtm_gases
        use parkind, only : im => kind_im, rb => kind_rb
        implicit none
    
      contains

!*****************************************************************************************
        subroutine rrtm_gases_init(axes,Time,ncols,nlay,lonb,latb)
          use rrtm_vars
          use parrrtm, only : nbndlw
          use parrrsw, only:  nbndsw
          use diag_manager_mod, only: register_diag_field, send_data
          use interpolator_mod, only: interpolate_type, interpolator_init, &
               CONSTANT, INTERP_WEIGHTED_P
          use fms_mod, only: open_namelist_file, check_nml_error,  &
                                    mpp_pe, mpp_root_pe, close_file, &
                                    write_version_number, stdlog, &
                                    error_mesg, NOTE, WARNING
          use time_manager_mod,only: time_type,get_time
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
! 
!------------ make sure namelist choices are consistent -------
          call get_time(Time,seconds)
          if(dt_rad .le. seconds .and. store_intermediate_rad)then
             store_intermediate_rad =.false.
             call error_mesg ( 'rrtm_gases_init', &
                  ' dt_rad <= dt_atmos, setting store_intermediate_rad=.false.', &
                  NOTE)
          endif
          if(dt_rad .gt. seconds .and. .not.store_intermediate_rad)then
             call error_mesg( 'rrtm_gases_init', &
                  ' dt_rad > dt_atmos, but store_intermediate_rad=.false. might cause time steps with zero radiative forcing!', &
                  WARNING)
          endif
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
          if(store_intermediate_rad .or. id_tdt_rad > 0)&
               allocate(tdt_rad(size(lonb,1)-1,size(latb,1)-1,nlay))
          if(id_tdt_sw > 0)allocate(tdt_sw_rad(size(lonb,1)-1,size(latb,1)-1,nlay)) 
          if(id_tdt_lw > 0)allocate(tdt_lw_rad(size(lonb,1)-1,size(latb,1)-1,nlay)) 

          h2o   = 0. !this will be set by the water vapor tracer
          o3    = 0. !this will be set by an input file
          co2   = co2ppmv*1.e-6
          if(include_secondary_gases)then
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

          rrtm_init=.true.

        end subroutine rrtm_gases_init
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

        subroutine run_rrtmg(is,js,Time,lat,lon,p_full,p_half,albedo,q,t,t_surf_rad,tdt,coszen,flux_sw,flux_lw)
          use fms_mod, only:  error_mesg, FATAL
          use mpp_mod, only: mpp_pe,mpp_root_pe
          use rrtmg_lw_rad, only: rrtmg_lw
          use rrtmg_sw_rad, only: rrtmg_sw
          use rrtm_astro, only: compute_zenith
          use rrtm_vars
          use time_manager_mod,only: time_type,get_time
          use interpolator_mod,only: interpolator
          implicit none

          integer, intent(in)                 :: is, js
          type(time_type),intent(in) :: Time
          real(kind=rb),dimension(:,:,:),intent(in) :: p_full,p_half,q,t
          real(kind=rb),dimension(:,:),intent(in)   :: lat,lon,&
               albedo,t_surf_rad
          real(kind=rb),dimension(:,:,:),intent(inout) :: tdt
          real(kind=rb),dimension(:,:),intent(out)     :: coszen
          real(kind=rb),dimension(:,:),intent(out),optional    :: flux_sw,flux_lw !need to have both or none!
          
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
          logical   :: do_rad

          if(.not. rrtm_init)&
               call error_mesg('run_rrtm','module not initialized', FATAL)

!check if we really want to recompute radiation
          call get_time(Time,seconds)
          if(seconds < dt_last) dt_last=dt_last-86400 !it's a new day
          do_rad=.false.
          if(seconds - dt_last .ge. dt_rad) do_rad=.true.
!compute zenith angle at every time step anyway
          if(do_rad_time_avg .and. do_rad) then
             call compute_zenith(Time,dt_rad,lat,lon,coszen,dyofyr)
          elseif(do_rad) then
             call compute_zenith(Time,0     ,lat,lon,coszen,dyofyr)
          end if

          if( do_rad )then
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

          si=size(tdt,1)
          sj=size(tdt,2)
          sk=size(tdt,3)

          if(.not. use_dyofyr) dyofyr=0 !use solrad instead of day of year

          !get ozone 
          if(do_read_ozone)then
             call interpolator( o3_interp, Time, p_half, o3f, trim(ozone_file))
          endif

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
          albedo_rr = reshape(albedo    (1:si:lonstep,:),(/ si*sj/lonstep /))
          tsrf      = reshape(t_surf_rad(1:si:lonstep,:),(/ si*sj/lonstep /))
          
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
!!$                  zeros       , o3       , co2      , zeros    , zeros, zeros, &
                  albedo_rr , albedo_rr, albedo_rr, albedo_rr, &
                  cosz_rr   , solrad   , dyofyr   , solr_cnst, &
                  inflglw   , iceflglw , liqflglw , cldfr    , &
                  taucld    , ssacld   , asmcld   , fsfcld   , &
                  cicewp    , cliqwp   , reice    , reliq    , &
                  tauaer    , ssaaer   , asmaer   , ecaer    , &
                  swuflx    , swdflx   , swhr     , swuflxc  , swdflxc, swhrc)
          endif
          
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
!!$                  zeros     , o3      , co2     , zeros, zeros, zeros, &
                  zeros     , zeros   , zeros   , zeros, emis , &
                  inflglw   , iceflglw, liqflglw, cldfr,  &
                  taucld    , cicewp  , cliqwp  , reice, reliq, &
                  tauaer    , &
                  uflx      , dflx    , hr      , uflxc, dflxc, hrc)
          endif

          lwijk   = reshape(hr(:,sk:1:-1),(/ si/lonstep,sj,sk /))*daypersec

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
             swflxijk = reshape(swdflx(:,1)-swuflx(:,1),(/ si/lonstep,sj /))
!!$             lwflxijk = reshape(  uflx(:,1)-  dflx(:,1),(/ si/lonstep,sj /))
             lwflxijk = reshape(  dflx(:,1)            ,(/ si/lonstep,sj /))
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
          
          call write_diag_rrtm(Time,is,js)
        end subroutine run_rrtmg

!*****************************************************************************************
        subroutine write_diag_rrtm(Time,is,js)
          use rrtm_vars,only: sw_flux,lw_flux,zencos,tdt_rad,tdt_sw_rad,tdt_lw_rad&
               &,id_tdt_rad,id_tdt_sw,id_tdt_lw,id_coszen,id_flux_sw,id_flux_lw
          use diag_manager_mod, only: register_diag_field, send_data
          use time_manager_mod,only: time_type
          implicit none
          type(time_type),intent(in) :: Time
          integer, intent(in)        :: is, js
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
        end subroutine write_diag_rrtm
!*****************************************************************************************

        subroutine rrtm_radiation_end
          use rrtm_vars
          use interpolator_mod, only: interpolator_end
          implicit none

          if(do_read_ozone)call interpolator_end(o3_interp)
        end subroutine rrtm_radiation_end


!*****************************************************************************************
      end module rrtm_gases
       


        
