      module rrtm_astro
        use parkind, only : im => kind_im, rb => kind_rb
        implicit none

!to be included in namelist
!
        real(kind=rb)    :: obliq=23.439

        contains
!--------------------------------------------------------------------------------------
! parts of this are taken from GFDL's astronomy.f90      
          subroutine compute_zenith(Time, dt, lat, lon, cosz, dyofyr)
            use time_manager_mod,only: time_type,get_time,length_of_year
            use constants_mod, only: PI
            implicit none

            type(time_type),             intent(in) :: Time
            integer(kind=im),            intent(in) :: dt
            real(kind=rb),dimension(:,:),intent(in) :: lat,lon
            real(kind=rb),dimension(:,:),intent(out):: cosz !,fracday
            integer(kind=im)            ,intent(out):: dyofyr
            
            real(kind=rb),dimension(size(lat,1),size(lat,2)) :: h,cos_h, &
                 lat_h

            real dec_sin,dec_tan,dec,dec_cos,twopi,dt_pi

            integer  :: seconds,sec2,days,daysperyear
            real,dimension(size(lon,1),size(lon,2)) :: time_pi,aa,bb,tt,st,stt,sh
            real     :: radsec,radday

            integer  :: i,j
            real     :: radpersec, radperday
            real     :: eps=1.0E-05,deg2rad
            
            deg2rad = PI/180.
            twopi = 2*PI

            call get_time(length_of_year(),sec2,daysperyear)

            radpersec=2*PI/86400.
            radperday=2*PI/daysperyear

            !get the time for origin
            call get_time(Time,seconds,days)
            !convert into radians
            radsec = seconds*radpersec
            dt_pi  = dt*radpersec
            
            !set local time throughout the globe
            do i=1,size(lon,1)
               !move it into interval [-PI,PI]
               time_pi(i,:) = modulo(radsec + lon(i,:),2*PI) - PI
            enddo
            where(time_pi >= PI) time_pi = time_pi - twopi
            where(time_pi < -PI) time_pi = time_pi + twopi 
            !time_pi now contains local time at each grid point
            !get day of the year relative to equinox. We set equinox at 0.25,0.75*daysperyear
            days = days - int(0.25*daysperyear)
            dyofyr   = modulo(days,daysperyear) !NH winter solstice at day 0
            !convert into radians
            radday = dyofyr*radperday
            !get declination
            dec_sin = sin(obliq*deg2rad)*sin(radday) !check sign ("-" in GFDL's code)
            dec = asin(dec_sin)
            dec_cos = cos(dec)

            !now compute the half day, to determine if it's day or night
            dec_tan = tan(dec)
            lat_h = lat
            where(lat_h == 0.5*PI) lat_h = lat - eps
            where(lat_h ==-0.5*PI) lat_h = lat + eps
            cos_h = -tan(lat_h)*dec_tan
            where(cos_h <= -1.0) h = PI
            where(cos_h >=  1.0) h = 0.0
            where(cos_h > -1.0 .and. cos_h < 1.0)&
                 h = acos(cos_h)
!---------------------------------------------------------------------
!    define terms needed in the cosine zenith angle equation.
!--------------------------------------------------------------------
            aa = sin(lat)*dec_sin
            bb = cos(lat)*dec_cos

            !finally, compute the zenith angle
            if(dt > 0)then !average over some given time interal dt
               tt  = time_pi + dt_pi
               st  = sin(time_pi)
               stt = sin(tt)
               sh  = sin(h)
               cosz = 0.0
!-------------------------------------------------------------------
!    case 1: entire averaging period is before sunrise.
!-------------------------------------------------------------------
               where (time_pi < -h .and. tt < -h) cosz = 0.0

!-------------------------------------------------------------------
!    case 2: averaging period begins before sunrise, ends after sunrise
!    but before sunset
!-------------------------------------------------------------------
               where ( (tt+h) /= 0.0 .and.   time_pi < -h .and. abs(tt) <= h)   &
                    cosz = aa + bb*(stt + sh)/ (tt + h)
!!$                    cosz=2.
!-------------------------------------------------------------------
!    case 3: averaging period begins before sunrise, ends after sunset,
!    but before the next sunrise. modify if averaging period extends 
!    past the next day's sunrise, but if averaging period is less than 
!    a half- day (pi) that circumstance will never occur.
!-------------------------------------------------------------------
               where (time_pi < -h .and. h /= 0.0 .and. h < tt)    &
                    cosz = aa + bb*( sh + sh)/(h+h)
!!$                    cosz=3.
!-------------------------------------------------------------------
!    case 4: averaging period begins after sunrise, ends before sunset.
!-------------------------------------------------------------------
               where ( abs(time_pi) <= h .and. abs(tt) <= h)    &
                    cosz = aa + bb*(stt - st)/ (tt - time_pi)
!!$                    cosz=4.
!-------------------------------------------------------------------
!    case 5: averaging period begins after sunrise, ends after sunset. 
!    modify when averaging period extends past the next day's sunrise.  
!------------------------------------------------------------------- 
               where ((h-time_pi) /= 0.0 .and. abs(time_pi) <= h .and.  h < tt)    &
                    cosz = aa + bb*(sh - st)/(h-time_pi)
!!$                    cosz=5.
!-------------------------------------------------------------------
!    case 6: averaging period begins after sunrise , ends after the
!    next day's sunrise. note that this includes the case when the
!    day length is one day (h = pi).
!-------------------------------------------------------------------
               where (twopi - h < tt .and. (tt+h-twopi) /= 0.0 .and. time_pi <= h ) &
                    cosz = (cosz*(h - time_pi) + (aa*(tt + h - twopi) +     &
                    bb*(stt + sh))) / ((h - time_pi) + (tt + h - twopi))
!!$                    cosz=6.

!-------------------------------------------------------------------
!    case 7: averaging period begins after sunset and ends before the
!    next day's sunrise
!-------------------------------------------------------------------
               where(  h <  time_pi .and. twopi - h >= tt  ) cosz = 0.0

!-------------------------------------------------------------------
!    case 8: averaging period begins after sunset and ends after the
!    next day's sunrise but before the next day's sunset. if the
!    averaging period is less than a half-day (pi) the latter
!    circumstance will never occur.
!-----------------------------------------------------------------
               where(  h <  time_pi .and. twopi - h < tt  ) &
                  cosz = aa + bb*(stt + sh) / (tt + h - twopi)
!!$                  cosz=8.0

!mj-----------------------------------------------------------------
!    case 9: averaging period begins after sunset and ends after the
!    next day's sunset but before the day after's sunrise. Typically,
!    this is daily average.
!mj-----------------------------------------------------------------
               where( h < time_pi .and. h + twopi < tt .and. h /= 0. )
                  cosz = aa + bb*(sh + sh) / (h + h)
!!$                  cosz = 9.0
               end where
    
!----------------------------------------------------------------------
!    if instantaneous values are desired, define cosz at time t.
!----------------------------------------------------------------------
            else !version w/o time averaging
               where(abs(time_pi) <= h)
                  cosz    = aa + bb*cos(time_pi)
               elsewhere
                  cosz    = 0.0
               endwhere
            endif
!----------------------------------------------------------------------
!    be sure that cosz is not negative.
!----------------------------------------------------------------------
          cosz = max(0.0, cosz)

          end subroutine compute_zenith


!!$        subroutine get_astronomy(lat, lon, gmt, time_since_ae, cosz, &
!!$             fracday, rrsun, dt) 
!!$          use constants_mod, only: PI
!!$
!!$!--stolen from astronomy.f90 - diurnal_solar_2d-----------------------
!!$
!!$          real, dimension(:,:), intent(in)           :: lat, lon
!!$          real,                 intent(in)           :: gmt, time_since_ae
!!$          real, dimension(:,:), intent(out)          :: cosz, fracday
!!$          real,                 intent(out)          :: rrsun
!!$          real,                 intent(in), optional :: dt
!!$
!!$          real, dimension(size(lat,1),size(lat,2)) :: t, tt, h, aa, bb,  &
!!$               st, stt, sh
!!$          real                                     :: twopi,ang, dec
!!$          twopi = 2.*PI
!!$
!!$!---------------------------------------------------------------------
!!$!    define the orbital angle (location in year), solar declination and
!!$!    earth sun distance factor. use functions contained in this module.
!!$!---------------------------------------------------------------------
!!$          !ang = angle(time_since_ae)
!!$          ang = time_since_ae
!!$          call declination_loc(ang,dec)
!!$          rrsun  = 1.0!r_inv_squared(ang)
!!$          
!!$!---------------------------------------------------------------------
!!$!    define terms needed in the cosine zenith angle equation.
!!$!--------------------------------------------------------------------
!!$          aa = sin(lat)*sin(dec)
!!$          bb = cos(lat)*cos(dec)
!!$
!!$!---------------------------------------------------------------------
!!$!    define local time. force it to be between -pi and pi.
!!$!--------------------------------------------------------------------
!!$          t = gmt + lon - PI
!!$          where(t >= PI) t = t - twopi  
!!$          where(t < -PI) t = t + twopi   
!!$      
!!$!---------------------------------------------------------------------
!!$!    perform a time integration to obtain cosz and fracday if desired.
!!$!    output is valid over the period from t to t + dt.
!!$!--------------------------------------------------------------------
!!$          h   = half_day_loc(lat,dec)
!!$          if ( present(dt) ) then   ! (perform time averaging)
!!$             tt = t + dt
!!$             st  = sin(t)
!!$             stt = sin(tt)
!!$             sh  = sin(h)
!!$             cosz = 0.0
!!$
!!$!-------------------------------------------------------------------
!!$!    case 1: entire averaging period is before sunrise.
!!$!-------------------------------------------------------------------
!!$             where (t < -h .and. tt < -h) cosz = 0.0
!!$
!!$!-------------------------------------------------------------------
!!$!    case 2: averaging period begins before sunrise, ends after sunrise
!!$!    but before sunset
!!$!-------------------------------------------------------------------
!!$             where ( (tt+h) /= 0.0 .and.   t < -h .and. abs(tt) <= h)   &
!!$                  cosz = aa + bb*(stt + sh)/ (tt + h)
!!$
!!$!-------------------------------------------------------------------
!!$!    case 3: averaging period begins before sunrise, ends after sunset,
!!$!    but before the next sunrise. modify if averaging period extends 
!!$!    past the next day's sunrise, but if averaging period is less than 
!!$!    a half- day (pi) that circumstance will never occur.
!!$!-------------------------------------------------------------------
!!$             where (t < -h .and. h /= 0.0 .and. h < tt)    &
!!$                  cosz = aa + bb*( sh + sh)/(h+h)
!!$
!!$!-------------------------------------------------------------------
!!$!    case 4: averaging period begins after sunrise, ends before sunset.
!!$!-------------------------------------------------------------------
!!$             where ( abs(t) <= h .and. abs(tt) <= h)    &
!!$                  cosz = aa + bb*(stt - st)/ (tt - t)
!!$             
!!$!-------------------------------------------------------------------
!!$!    case 5: averaging period begins after sunrise, ends after sunset. 
!!$!    modify when averaging period extends past the next day's sunrise.  
!!$!------------------------------------------------------------------- 
!!$             where ((h-t) /= 0.0 .and. abs(t) <= h .and.  h < tt)    &
!!$                  cosz = aa + bb*(sh - st)/(h-t)
!!$
!!$!-------------------------------------------------------------------
!!$!    case 6: averaging period begins after sunrise , ends after the
!!$!    next day's sunrise. note that this includes the case when the
!!$!    day length is one day (h = pi).
!!$!-------------------------------------------------------------------
!!$             where (twopi - h < tt .and. (tt+h-twopi) /= 0.0 .and. t <= h ) &
!!$                  cosz = (cosz*(h - t) + (aa*(tt + h - twopi) +     &
!!$                  bb*(stt + sh))) / ((h - t) + (tt + h - twopi))
!!$
!!$!-------------------------------------------------------------------
!!$!    case 7: averaging period begins after sunset and ends before the
!!$!    next day's sunrise
!!$!-------------------------------------------------------------------
!!$             where(  h <  t .and. twopi - h >= tt  ) cosz = 0.0
!!$
!!$!-------------------------------------------------------------------
!!$!    case 8: averaging period begins after sunset and ends after the
!!$!    next day's sunrise but before the next day's sunset. if the
!!$!    averaging period is less than a half-day (pi) the latter
!!$!    circumstance will never occur.
!!$!-----------------------------------------------------------------
!!$             where(  h <  t .and. twopi - h < tt  ) 
!!$                cosz = aa + bb*(stt + sh) / (tt + h - twopi)
!!$             end where
!!$
!!$!-------------------------------------------------------------------
!!$!    day fraction is the fraction of the averaging period contained 
!!$!    within the (-h,h) period.
!!$!-------------------------------------------------------------------
!!$             where (t < -h .and.      tt < -h)      fracday = 0.0
!!$             where (t < -h .and. abs(tt) <= h)      fracday = (tt + h )/dt
!!$             where (t < -h .and.       h < tt)      fracday = ( h + h )/dt
!!$             where (abs(t) <= h .and. abs(tt) <= h) fracday = (tt - t )/dt
!!$             where (abs(t) <= h .and.       h < tt) fracday = ( h - t )/dt
!!$             where (      h <  t                 )  fracday = 0.0
!!$             where (twopi - h < tt)                 fracday = fracday +  &
!!$                  (tt + h - &
!!$                  twopi)/dt      
!!$!----------------------------------------------------------------------
!!$!    if instantaneous values are desired, define cosz at time t.
!!$!----------------------------------------------------------------------
!!$          else  ! (no time averaging)
!!$             where (abs(t) < h)
!!$                cosz = aa + bb*cos(t)
!!$                fracday = 1.0
!!$             else where
!!$                cosz = 0.0
!!$                fracday = 0.0
!!$             end where
!!$          end if
!!$
!!$!----------------------------------------------------------------------
!!$!    be sure that cosz is not negative.
!!$!----------------------------------------------------------------------
!!$          cosz = max(0.0, cosz)
!!$
!!$        end subroutine get_astronomy
!!$
!!$
!!$        function half_day_loc(latitude, dec) result(h)
!!$          use constants_mod, only: PI
!!$!copied from astronomy.f90
!!$          real, dimension(:,:), intent(in)                     :: latitude
!!$          real,                 intent(in)                     :: dec
!!$          real, dimension(size(latitude,1),size(latitude,2))   :: h
!!$
!!$          real, dimension (size(latitude,1),size(latitude,2)):: & 
!!$               cos_half_day, lat
!!$          real :: tan_dec 
!!$          real :: eps = 1.0E-05
!!$   
!!$!--------------------------------------------------------------------
!!$!    define tangent of the declination.
!!$!--------------------------------------------------------------------
!!$          tan_dec = tan(dec)
!!$
!!$!--------------------------------------------------------------------
!!$!    adjust latitude so that its tangent will be defined.
!!$!--------------------------------------------------------------------
!!$          lat = latitude
!!$          where (latitude ==  0.5*PI) lat= latitude - eps
!!$          where (latitude == -0.5*PI) lat= latitude + eps
!!$
!!$!--------------------------------------------------------------------
!!$!    define the cosine of the half-day length. adjust for cases of 
!!$!    all daylight or all night.
!!$!--------------------------------------------------------------------
!!$          cos_half_day = -tan(lat)*tan_dec
!!$          where (cos_half_day <= -1.0)  h = PI
!!$          where (cos_half_day >= +1.0)  h = 0.0
!!$          where(cos_half_day > -1.0 .and. cos_half_day < 1.0) &
!!$                                               h = acos(cos_half_day)
!!$
!!$
!!$        end function half_day_loc
!!$!
!!$        subroutine declination_loc(ang,declination)
!!$          use rrtm_vars
!!$
!!$!---stolen from astronomy.f90 ---------------------------------------
!!$          real, intent(in) :: ang
!!$!--------------------------------------------------------------------
!!$
!!$!--------------------------------------------------------------------
!!$!   local variables
!!$
!!$          real :: declination
!!$          real :: rad_obliq, sin_dec
!!$          real :: obliq=23.439
!!$          rad_obliq   =   obliq*deg2rad
!!$          sin_dec     = - sin(rad_obliq)*sin(ang)
!!$          declination =   asin(sin_dec)
!!$
!!$
!!$        end subroutine declination_loc


        end module rrtm_astro
