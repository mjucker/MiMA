!FDOC_TAG_GFDL

                 module donner_deep_clouds_W_mod
! <CONTACT EMAIL="fei.liu@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!          donner deep cloud radiative properties module
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use time_manager_mod,       only: time_type
use donner_deep_mod,        only: donner_deep_avg, donner_deep_init
use       fms_mod,          only: open_namelist_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  close_file, FATAL, NOTE, &
                                  WARNING, mpp_pe, mpp_root_pe, &
                                  write_version_number, stdlog
use rad_utilities_mod,      only: longwave_control_type, Lw_control, &
                                  shortwave_control_type, Sw_control,&
                                  microphysics_type,  &
                                  microrad_properties_type, &
                                  cld_specification_type, &
                                  cloudrad_control_type, Cldrad_control

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!          donner deep cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: donner_deep_clouds_W.f90,v 12.0 2005/04/14 15:45:00 fms Exp $'
   character(len=128)  :: tagname =  '$Name: lima $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          donner_deep_clouds_W_init, donner_deep_clouds_calc,  &
          donner_deep_clouds_W_end , donner_deep_clouds_amt

!---------------------------------------------------------------------
!-------- namelist  ---------

logical   :: using_dge_lw = .true.
logical   :: using_dge_sw = .true.



namelist /donner_deep_clouds_W_nml /     &
       using_dge_sw, using_dge_lw


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


  logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





! <SUBROUTINE NAME="donner_deep_clouds_W_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_W_init  (pref, lonb, latb, axes, Time)
!		
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
! 
!  </IN>
!  <IN NAME="axes" TYPE="integer">
! 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_W_init  (pref, lonb, latb, axes, Time)

real, dimension(:,:), intent(in) :: pref
real, dimension(:), intent(in) :: lonb, latb
integer, dimension(4), intent(in)      :: axes
type(time_type),       intent(in)      :: Time

      integer            :: unit, ierr, io
      integer            :: ix, jx, kx

     if (module_is_initialized) return
!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=donner_deep_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'donner_deep_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif

      if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
         write (stdlog(),nml=donner_deep_clouds_W_nml)
      endif

       ix = size(lonb,1)-1
       jx = size(latb,1)-1
       kx = size(pref,1) - 1

!---------------------------------------------------------------------
       call donner_deep_init(lonb, latb, pref(:,1), axes, Time)

       module_is_initialized = .true.


end subroutine donner_deep_clouds_W_init

! <SUBROUTINE NAME="donner_deep_clouds_W_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_W_end
       
!----------------------------------------------------------------------
!    diag_clouds_W_end is the destructor for diag_clouds_W_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!--------------------------------------------------------------------


end subroutine donner_deep_clouds_W_end


!#################################################################

! <SUBROUTINE NAME="donner_deep_clouds_amt2">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_amt2  (             &
!		is, ie, js, je,    &
!		cld_cell,        &
!		Cell_microphys,  &
!		cld_meso,        &
!		Meso_microphys  )
!		
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <OUT NAME="cld_cell" TYPE="real">
! 
!  </OUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
!  <OUT NAME="cld_meso" TYPE="real">
! 
!  </OUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_amt2  (             &
                       is, ie, js, je,    &
                                cld_cell,        &
                          Cell_microphys,  &
                       cld_meso,        &
                        Meso_microphys  )


integer, intent(in) :: is,ie,js,je
type(microphysics_type), intent(inout) :: Cell_microphys, Meso_microphys
real, dimension(:,:,:), intent(  out) :: cld_cell, cld_meso


!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------



!---------------------------------------------------------------------
!     obtain the deep cloud areas and properties
!---------------------------------------------------------------------


      call donner_deep_avg(  is, ie, js, je,           &
                             cell_cloud_frac_out  = cld_cell,        &
                      cell_liquid_amt_out  = Cell_microphys%conc_drop, &
                 cell_liquid_size_out = Cell_microphys%size_drop,&
               cell_ice_amt_out     = Cell_microphys%conc_ice,    &
              cell_ice_size_out    = Cell_microphys%size_ice,   &
                             meso_cloud_frac_out  = cld_meso,        &
               meso_liquid_amt_out  = Meso_microphys%conc_drop , &
               meso_liquid_size_out = Meso_microphys%size_drop,&
               meso_ice_amt_out     = Meso_microphys%conc_ice,   &
            meso_ice_size_out    = Meso_microphys%size_ice    )


end subroutine donner_deep_clouds_amt2 

!---------------------------------------------------------------------

! <SUBROUTINE NAME="donner_deep_clouds_amt">
!  <OVERVIEW>
!    donner_deep_clouds_amt defines the distribution of cloud water and
!    cloud ice concentration and particle size and total cloud fraction
!    in both the mesoscale and convective cell-scale components of the
!    clouds associated with donner_deep convection. these values will
!    be combined with the large-scale cloud fields to produce the dist-
!    ribution of cloud radiative properties that will be seen by the
!    radiation package.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    donner_deep_clouds_amt defines the distribution of cloud water and
!    cloud ice concentration and particle size and total cloud fraction
!    in both the mesoscale and convective cell-scale components of the
!    clouds associated with donner_deep convection. these values will
!    be combined with the large-scale cloud fields to produce the dist-
!    ribution of cloud radiative properties that will be seen by the
!    radiation package.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_amt (is, ie, js, je, Cell_microphys,  &
!		Meso_microphys)
!		
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_amt (is, ie, js, je, Cell_microphys,  &
                                   Meso_microphys)

!---------------------------------------------------------------------
!    donner_deep_clouds_amt defines the distribution of cloud water and
!    cloud ice concentration and particle size and total cloud fraction
!    in both the mesoscale and convective cell-scale components of the
!    clouds associated with donner_deep convection. these values will
!    be combined with the large-scale cloud fields to produce the dist-
!    ribution of cloud radiative properties that will be seen by the
!    radiation package.
!----------------------------------------------------------------------

integer,                 intent(in)    :: is,ie,js,je
type(microphysics_type), intent(inout) :: Cell_microphys, Meso_microphys

!---------------------------------------------------------------------
!     call donner_deep_avg to obtain the specification fields for both
!     the mesoscale and convective cellscale clouds assocated with 
!     donner_deep convection.
!---------------------------------------------------------------------
      call donner_deep_avg (                           &
                      is, ie, js, je,           &
                      cell_cloud_frac_out  = Cell_microphys%cldamt, &
                      cell_liquid_amt_out  = Cell_microphys%conc_drop, &
                      cell_liquid_size_out = Cell_microphys%size_drop,&
                      cell_ice_amt_out     = Cell_microphys%conc_ice, &
                      cell_ice_size_out    = Cell_microphys%size_ice, &
                      meso_cloud_frac_out  = Meso_microphys%cldamt,   &
                      meso_liquid_amt_out  = Meso_microphys%conc_drop, &
                      meso_liquid_size_out = Meso_microphys%size_drop,&
                      meso_ice_amt_out     = Meso_microphys%conc_ice, &
                      meso_ice_size_out    = Meso_microphys%size_ice )

!---------------------------------------------------------------------



end subroutine donner_deep_clouds_amt  


!#####################################################################


! <SUBROUTINE NAME="donner_deep_clouds_calc">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_calc (             &
!		is,ie,js,je,deltaz,press,temp,                 &
!		cld_cell,               &
!		cldext_cell, cldsct_cell, cldasymm_cell,  &
!		abscoeff_cell,          &
!		cld_meso,               &
!		cldext_meso, cldsct_meso, cldasymm_meso,  &
!		abscoeff_meso)
!		
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <IN NAME="deltaz" TYPE="real">
! 
!  </IN>
!  <IN NAME="press" TYPE="real">
! 
!  </IN>
!  <IN NAME="temp" TYPE="real">
! 
!  </IN>
!  <INOUT NAME="cld_cell" TYPE="real">
! 
!  </INOUT>
!  <OUT NAME="cldext_cell" TYPE="real">
! 
!  </OUT>
!  <OUT NAME="cldsct_cell" TYPE="real">
! 
!  </OUT>
!  <OUT NAME="cldasymm_cell" TYPE="real">
! 
!  </OUT>
!  <OUT NAME="abscoeff_cell" TYPE="real">
! 
!  </OUT>
!  <INOUT NAME="cld_meso" TYPE="real">
! 
!  </INOUT>
!  <OUT NAME="cldext_meso" TYPE="real">
! 
!  </OUT>
!  <OUT NAME="cldsct_meso" TYPE="real">
! 
!  </OUT>
!  <OUT NAME="cldasymm_meso" TYPE="real">
! 
!  </OUT>
!  <OUT NAME="abscoeff_meso" TYPE="real">
! 
!  </OUT>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_calc (             &
                  is,ie,js,je,deltaz,press,temp,                 &
                                    cld_cell,               &
 cldext_cell, cldsct_cell, cldasymm_cell,  &
    abscoeff_cell,          &
                                    cld_meso,               &
 cldext_meso, cldsct_meso, cldasymm_meso,  &
    abscoeff_meso)


integer, intent(in) :: is,ie,js,je
real, dimension(:,:,:), intent(in) :: deltaz, press, temp
real, dimension(:,:,:), intent(inout) :: cld_cell, cld_meso
real, dimension(:,:,:,:), intent(out) :: cldext_cell, cldsct_cell,  &
                                         cldasymm_cell, abscoeff_cell
real, dimension(:,:,:,:), intent(out) :: cldext_meso, cldsct_meso,  &
                                         cldasymm_meso, abscoeff_meso


!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------

integer   :: idim, jdim, kdim
real, dimension(size(cld_cell,1),size(cld_cell,2),size(cld_cell,3)) :: &
      cell_liquid_amt, cell_ice_amt, cell_liquid_size, cell_ice_size,  &
      meso_liquid_amt, meso_ice_amt, meso_liquid_size, meso_ice_size
integer :: unit

     idim = size(cld_cell,1)
     jdim = size(cld_cell,2)
     kdim = size(cld_cell,3)

!--------------------------------------------------------------------


!---------------------------------------------------------------------
!     obtain the deep cloud areas and properties
!---------------------------------------------------------------------


      call donner_deep_avg(  is, ie, js, je,           &
                             cell_cloud_frac_out  = cld_cell,        &
                             cell_liquid_amt_out  = cell_liquid_amt, &
                             cell_liquid_size_out = cell_liquid_size,&
                             cell_ice_amt_out     = cell_ice_amt,    &
                             cell_ice_size_out    = cell_ice_size,   &
                             meso_cloud_frac_out  = cld_meso,        &
                             meso_liquid_amt_out  = meso_liquid_amt, &
                             meso_liquid_size_out = meso_liquid_size,&
                             meso_ice_amt_out     = meso_ice_amt,    &
                             meso_ice_size_out    = meso_ice_size    )
 
!      unit = open_file ('fort.149', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' donner_deep_clouds_w'
! write (unit,*) idim,jdim,kdim
! write (unit,*) ' jabs = ', jabs
!        write (unit,*) ' cld_cell'
! write (unit,*) cld_cell
!        write (unit,*) ' cell_liquid_amt'
! write (unit,*) cell_liquid_amt
!        write (unit,*) ' cell_liquid_size'
! write (unit,*) cell_liquid_size
!        write (unit,*) ' cell_ice_amt'
! write (unit,*) cell_ice_amt
!        write (unit,*) ' cell_ice_size'
! write (unit,*) cell_ice_size
!      call close_file (unit)
  
!
!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------

!         call microphys_rad_driver(                         &
!                 is,ie,js,je,deltaz,press,temp,                  &
!                   conc_drop_in=cell_liquid_amt,        &
!   conc_ice_in=cell_ice_amt,            &
!   size_drop_in=cell_liquid_size,       &
!   size_ice_in=cell_ice_size,           &
!   cldext=cldext_cell,               &
!   cldsct=cldsct_cell,               &
!   cldasymm=cldasymm_cell,           &
!   abscoeff=abscoeff_cell,           &
!   using_dge_sw=using_dge_sw,        &
!   using_dge_lw=using_dge_lw)

!         call microphys_rad_driver(                         &
!                 is,ie,js,je,deltaz,press,temp,                  &
!                   conc_drop_in=meso_liquid_amt,        &
!   conc_ice_in=meso_ice_amt,            &
!   size_drop_in=meso_liquid_size,       &
!   size_ice_in=meso_ice_size,           &
!   cldext=cldext_meso,               &
!   cldsct=cldsct_meso,               &
!   cldasymm=cldasymm_meso,           &
!   abscoeff=abscoeff_meso,           &
!   using_dge_sw=using_dge_sw,        &
!   using_dge_lw=using_dge_lw)



 

!
!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------
!         call microphys_rad_driver(                         &
!                 is,ie,js,je,deltaz,press,temp,                  &
!                   conc_drop_in=cell_liquid_amt,        &
!   conc_ice_in=cell_ice_amt,            &
!   size_drop_in=cell_liquid_size,       &
!   size_ice_in=cell_ice_size,           &
!   cldext=cldext_cell,               &
!   cldsct=cldsct_cell,               &
!   cldasymm=cldasymm_cell,           &
!   abscoeff=abscoeff_cell,           &
!   using_dge_sw=using_dge_sw,        &
!   using_dge_lw=using_dge_lw)

!         call microphys_rad_driver(                         &
!                 is,ie,js,je,deltaz,press,temp,                  &
!                   conc_drop_in=meso_liquid_amt,        &
!   conc_ice_in=meso_ice_amt,            &
!   size_drop_in=meso_liquid_size,       &
!   size_ice_in=meso_ice_size,           &
!   cldext=cldext_meso,               &
!   cldsct=cldsct_meso,               &
!   cldasymm=cldasymm_meso,           &
!   abscoeff=abscoeff_meso,           &
!   using_dge_sw=using_dge_sw,        &
!   using_dge_lw=using_dge_lw)

  





end subroutine donner_deep_clouds_calc


!####################################################################


       end module donner_deep_clouds_W_mod



