module cloud_generator_mod

use fms_mod, only: error_mesg, FATAL, write_version_number

use random_numbers_mod, only: randomNumberStream
!--------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloud_generator.f90,v 1.1 2004/08/10 18:48:55 fms Exp $'
character(len=128)  :: tagname =  '$Name: lima $'

!---------------------------------------------------------------------
  public :: cloud_generator_init, &
            cloud_generator_end,  &
            generate_stochastic_clouds, &
            do_cloud_generator,   &
            compute_overlap_weighting
  
!---------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.  ! module is initialized ?

!----------------------------------------------------------------------

                              contains 
                              
!######################################################################
subroutine cloud_generator_init


!---------------------------------------------------------------------
!    cloud_generator_init is the constructor for 
!    cloud_generator_mod.

!----------------------------------------------------------------------

     if (.not. module_is_initialized) then
       call write_version_number ('Null module: '//version, tagname)
       module_is_initialized = .true.
     end if

     call error_mesg('subroutine cloud_generator_init in cloud_generator_mod', &
     'This module is not supported as part of the public release', FATAL)

end subroutine cloud_generator_init
!----------------------------------------------------------------------
subroutine generate_stochastic_clouds(streams, ql, qi, qa,         &
                                      overlap, pFull, temperature, &
                                      cld_thickness, ql_stoch, qi_stoch, qa_stoch)
!--------------------------------------------------------------------
!   intent(in) variables:
!
  type(randomNumberStream), &
           dimension(:, :),     intent(inout) :: streams
  ! Dimension nx, ny, nz
  real,    dimension(:, :, :),    intent( in) :: ql, qi, qa
  integer,                     optional, &
                                  intent( in) :: overlap
  real,    dimension(:, :, :), optional, &
                                 intent( in)  :: pFull, temperature
  ! Dimension nx, ny, nz, nCol = nBands
  integer, dimension(:, :, :, :), intent(out) :: cld_thickness 
  real,    dimension(:, :, :, :), intent(out) :: ql_stoch, &
                                                 qi_stoch, qa_stoch
  ! ---------------------------------------------------------

  call error_mesg('subroutine generate_stochastic_clouds in cloud_generator_mod', &
  'This module is not supported as part of the public release', FATAL)
           
end subroutine generate_stochastic_clouds
!----------------------------------------------------------------------

function compute_overlap_weighting(qaPlus, qaMinus, pPlus, pMinus) result(weighting)
  real, dimension(:, :), intent( in) :: qaPlus, qaMinus, pPlus, pMinus
  real, dimension(size(pPlus,1),size(pPlus,2)) :: weighting
        
  call error_mesg("function compute_overlap_weighting in cloud_generator_mod", &
  'This module is not supported as part of the public release', FATAL)

  weighting = 0. ! This line of code exists only to prevent compiler warnings
      
end function compute_overlap_weighting

! ---------------------------------------------------------

 subroutine cloud_generator_end       
          
   call error_mesg ('subroutine cloud_generator_end in cloud_generator_mod',   &
   'This module is not supported as part of the public release', FATAL)
        
 end subroutine cloud_generator_end
!--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !
  !  Function to report if the cloud generator is being used.
  !
  function do_cloud_generator()
    logical :: do_cloud_generator
    
    do_cloud_generator = .false.
  end function do_cloud_generator
  !--------------------------------------------------------------------

end module cloud_generator_mod
