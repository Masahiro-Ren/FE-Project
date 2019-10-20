!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_fieldutil
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_const, only: &
    PI => CONST_PI

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  public :: fieldutil_get_profile1d_tracer
  public :: fieldutil_get_upwind_pos1d
  
  public :: fieldutil_get_profile2d_tracer
  public :: fieldutil_get_profile2d_flow

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !


  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------
contains  

!--- 1d -----------------------------------------------------

function fieldutil_get_upwind_pos1d(pos, ADV_VEL, nowtime, dom_min, dom_max) result(upos)
  real(RP), intent(in) :: pos(:)
  real(RP), intent(in) :: ADV_VEL
  real(RP), intent(in) :: nowtime
  real(RP), intent(in) :: dom_min, dom_max
  real(RP) :: upos(size(pos))

  integer :: period
  !-------

  period = ADV_VEL*nowtime/(dom_max - dom_min)
  upos(:) = pos(:) - (ADV_VEL*nowtime - dble(period)*(dom_max - dom_min))
  where (upos < dom_min)
    upos = dom_max + (upos - dom_min)
  end where
end function fieldutil_get_upwind_pos1d

!> Get 1D data whose value is set based on specified pattern. 
!
! Assume that range of domain is 0 <= x <= 1.
! If the profile_name is 'sin', param1 indicates the wavenumber. 
! If the profile_name is 'cosbell', param1 indicates the half of width. 
! If the profile_name is 'top-hat', param1 indicates the half of width. 
!
subroutine fieldutil_get_profile1d_tracer( profile,  &  
  profile_name, x, params, N ) 
  implicit none
  !--------------------------------
  integer, intent(in) :: N
  real(RP), intent(out) :: profile(N)
  character(*), intent(in) :: profile_name
  real(RP), intent(in) :: x(N)
  real(RP), intent(in) :: params(:)

  real(RP) :: dist(N)
  !--------------------------------

  profile(:) = 0.0_RP

  select case(profile_name)
  case ('constant')
    profile(:) = params(1)    
  case ('sin')
    profile(:) = sin( params(1)*2.0_RP*PI*x(:) )
  case ('gaussian-hill')
    profile(:) = exp( - 0.5_RP * ((x(:) - params(1))**2)/params(2)**2 )
  case ('cosine-bell')
    dist(:) = abs(x(:) - params(1))
    where( dist <= params(2) )
      profile(:) = (1.0_RP + cos(PI*dist/params(2)))*0.5_RP
    end where
  case ('top-hat')
    where( abs(x - params(1)) <= params(2) )
      profile(:) = 1.0_RP
    end where
  case default
    LOG_ERROR('fieldutil_get_profile2d',*) trim(profile_name)//' is not supported. Check!'
    call PRC_abort    
  end select

  return
end subroutine fieldutil_get_profile1d_tracer

!-- 2d--------------------------------------------------------------

!> Get 2D data whose value is set based on specified pattern. 
!
! Assume that range of domain is 0 <= x <= 1.
! If the profile_name is 'sin', param1 and param2 are the wavenumbers in x- and y- directions, respectively. 
! If the profile_name is 'gaussian-hill', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
! If the profile_name is 'cosine-bell', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
! If the profile_name is 'top-hat', (param1,param2) is the coordinate of center position, param3 and param4 are the half of width in x- and y- directions, respectively. 
!
subroutine fieldutil_get_profile2d_tracer( profile, &
  profile_name, x, y, params, N )
  implicit none
  !--------------------------------
  integer, intent(in) :: N
  real(RP), intent(out) :: profile(N)
  character(*), intent(in) :: profile_name
  real(RP), intent(in) :: x(N)
  real(RP), intent(in) :: y(N)
  real(RP), intent(in) :: params(:)

  real(RP) :: dist(N)
  !------------------------------------------------------------------------

  profile(:) = 0.0_RP

  select case(profile_name)
  case ('constant')
    profile(:) = params(1)
  case ('sin')
    profile(:) = sin( params(1)*2.0_RP*PI*x(:) ) * sin( params(2)*2.0_RP*PI*y(:) )
  case ('gaussian-hill')
    profile(:) = exp( - 0.5_RP * ((x(:) - params(1))**2 + (y(:) - params(2))**2)/params(3)**2 )
  case ('cosine-bell')
    dist(:) = sqrt( (x(:) - params(1))**2 + (y(:) - params(2))**2 )
    where( dist <= params(3) )
      profile(:) = (1.0_RP + cos(PI*dist(:)/params(3)))*0.5_RP
    end where
  case ('top-hat')
    where( abs(x-params(1)) <= params(3) .and. abs(y-params(2)) <= params(4) )
      profile(:) = 1.0_RP
    end where
  case default
    LOG_ERROR('fieldutil_get_profile2d',*) trim(profile_name)//' is not supported. Check!'
    call PRC_abort
  end select

  return
end subroutine fieldutil_get_profile2d_tracer

!-------------------------------

subroutine fieldutil_get_profile2d_flow( flow_x, flow_y, &
  profile_name, x, y, params, N)

  implicit none
  !--------------------------------
  integer, intent(in) :: N
  real(RP), intent(out) :: flow_x(N)
  real(RP), intent(out) :: flow_y(N)
  character(*), intent(in) :: profile_name
  real(RP), intent(in) :: x(N)
  real(RP), intent(in) :: y(N)
  real(RP), intent(in) :: params(:)

  real(RP) :: dist(N)
  real(RP) :: fac
  !------------------------------------------------------------------------

  flow_x(:) = 0.0_RP
  flow_y(:) = 0.0_RP

  select case(profile_name)
  case ('constant')
    flow_x(:) = params(1)
    flow_y(:) = params(2)
  case ('rigid-body-rot')
    flow_x(:) = - 2.0_RP*PI/params(3)*(y(:) - params(1))
    flow_y(:) = + 2.0_RP*PI/params(3)*(x(:) - params(2))
  case ('swirling')
    fac = cos(PI*params(4)/params(3))
    flow_x(:) = + sin(PI*x(:))**2 * sin(2.0_RP*PI*y(:)) * fac
    flow_y(:) = - sin(PI*y(:))**2 * sin(2.0_RP*PI*x(:)) * fac
  case default
    LOG_ERROR('fieldutil_get_flow2d',*) trim(profile_name)//' is not supported. Check!'
    call PRC_abort
  end select

  return
end subroutine fieldutil_get_profile2d_flow

end module mod_fieldutil
