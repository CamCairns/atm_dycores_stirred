
! ==========================================================================================

module barotropic_diagnostics_mod

use              fms_mod, only: write_version_number, &
                                mpp_pe,      &
                                mpp_root_pe, &
                                error_mesg, FATAL

use       transforms_mod, only: get_grid_domain,     &
                                get_spec_domain,     & 
                                grid_domain,         &
                                trans_spherical_to_grid

use     diag_manager_mod, only: register_diag_field, &
                                register_static_field, &
                                send_data

use     time_manager_mod, only: time_type, &
                                get_time

use barotropic_physics_mod,  only:  phys_type
use barotropic_dynamics_mod, only:  grid_type, dynamics_type

implicit none
private

public :: barotropic_diagnostics_init, &
          barotropic_diagnostics

character(len=84), parameter :: version = '$Id: barotropic_diagnostics.F90,v 17.0 2009/07/21 03:00:18 fms Exp $'
character(len=84), parameter :: tagname = '$Name: riga $'

logical :: module_is_initialized = .false.

character(len=8)  :: axiset   = 'barotropic'
character(len=84) :: mod_name = 'barotropic_diagnostics'

integer :: id_vor, id_stream, id_pv, id_u, id_v, id_trs, id_tr, id_eddy_vor, id_delta_u
integer :: is, ie, js, je, ms, me, ns, ne

contains

!-----------------------------------------------------------------------------------------------------------------
subroutine barotropic_diagnostics_init(Time, lon_max, lat_max, id_lon, id_lat, id_lonb, id_latb)

type(time_type), intent(in) :: Time
integer, intent(in) :: lon_max, lat_max, id_lon, id_lat, id_lonb, id_latb

integer, dimension(2) :: axis_2d

integer :: log_unit
integer :: namelist_unit, ierr, io
logical :: used

call write_version_number(version, tagname)

call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)

axis_2d(1) = id_lon
axis_2d(2) = id_lat 

id_u       = register_diag_field(mod_name, 'ucomp'   , axis_2d, Time, 'u_wind'            , 'm/s'  ) 
id_v       = register_diag_field(mod_name, 'vcomp'   , axis_2d, Time, 'v_wind'            , 'm/s'  ) 
id_vor     = register_diag_field(mod_name, 'vor'     , axis_2d, Time, 'relative vorticity', '1/s'  )
id_pv      = register_diag_field(mod_name, 'pv'      , axis_2d, Time, 'absolute vorticity', '1/s'  )
id_stream  = register_diag_field(mod_name, 'stream'  , axis_2d, Time, 'streamfunction'    , 'm^2/s')
id_trs     = register_diag_field(mod_name, 'trs'     , axis_2d, Time, 'spectral tracer'   , 'none' )
id_tr      = register_diag_field(mod_name, 'tr'      , axis_2d, Time, 'grid tracer'       , 'none' )
id_eddy_vor= register_diag_field(mod_name, 'eddy_vor', axis_2d, Time, 'eddy vorticity'    , '1/s'  )
id_delta_u = register_diag_field(mod_name, 'delta_u' , axis_2d, Time, 'change in zonal wind','m/s' )

module_is_initialized = .true.

return
end subroutine barotropic_diagnostics_init

!--------------------------------------------------------------------------------------------

subroutine barotropic_diagnostics(Time, Dyn, Phys, time_index)

type(time_type), intent(in) :: Time
type(phys_type), intent(in) :: Phys
type(dynamics_type), intent(in) :: Dyn
integer,         intent(in) :: time_index

real,    dimension(is:ie,js:je) :: grid_tmp
real,    dimension(is:ie,js:je) :: delta_zonal_u
complex, dimension(ms:me,ns:ne) :: spec_tmp
logical :: used
integer :: j

if(.not.module_is_initialized) call error_mesg('barotropic_diagnostics', &
                              'barotropic_diagnostics_init has not been called', FATAL)
   
if(id_u        > 0) used = send_data(id_u       , Dyn%Grid%u       (:,:,time_index) , time)
if(id_v        > 0) used = send_data(id_v       , Dyn%Grid%v       (:,:,time_index) , time)
if(id_vor      > 0) used = send_data(id_vor     , Dyn%Grid%vor     (:,:,time_index) , time)
if(id_pv       > 0) used = send_data(id_pv      , Dyn%Grid%pv      (:,:)            , time)
if(id_stream   > 0) used = send_data(id_stream  , Dyn%Grid%stream  (:,:)            , time)
if(id_tr       > 0) used = send_data(id_tr      , Dyn%Grid%tr      (:,:,time_index) , time)
if(id_trs      > 0) used = send_data(id_trs     , Dyn%Grid%trs     (:,:,time_index) , time)
if(id_eddy_vor > 0) then
  spec_tmp = Dyn%Spec%vor(:,:,time_index)
  if(ms == 0) spec_tmp(0,:) = cmplx(0.0,0.0)
  call trans_spherical_to_grid(spec_tmp, grid_tmp)
  used = send_data(id_eddy_vor, grid_tmp, time)
endif
if(id_delta_u > 0) then
  do j=js,je
    delta_zonal_u(:,j) = sum(Dyn%Grid%u(:,j,time_index))/Dyn%num_lon - Dyn%Grid%zonal_u_init(j)
  enddo
  used = send_data(id_delta_u, delta_zonal_u, time)
endif

return
end subroutine barotropic_diagnostics
!--------------------------------------------------------------------------------------------

subroutine barotropic_diagnostics_end(Time)

type(time_type), intent(in) :: Time

if(.not.module_is_initialized) call error_mesg('barotropic_diagnostics_end', &
                              'barotropic_diagnostics_init has not been called', FATAL)

module_is_initialized = .false.

return
end subroutine barotropic_diagnostics_end
!--------------------------------------------------------------------------------------------

end module barotropic_diagnostics_mod
