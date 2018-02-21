
! ==========================================================================================

module barotropic_diagnostics_testing_mod

use              fms_mod, only: write_version_number, &
                                mpp_pe,      &
                                mpp_root_pe, &
                                error_mesg, FATAL

use       transforms_mod, only: get_grid_boundaries, &
                                get_deg_lon,         &
                                get_deg_lat,         &
                                get_grid_domain,     &
                                get_spec_domain,     & 
                                grid_domain

use     diag_manager_mod, only: diag_axis_init, &
                                register_diag_field, &
                                register_static_field, &
                                send_data

use     time_manager_mod, only: time_type, &
                                get_time

use barotropic_physics_mod,  only:  phys_type
use barotropic_dynamics_mod, only:  testing_type

implicit none
private

public :: barotropic_diagnostics_init, &
          barotropic_diagnostics

character(len=84), parameter :: version = '$Id: barotropic_diagnostics.f90,v 10.0 2003/10/24 22:00:57 fms Exp $'
character(len=84), parameter :: tagname = '$Name: latest $'

logical :: module_is_initialized = .false.

character(len=8)  :: axiset   = 'barotropic'
character(len=84) :: mod_name = 'barotropic_diagnostics'

integer :: id_f_test
integer :: is, ie, js, je

contains

!-----------------------------------------------------------------------------------------------------------------
subroutine barotropic_diagnostics_testing_init(Time, lon_max, lat_max, id_lon, id_lat, id_lonb, id_latb)

type(time_type), intent(in) :: Time
!integer, intent(in) :: lon_max, lat_max
integer, intent(in) :: lon_max, lat_max, id_lon, id_lat, id_lonb, id_latb

real, dimension(lon_max  ) :: lon
real, dimension(lon_max+1) :: lonb
real, dimension(lat_max  ) :: lat
real, dimension(lat_max+1) :: latb

integer, dimension(3) :: axis_3d

integer :: log_unit
integer :: namelist_unit, ierr, io
real    :: rad_to_deg
logical :: used

 call write_version_number(version, tagname)

 call get_grid_domain(is, ie, js, je)

rad_to_deg = 45./atan(1.)
 call get_grid_boundaries(lonb,latb,global=.true.)
 call get_deg_lon(lon)
 call get_deg_lat(lat)

!id_lonb=diag_axis_init('lonb', rad_to_deg*lonb, 'degrees_E', 'x', 'longitude edges', set_name=axiset, Domain2=grid_domain)
!id_latb=diag_axis_init('latb', rad_to_deg*latb, 'degrees_N', 'y', 'latitude edges',  set_name=axiset, Domain2=grid_domain)
!id_lon =diag_axis_init('lon', lon, 'degrees_E', 'x', 'longitude', set_name=axiset, Domain2=grid_domain, edges=id_lonb)
!id_lat =diag_axis_init('lat', lat, 'degrees_N', 'y', 'latitude',  set_name=axiset, Domain2=grid_domain, edges=id_latb)

axis_3d(1) = id_lon
axis_3d(2) = id_lon
axis_3d(3) = id_lat 

! BUDGET TERMS
id_f_test	    		= register_diag_field(mod_name, 'f_test' , 			axis_3d, Time, 'f_test'              			, 'm/s'      ) 

module_is_initialized = .true.

return
end subroutine barotropic_diagnostics_testing_init

!--------------------------------------------------------------------------------------------

subroutine barotropic_diagnostics_testing(Time, Testing, Phys, time_index)

type(time_type), intent(in) :: Time
type(phys_type), intent(in) :: Phys
type(testing_type), intent(in) :: Testing
integer,         intent(in) :: time_index

logical :: used

if(.not.module_is_initialized) call error_mesg('barotropic_diagnostics_testing', &
                              'barotropic_diagnostics_testing_init has not been called', FATAL)

if(id_f_array      > 0) used = send_data(id_f_array      , Testing%f_array     (:,:,:,time_index) , time)

return
end subroutine barotropic_diagnostics_testing

!--------------------------------------------------------------------------------------------

subroutine barotropic_diagnostics_testing_end(Time)

type(time_type), intent(in) :: Time

if(.not.module_is_initialized) call error_mesg('barotropic_diagnostics_end', &
                              'barotropic_diagnostics_init has not been called', FATAL)

module_is_initialized = .false.

return
end subroutine barotropic_diagnostics_testing_end
!--------------------------------------------------------------------------------------------

end module barotropic_diagnostics_testing_mod
