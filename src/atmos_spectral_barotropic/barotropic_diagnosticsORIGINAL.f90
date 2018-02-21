
! ==========================================================================================

module barotropic_diagnostics_mod

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
use barotropic_dynamics_mod, only:  grid_type

implicit none
private

public :: barotropic_diagnostics_init, &
          barotropic_diagnostics

character(len=84), parameter :: version = '$Id: barotropic_diagnostics.f90,v 10.0 2003/10/24 22:00:57 fms Exp $'
character(len=84), parameter :: tagname = '$Name: latest $'

logical :: module_is_initialized = .false.

character(len=8)  :: axiset   = 'barotropic'
character(len=84) :: mod_name = 'barotropic_diagnostics'

integer :: id_vor, id_stream, id_pv, id_u, id_utend, id_vtend, id_sink, id_sink_beta, id_dens_dt, id_robertssink, id_source, id_source_beta, id_ensflux, id_ensflux_div, id_ensflux_div_beta, id_vortend, id_vvor_beta, id_v, id_trs, id_tr
integer :: is, ie, js, je

contains

!-----------------------------------------------------------------------------------------------------------------
!subroutine barotropic_diagnostics_init(Time, lon_max, lat_max)
subroutine barotropic_diagnostics_init(Time, lon_max, lat_max, id_lon, id_lat, id_lonb, id_latb)

type(time_type), intent(in) :: Time
!integer, intent(in) :: lon_max, lat_max
integer, intent(in) :: lon_max, lat_max, id_lon, id_lat, id_lonb, id_latb

real, dimension(lon_max  ) :: lon
real, dimension(lon_max+1) :: lonb
real, dimension(lat_max  ) :: lat
real, dimension(lat_max+1) :: latb

integer, dimension(2) :: axis_2d

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

axis_2d(1) = id_lon
axis_2d(2) = id_lat 

id_sink              = register_diag_field(mod_name, 'sink' ,             axis_2d, Time, 'sink'                , 'm/s'      ) 
id_sink_beta         = register_diag_field(mod_name, 'sink_beta' ,        axis_2d, Time, 'sink_beta'           , 'm/s'      ) 
id_robertssink       = register_diag_field(mod_name, 'robertssink' ,      axis_2d, Time, 'robertssink'         , 'm/s'      ) 
id_source            = register_diag_field(mod_name, 'source' ,           axis_2d, Time, 'source'              , 'm/s'      ) 
id_source_beta       = register_diag_field(mod_name, 'source_beta' ,      axis_2d, Time, 'source_beta'         , 'm/s'      ) 
id_ensflux           = register_diag_field(mod_name, 'ensflux' ,          axis_2d, Time, 'ensflux'             , 'm/s'      ) 
id_ensflux_div       = register_diag_field(mod_name, 'ensflux_div' ,      axis_2d, Time, 'ensflux_div'         , 'm/s'      ) 
id_ensflux_div_beta  = register_diag_field(mod_name, 'ensflux_div_beta' , axis_2d, Time, 'ensflux_div_beta'    , 'm/s'      ) 
id_dens_dt           = register_diag_field(mod_name, 'dens_dt' ,          axis_2d, Time, 'dens_dt'             , 'm/s'      ) 
id_utend             = register_diag_field(mod_name, 'ucomp_tend' ,       axis_2d, Time, 'u_wind_tend'         , 'm/s'      ) 
id_vtend             = register_diag_field(mod_name, 'vcomp_tend' ,       axis_2d, Time, 'v_wind_tend'         , 'm/s'      ) 
id_vortend           = register_diag_field(mod_name, 'vor_tend' ,         axis_2d, Time, 'vor_tend'            , 'm/s'      ) 
id_vvor_beta         = register_diag_field(mod_name, 'vvor_beta' ,        axis_2d, Time, 'vvor_beta'           , 'm/s'      ) 
id_u                 = register_diag_field(mod_name, 'ucomp' ,            axis_2d, Time, 'u_wind'              , 'm/s'      ) 
id_v                 = register_diag_field(mod_name, 'vcomp' ,            axis_2d, Time, 'v_wind'              , 'm/s'      ) 
id_vor               = register_diag_field(mod_name, 'vor'   ,            axis_2d, Time, 'relative vorticity'  , '1/s'      )
id_pv                = register_diag_field(mod_name, 'pv'    ,            axis_2d, Time, 'absolute vorticity'  , '1/s'      )
id_stream            = register_diag_field(mod_name, 'stream',            axis_2d, Time, 'streamfunction'      , 'm^2/s'    )
id_trs               = register_diag_field(mod_name, 'trs'   ,            axis_2d, Time, 'spectral tracer'     , 'none'     )
id_tr                = register_diag_field(mod_name, 'tr'    ,            axis_2d, Time, 'grid tracer'         , 'none'     )

module_is_initialized = .true.

return
end subroutine barotropic_diagnostics_init

!--------------------------------------------------------------------------------------------

subroutine barotropic_diagnostics(Time, Grid, Phys, time_index)

type(time_type), intent(in) :: Time
type(phys_type), intent(in) :: Phys
type(grid_type), intent(in) :: Grid
integer,         intent(in) :: time_index

logical :: used

if(.not.module_is_initialized) call error_mesg('barotropic_diagnostics', &
                              'barotropic_diagnostics_init has not been called', FATAL)
   
if(id_sink                  > 0) used = send_data(id_sink      , Grid%sink     (:,:,time_index) , time)
if(id_sink_beta             > 0) used = send_data(id_sink_beta      , Grid%sink_beta     (:,:,time_index) , time)
if(id_robertssink           > 0) used = send_data(id_robertssink      , Grid%robertssink     (:,:,time_index) , time)
if(id_source                > 0) used = send_data(id_source      , Grid%source     (:,:,time_index) , time)
if(id_source_beta           > 0) used = send_data(id_source_beta      , Grid%source_beta     (:,:,time_index) , time)
if(id_ensflux               > 0) used = send_data(id_ensflux      , Grid%ensflux     (:,:,time_index) , time)
if(id_ensflux_div           > 0) used = send_data(id_ensflux_div      , Grid%ensflux_div     (:,:,time_index) , time)
if(id_ensflux_div_beta      > 0) used = send_data(id_ensflux_div_beta      , Grid%ensflux_div_beta     (:,:,time_index) , time)
if(id_dens_dt               > 0) used = send_data(id_dens_dt      , Grid%dens_dt     (:,:,time_index) , time)
if(id_utend                 > 0) used = send_data(id_utend      , Grid%utend     (:,:,time_index) , time)
if(id_vtend                 > 0) used = send_data(id_vtend      , Grid%vtend     (:,:,time_index) , time)
if(id_vortend      > 0) used = send_data(id_vortend      , Grid%vortend     (:,:,time_index) , time)
if(id_vvor_beta      > 0) used = send_data(id_vvor_beta      , Grid%vvor_beta     (:,:,time_index) , time)
if(id_u      > 0) used = send_data(id_u      , Grid%u     (:,:,time_index) , time)
if(id_v      > 0) used = send_data(id_v      , Grid%v     (:,:,time_index) , time)
if(id_vor    > 0) used = send_data(id_vor    , Grid%vor   (:,:,time_index) , time)
if(id_pv     > 0) used = send_data(id_pv     , Grid%pv    (:,:)            , time)
if(id_stream > 0) used = send_data(id_stream , Grid%stream(:,:)            , time)
if(id_tr     > 0) used = send_data(id_tr     , Grid%tr    (:,:,time_index) , time)
if(id_trs    > 0) used = send_data(id_trs    , Grid%trs   (:,:,time_index) , time)

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
