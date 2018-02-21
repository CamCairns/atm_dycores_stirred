
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
!use barotropic_dynamics_mod, only:  testing_type

implicit none
private

public :: barotropic_diagnostics_init, &
          barotropic_diagnostics

 character(len=84), parameter :: version = '$Id: barotropic_diagnostics.f90,v 10.0 2003/10/24 22:00:57 fms Exp $'
 character(len=84), parameter :: tagname = '$Name: latest $'

logical :: module_is_initialized = .false.

 character(len=8)  :: axiset   = 'barotropic'
 character(len=84) :: mod_name = 'barotropic_diagnostics'

integer :: id_vor, id_stream, id_pv, id_u, id_v, id_utend, id_vtend, id_vortend, id_utrunc, id_vtrunc, id_vortrunc, id_ulindamp_m, id_ulindamp_e, id_vlindamp, id_vorlindamp_m, id_vorlindamp_e, id_ustir, id_vstir, id_vorstir, id_uspecdamp, id_vspecdamp, id_vorspecdamp, id_urobdamp, id_vrobdamp, id_vq, id_uq, id_dt_vor, id_vorrobdamp, id_trs, id_tr, id_dx_u2_v2, id_dy_u2_v2, id_pvor_advec, id_rvor_advec, id_rvormean_advec, id_rvorprime_advec, id_beta, id_beta_star, id_vormean_x, id_vormean_y, id_vorprime_x, id_vorprime_y, id_vormean_y_tend, id_vorprime, id_vorprime_prev, id_vorprime_curr, id_vork, id_vork_prev, id_vork_curr, id_vortendk, id_vor_x, id_vor_y, id_vormeantend, id_vorprimetend, id_enstrophy_rvormean_adevc_Matlab, id_enstrophy_rvorprime_adevc_Matlab
integer :: id_uprev, id_ucurr, id_vprev, id_vcurr,  id_vorprev, id_vorcurr, id_v_y, id_rvor_adveck, id_u_dvordx, id_v_dvordy, id_vor_dudx, id_vor_dvdy
integer :: id_energy, id_energy_tend, id_energy_voradvec, id_energy_gradterm, id_energy_trunc, id_energy_lindamp_m, id_energy_lindamp_e, id_energy_specdamp, id_energy_stir, id_energy_robdamp, id_energy_mean, id_energy_tend_mean, id_energy_voradvec_mean, id_energy_gradterm_mean, id_energy_trunc_mean, id_energy_lindamp_m_mean, id_energy_lindamp_e_mean, id_energy_specdamp_mean, id_energy_stir_mean, id_energy_robdamp_mean, id_tester, id_tester2, id_tester3, id_tester4, id_tester5, id_tester6, id_vorprime_yk, id_v_yk, id_vorprimek_curr, id_vk
integer :: id_enstrophy, id_enstrophy_tend, id_enstrophy_dt_vor, id_enstrophy_pvor_advec, id_enstrophy_rvormean_advec, id_enstrophy_rvorprime_advec, id_enstrophy_trunc, id_enstrophy_lindamp_m, id_enstrophy_lindamp_e, id_enstrophy_specdamp, id_enstrophy_stir, id_enstrophy_robdamp
integer :: id_enstrophy_tendk, id_enstrophy_pvor_adveck, id_enstrophy_rvormean_adveck, id_enstrophy_trunck, id_enstrophy_lindamp_ek, id_enstrophy_lindamp_mk, id_enstrophy_specdampk, id_enstrophy_robdampk, id_enstrophy_stirk, id_enstrophy_rvorprime_adveck_SUM, id_vorprime_fluxk, id_I_higher_freq, id_I1(1:12,1:12), id_I2(1:12,1:12), id_I3(1:12,1:12), id_I4(1:12,1:12), id_spec_deriv_vq_real(1:128,1), id_spec_deriv_vq_imag(1:128,1)
integer :: id_I_pp, id_I_sp, id_I_ps, id_I_hp, id_I_ph, id_I_ss, id_I_hh, id_I_sh, id_I_hs, id_I_0p, id_I_p0, id_I_0s, id_I_s0, id_I_0h, id_I_h0, id_I_00, id_I_XX
!integer :: id_I2_pp, id_I2_sp, id_I2_ps, id_I2_hp, id_I2_ph, id_I2_ss, id_I2_hh, id_I2_sh, id_I2_hs
!integer :: id_sink, id_sink_beta, id_dens_dt, id_robertssink, id_source, id_source_beta, id_ensflux, id_ensflux_div, id_ensflux_div_beta, id_vvor_beta
integer :: is, ie, js, je, p, q

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
 integer :: p, q
 character(LEN=20) :: output
! character(LEN=3) :: output2
! character(LEN=3) :: output3
! character(LEN=3) :: output4
  character(LEN=2) :: pstr
 character(LEN=2) :: qstr
 character(LEN=3) :: rstr
 character(LEN=8) :: filestr
 character(LEN=7) :: filestr2
 character(LEN=25) :: filestr3
 character(LEN=24) :: filestr4
! character(LEN=8) :: filestrII
! character(LEN=7) :: filestrII2
! character(LEN=8) :: filestrIII
! character(LEN=7) :: filestrIII2
! character(LEN=8) :: filestrIV
! character(LEN=7) :: filestrIV2

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

! BUDGET TERMS
id_u		    			= register_diag_field(mod_name, 'ucomp' , 		axis_2d, Time, 'u_wind'              			, 'm/s'      ) 
id_v     				= register_diag_field(mod_name, 'vcomp' , 		axis_2d, Time, 'v_wind'              			, 'm/s'      ) 
id_uprev				= register_diag_field(mod_name, 'uprev' , 		axis_2d, Time, 'u wind tendency'              		, 'm/s^2'    ) 
id_ucurr				= register_diag_field(mod_name, 'ucurr' , 		axis_2d, Time, 'u wind tendency'              		, 'm/s^2'    )
id_vprev				= register_diag_field(mod_name, 'vprev' , 		axis_2d, Time, 'u wind tendency'              		, 'm/s^2'    )  
id_vcurr				= register_diag_field(mod_name, 'vcurr' ,		axis_2d, Time, 'v wind tendency'              		, 'm/s^2'    ) 
id_vorprev				= register_diag_field(mod_name, 'vorprev' , 		axis_2d, Time, 'u wind tendency'              		, 'm/s^2'    )  
id_vorcurr				= register_diag_field(mod_name, 'vorcurr' ,		axis_2d, Time, 'v wind tendency'              		, 'm/s^2'    ) 
id_utend				= register_diag_field(mod_name, 'utend' , 		axis_2d, Time, 'u wind tendency'              		, 'm/s^2'    ) 
id_vtend				= register_diag_field(mod_name, 'vtend' ,		axis_2d, Time, 'v wind tendency'              		, 'm/s^2'    ) 
id_vortend				= register_diag_field(mod_name, 'vortend' , 		axis_2d, Time, 'vorticity tendency'           		, '1/s^2'    )
id_dt_vor				= register_diag_field(mod_name, 'dt_vor' , 		axis_2d, Time, 'vorticity forcing'        		, '1/s^2'    )  
id_utrunc				= register_diag_field(mod_name, 'utrunc', 		axis_2d, Time, 'u truncation'        			, 'm/s^2'    ) 
id_vtrunc				= register_diag_field(mod_name, 'vtrunc', 		axis_2d, Time, 'v truncation'        			, 'm/s^2'    )
id_vortrunc				= register_diag_field(mod_name, 'vortrunc', 		axis_2d, Time, 'vorticity truncation'        		, '1/s^2'    )
id_ulindamp_m				= register_diag_field(mod_name, 'ulindamp_m' , 		axis_2d, Time, 'u mean flow linear damping'   		, 'm/s^2'    ) 
id_ulindamp_e				= register_diag_field(mod_name, 'ulindamp_e' , 		axis_2d, Time, 'u eddy linear damping'        		, 'm/s^2'    )
id_vlindamp				= register_diag_field(mod_name, 'vlindamp' , 		axis_2d, Time, 'v (eddy) linear damping'      		, 'm/s^2'    )
id_vorlindamp_m				= register_diag_field(mod_name, 'vorlindamp_m' , 	axis_2d, Time, 'vorticity mean flow linear damping' 	, '1/s^2'    ) 
id_vorlindamp_e				= register_diag_field(mod_name, 'vorlindamp_e' , 	axis_2d, Time, 'vorticity eddy linear damping' 		, '1/s^2'    )
id_dx_u2_v2				= register_diag_field(mod_name, 'dx_u2_v2' ,	 	axis_2d, Time, 'dx_(u^2+v^2)'		 		, 'm/s^2'    )
id_dy_u2_v2				= register_diag_field(mod_name, 'dy_u2_v2' ,	 	axis_2d, Time, 'dy_(u^2+v^2)'		 		, 'm/s^2'    )
id_ustir				= register_diag_field(mod_name, 'ustir' , 		axis_2d, Time, 'u tendency due to stirring'   		, 'm/s^2'    ) 
id_vstir				= register_diag_field(mod_name, 'vstir' , 		axis_2d, Time, 'v tendency due to stirring'        	, 'm/s^2'    )
id_vorstir				= register_diag_field(mod_name, 'vorstir' , 		axis_2d, Time, 'vorticity tendency due to stirring'	, '1/s^2'    )
id_uspecdamp				= register_diag_field(mod_name, 'uspecdamp' , 		axis_2d, Time, 'u spectral damping'          		, 'm/s^2'    ) 
id_vspecdamp				= register_diag_field(mod_name, 'vspecdamp' , 		axis_2d, Time, 'v spectral damping'           		, 'm/s^2'    )
id_vorspecdamp				= register_diag_field(mod_name, 'vorspecdamp' , 	axis_2d, Time, 'vorticity spectral damping'        	, '1/s^2'    ) 
id_urobdamp				= register_diag_field(mod_name, 'urobdamp' , 		axis_2d, Time, 'u roberts damping'          		, 'm/s^2'    ) 
id_vrobdamp				= register_diag_field(mod_name, 'vrobdamp' , 		axis_2d, Time, 'v roberts damping'           		, 'm/s^2'    )
id_vorrobdamp				= register_diag_field(mod_name, 'vorrobdamp' , 		axis_2d, Time, 'vorticity roberts damping'           	, '1/s^2'    ) 
id_pvor_advec				= register_diag_field(mod_name, 'pvor_advec' , 		axis_2d, Time, 'planetary vorticity advection' 		, '1/s^2'    )
id_rvor_advec				= register_diag_field(mod_name, 'rvor_advec' , 		axis_2d, Time, 'relative vorticity advection'          	, '1/s^2'    )
id_rvormean_advec			= register_diag_field(mod_name, 'rvormean_advec' ,	axis_2d, Time, 'relative mean vorticity advection'    	, '1/s^2'    )
id_rvorprime_advec			= register_diag_field(mod_name, 'rvorprime_advec' ,	axis_2d, Time, 'relative prime vorticity advection'    	, '1/s^2'    )
id_beta					= register_diag_field(mod_name, 'beta' , 		axis_2d, Time, 'beta = df/dy' 		         	, '1/(s^2*m)')
id_beta_star				= register_diag_field(mod_name, 'beta_star' , 		axis_2d, Time, 'beta_star = beta - u_{yy}'         	, '1/(s^2*m)')
id_vormean_x				= register_diag_field(mod_name, 'vormean_x' , 		axis_2d, Time, 'meridional derivative of zonal mean vor', '1/s^2*m'  )
id_vormean_y				= register_diag_field(mod_name, 'vormean_y' , 		axis_2d, Time, 'meridional derivative of zonal mean vor', '1/s^2*m'  )
id_vorprime_x				= register_diag_field(mod_name, 'vorprime_x' , 		axis_2d, Time, 'meridional derivative of zonal mean vor', '1/s^2*m'  )
id_vorprime_y				= register_diag_field(mod_name, 'vorprime_y' , 		axis_2d, Time, 'meridional derivative of zonal mean vor', '1/s^2*m'  )
id_vormean_y_tend			= register_diag_field(mod_name, 'vormean_y_tend' ,	axis_2d, Time, 'tend of y derivative of zonal mean vor' ,'1/(s^2*m*t)')
id_vor_x				= register_diag_field(mod_name, 'vor_x' ,		axis_2d, Time, 'tend of y derivative of zonal mean vor' ,'1/(s^2*m*t)')
id_vor_y				= register_diag_field(mod_name, 'vor_y' ,		axis_2d, Time, 'tend of y derivative of zonal mean vor' ,'1/(s^2*m*t)')
id_vorprime				= register_diag_field(mod_name, 'vorprime' ,		axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vorprime_prev			= register_diag_field(mod_name, 'vorprime_prev' ,	axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vorprime_curr			= register_diag_field(mod_name, 'vorprime_curr' ,	axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vork					= register_diag_field(mod_name, 'vork' ,		axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vork_prev				= register_diag_field(mod_name, 'vork_prev' ,		axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vork_curr				= register_diag_field(mod_name, 'vork_curr' ,		axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vortendk 				= register_diag_field(mod_name, 'vortendk' ,		axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vormeantend				= register_diag_field(mod_name, 'vormeantend' ,		axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vorprimetend				= register_diag_field(mod_name, 'vorprimetend' ,	axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_rvormean_adevc_Matlab	= register_diag_field(mod_name, 'enstrophy_rvormean_adevc_Matlab' ,	axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_rvorprime_adevc_Matlab	= register_diag_field(mod_name, 'enstrophy_rvorprime_adevc_Matlab' ,	axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_tester				= register_diag_field(mod_name, 'tester' , 		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_tester2				= register_diag_field(mod_name, 'tester2' , 		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_tester3				= register_diag_field(mod_name, 'tester3' , 		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_tester4				= register_diag_field(mod_name, 'tester4' , 		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_tester5				= register_diag_field(mod_name, 'tester5' , 		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_tester6				= register_diag_field(mod_name, 'tester6' , 		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_v_y					= register_diag_field(mod_name, 'v_y' , 		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_rvor_adveck				= register_diag_field(mod_name, 'rvor_adveck' ,		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_u_dvordx				= register_diag_field(mod_name, 'u_dvordx' ,		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_v_dvordy				= register_diag_field(mod_name, 'v_dvordy' ,		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_vor_dudx				= register_diag_field(mod_name, 'vor_dudx' ,		axis_2d, Time, 'Test Variable'		         	, '???'      )
id_vor_dvdy				= register_diag_field(mod_name, 'vor_dvdy' ,		axis_2d, Time, 'Test Variable'		         	, '???'      )

! ENERGY TERMS	
id_energy				= register_diag_field(mod_name, 'energy' , 			axis_2d, Time, 'energy'			           	, 'm^2/s^2'  )
id_energy_tend				= register_diag_field(mod_name, 'energy_tend' ,			axis_2d, Time, 'energy tendency'	           	, 'm^2/s^3'  )
id_energy_voradvec			= register_diag_field(mod_name, 'energy_voradvec' ,		axis_2d, Time, 'energy vorticity advection term'      	, 'm^2/s^3'  )
id_energy_gradterm			= register_diag_field(mod_name, 'energy_gradterm' ,		axis_2d, Time, 'energy gradient term term'      	, 'm^2/s^3'  )        
id_energy_trunc				= register_diag_field(mod_name, 'energy_trunc' ,		axis_2d, Time, 'energy truncation term'		      	, 'm^2/s^3'  )    
id_energy_lindamp_m			= register_diag_field(mod_name, 'energy_lindamp_m' ,		axis_2d, Time, 'energy mean flow linear damping term' 	, 'm^2/s^3'  )
id_energy_lindamp_e			= register_diag_field(mod_name, 'energy_lindamp_e' ,		axis_2d, Time, 'energy eddy linear damping term'      	, 'm^2/s^3'  )
id_energy_specdamp			= register_diag_field(mod_name, 'energy_specdamp' ,		axis_2d, Time, 'energy spectral damping term'      	, 'm^2/s^3'  )
id_energy_stir				= register_diag_field(mod_name, 'energy_stir' ,			axis_2d, Time, 'energy stirring term'		      	, 'm^2/s^3'  )
id_energy_robdamp			= register_diag_field(mod_name, 'energy_robdamp' ,		axis_2d, Time, 'energy Roberts damping term'	      	, 'm^2/s^3'  )
! MEAN ENERGY TERMS
id_energy_mean				= register_diag_field(mod_name, 'energy_mean' , 		axis_2d, Time, 'zonal mean energy'			        , 'm^2/s^2'  )
id_energy_tend_mean			= register_diag_field(mod_name, 'energy_tend_mean' ,		axis_2d, Time, 'zonal mean energy tendency'	           	, 'm^2/s^3'  )
id_energy_voradvec_mean			= register_diag_field(mod_name, 'energy_voradvec_mean' ,	axis_2d, Time, 'zonal mean energy vorticity advection term'     , 'm^2/s^3'  )
id_energy_gradterm_mean			= register_diag_field(mod_name, 'energy_gradterm_mean' ,	axis_2d, Time, 'zonal mean energy gradient term term'      	, 'm^2/s^3'  )        
id_energy_trunc_mean			= register_diag_field(mod_name, 'energy_trunc_mean' ,		axis_2d, Time, 'zonal mean energy truncation term'		, 'm^2/s^3'  )    
id_energy_lindamp_m_mean		= register_diag_field(mod_name, 'energy_lindamp_m_mean' ,	axis_2d, Time, 'zonal mean energy mean flow linear damping term', 'm^2/s^3'  )
id_energy_lindamp_e_mean		= register_diag_field(mod_name, 'energy_lindamp_e_mean' ,	axis_2d, Time, 'zonal mean energy eddy linear damping term'     , 'm^2/s^3'  )
id_energy_specdamp_mean			= register_diag_field(mod_name, 'energy_specdamp_mean' ,	axis_2d, Time, 'zonal mean energy spectral damping term'      	, 'm^2/s^3'  )
id_energy_stir_mean			= register_diag_field(mod_name, 'energy_stir_mean' ,		axis_2d, Time, 'zonal mean energy stirring term'		,'m^2/s^3'   )
id_energy_robdamp_mean			= register_diag_field(mod_name, 'energy_robdamp_mean' ,		axis_2d, Time, 'zonal mean energy Roberts damping term'	      	, 'm^2/s^3'  )
! ENSTROPHY TERMS
id_enstrophy				= register_diag_field(mod_name, 'enstrophy' , 			axis_2d, Time, 'enstrophy'			        , 'm^2/s^2'  )
id_enstrophy_tend			= register_diag_field(mod_name, 'enstrophy_tend' ,		axis_2d, Time, 'enstrophy tendency'	           	, 'm^2/s^3'  )
id_enstrophy_dt_vor			= register_diag_field(mod_name, 'enstrophy_dt_vor' , 		axis_2d, Time, 'enstrophy absolute vorticity advection'	, 'm^2/s^3'    )
id_enstrophy_pvor_advec			= register_diag_field(mod_name, 'enstrophy_pvor_advec' ,	axis_2d, Time, 'enstrophy planetary vorticity advection', 'm^2/s^3'  )
id_enstrophy_rvormean_advec		= register_diag_field(mod_name, 'enstrophy_rvormean_advec' ,	axis_2d, Time, 'enstrophy relative mean vorticity advection' , 'm^2/s^3'  )
id_enstrophy_rvorprime_advec		= register_diag_field(mod_name, 'enstrophy_rvorprime_advec' ,	axis_2d, Time, 'enstrophy relative prime vorticity advection' , 'm^2/s^3'  )                
id_enstrophy_trunc			= register_diag_field(mod_name, 'enstrophy_trunc' ,		axis_2d, Time, 'enstrophy truncation term'	      	, 'm^2/s^3'  )    
id_enstrophy_lindamp_m			= register_diag_field(mod_name, 'enstrophy_lindamp_m' ,		axis_2d, Time, 'enstrophy mean flow linear damping' 	, 'm^2/s^3'  )
id_enstrophy_lindamp_e			= register_diag_field(mod_name, 'enstrophy_lindamp_e' ,		axis_2d, Time, 'enstrophy eddy linear damping term'    	, 'm^2/s^3'  )
id_enstrophy_specdamp			= register_diag_field(mod_name, 'enstrophy_specdamp' ,		axis_2d, Time, 'enstrophy spectral damping term'      	, 'm^2/s^3'  )
id_enstrophy_stir			= register_diag_field(mod_name, 'enstrophy_stir' ,		axis_2d, Time, 'enstrophy stirring term'	      	, 'm^2/s^3'  )
id_enstrophy_robdamp			= register_diag_field(mod_name, 'enstrophy_robdamp' ,		axis_2d, Time, 'enstrophy Roberts damping term'	      	, 'm^2/s^3'  )
id_enstrophy_tendk			= register_diag_field(mod_name, 'enstrophy_tendk' ,		axis_2d, Time, 'tend of enstrophy cospec decomp along lon dim'	,'m^2/s^2'   )
id_enstrophy_pvor_adveck		= register_diag_field(mod_name, 'enstrophy_pvor_adveck'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_rvormean_adveck		= register_diag_field(mod_name, 'enstrophy_rvormean_adveck'	,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_trunck			= register_diag_field(mod_name, 'enstrophy_trunck'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_lindamp_ek			= register_diag_field(mod_name, 'enstrophy_lindamp_ek'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_lindamp_mk			= register_diag_field(mod_name, 'enstrophy_lindamp_mk'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_specdampk			= register_diag_field(mod_name, 'enstrophy_specdampk'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_robdampk			= register_diag_field(mod_name, 'enstrophy_robdampk'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_stirk			= register_diag_field(mod_name, 'enstrophy_stirk'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_enstrophy_rvorprime_adveck_SUM	= register_diag_field(mod_name, 'enstrophy_rvorprime_adveck_SUM',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vorprime_yk				= register_diag_field(mod_name, 'vorprime_yk'			,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_v_yk					= register_diag_field(mod_name, 'v_yk'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vk					= register_diag_field(mod_name, 'vk'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_vorprimek_curr			= register_diag_field(mod_name, 'vorprimek_curr'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_higher_freq			= register_diag_field(mod_name, 'I_higher_freq'			,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)

! FLUX TERMS
id_vorprime_fluxk			= register_diag_field(mod_name, 'vorprime_fluxk'		,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_pp					= register_diag_field(mod_name, 'I_pp'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_sp					= register_diag_field(mod_name, 'I_sp'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_ps					= register_diag_field(mod_name, 'I_ps'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_hp					= register_diag_field(mod_name, 'I_hp'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_ph					= register_diag_field(mod_name, 'I_ph'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_ss					= register_diag_field(mod_name, 'I_ss'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_hh					= register_diag_field(mod_name, 'I_hh'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_sh					= register_diag_field(mod_name, 'I_sh'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_hs					= register_diag_field(mod_name, 'I_hs'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_pp				= register_diag_field(mod_name, 'I_pp'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_sp				= register_diag_field(mod_name, 'I_sp'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_ps				= register_diag_field(mod_name, 'I_ps'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_hp				= register_diag_field(mod_name, 'I_hp'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_ph				= register_diag_field(mod_name, 'I_ph'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_ss				= register_diag_field(mod_name, 'I_ss'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_hh				= register_diag_field(mod_name, 'I_hh'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_sh				= register_diag_field(mod_name, 'I_sh'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2_hs				= register_diag_field(mod_name, 'I_hs'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_0p					= register_diag_field(mod_name, 'I_0p'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_p0					= register_diag_field(mod_name, 'I_p0'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_0s					= register_diag_field(mod_name, 'I_0s'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_s0					= register_diag_field(mod_name, 'I_s0'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_0h					= register_diag_field(mod_name, 'I_0h'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_h0					= register_diag_field(mod_name, 'I_h0'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_00					= register_diag_field(mod_name, 'I_00'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
id_I_XX					= register_diag_field(mod_name, 'I_XX'				,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)

!do p=1,12
!   do q = 1,12
!       id_I(p,q)			= register_diag_field(mod_name, 'I',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!   end do
!end do

!id_I					= register_diag_field(mod_name, 'I(1,1)',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I2					= register_diag_field(mod_name, 'I(1,2)',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)

!id_I(1,1)				= register_diag_field(mod_name, 'I(1,1)',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I(1,2)				= register_diag_field(mod_name, 'I(1,2)',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I(1,1)				= register_diag_field(mod_name, 'I_0101',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I(1,1)				= register_diag_field(mod_name, 'I_0102',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I(1,2)				= register_diag_field(mod_name, 'I_12',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!id_I(1,2)				= register_diag_field(mod_name, 'I(', p, ',', q, '2)',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)

!do p=1,1
!	id_spec_deriv_vq_real(p,1)		= register_diag_field(mod_name, 'spec_deriv_vq_real_p001',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!	id_spec_deriv_vq_imag(p,1)		= register_diag_field(mod_name, 'spec_deriv_vq_imag_p001',axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
!end do

output = 'spec_deriv_vq_real_p'

do p=1,128!12
     write(rstr ,"(I3.3)") p
!     print *, rstr
!     write(*,*) trim(output)//pstr//qstr
     write(filestr3,*) adjustl(trim(output))//rstr
!     print *, filestr3(2:24)
     filestr4 = filestr3(2:24)
!     print *, filestr4
     id_spec_deriv_vq_real(p,1)		= register_diag_field(mod_name, filestr4,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
end do

output = 'spec_deriv_vq_imag_p'

do p=1,128!12
     write(rstr ,"(I3.3)") p
!     print *, rstr
!     write(*,*) trim(output)//pstr//qstr
     write(filestr3,*) adjustl(trim(output))//rstr
!     print *, filestr3(2:24)
     filestr4 = filestr3(2:24)
!     print *, filestr4
     id_spec_deriv_vq_imag(p,1)		= register_diag_field(mod_name, filestr4,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
end do

output = 'I1_'

do p=1,12
   do q=1,12
     write(pstr ,"(I2.2)") p
     write(qstr ,"(I2.2)") q
!     print *, pstr, qstr
!     write(*,*) trim(output)//pstr//qstr
     write(filestr,*) adjustl(trim(output))//pstr//qstr
!     print *, filestr(2:8)
     filestr2 = filestr(2:8)
!     print *, filestr2
     id_I1(p,q)				= register_diag_field(mod_name, filestr2,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
end do
   end do

output = 'I2_'

do p=1,12
   do q=1,12
     write(pstr ,"(I2.2)") p
     write(qstr ,"(I2.2)") q
!     print *, pstr, qstr
!     write(*,*) trim(output)//pstr//qstr
     write(filestr,*) adjustl(trim(output))//pstr//qstr
!     print *, filestr(2:8)
     filestr2 = filestr(2:8)
!    print *, filestr2
     id_I2(p,q)				= register_diag_field(mod_name, filestr2,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
end do
   end do

output = 'I3_'

do p=1,12
   do q=1,12
     write(pstr ,"(I2.2)") p
     write(qstr ,"(I2.2)") q
!     print *, pstr, qstr
!     write(*,*) trim(output)//pstr//qstr
     write(filestr,*) adjustl(trim(output))//pstr//qstr
!     print *, filestr(2:8)
     filestr2 = filestr(2:8)
!     print *, filestr2
     id_I3(p,q)				= register_diag_field(mod_name, filestr2,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
  end do
end do

output = 'I4_'

do p=1,12
   do q=1,12
     write(pstr ,"(I2.2)") p
     write(qstr ,"(I2.2)") q
!     print *, pstr, qstr
!     write(*,*) trim(output)//pstr//qstr
     write(filestr,*) adjustl(trim(output))//pstr//qstr
!     print *, filestr(2:8)
     filestr2 = filestr(2:8)
!     print *, filestr2
     id_I4(p,q)				= register_diag_field(mod_name, filestr2,axis_2d, Time, 'vor with an fft taken along lon dim'	,'m/s'		)
   end do
end do

! JOE TERMS
!id_sink          		 	= register_diag_field(mod_name, 'sink' ,             	axis_2d, Time, 'sink'              			, 'm/s'      ) 
!id_sink_beta         			= register_diag_field(mod_name, 'sink_beta' ,        	axis_2d, Time, 'sink_beta'           			, 'm/s'      ) 
!id_robertssink       			= register_diag_field(mod_name, 'robertssink' ,      	axis_2d, Time, 'robertssink'         			, 'm/s'      ) 
!id_source            			= register_diag_field(mod_name, 'source' ,           	axis_2d, Time, 'source'              			, 'm/s'      ) 
!id_source_beta       			= register_diag_field(mod_name, 'source_beta' ,      	axis_2d, Time, 'source_beta'         			, 'm/s'      ) 
!id_ensflux           			= register_diag_field(mod_name, 'ensflux' ,          	axis_2d, Time, 'ensflux'             			, 'm/s'      ) 
!id_ensflux_div       			= register_diag_field(mod_name, 'ensflux_div' ,      	axis_2d, Time, 'ensflux_div'         			, 'm/s'      ) 
!id_ensflux_div_beta  			= register_diag_field(mod_name, 'ensflux_div_beta' , 	axis_2d, Time, 'ensflux_div_beta'    			, 'm/s'      ) 
!id_dens_dt        		   	= register_diag_field(mod_name, 'dens_dt' ,          	axis_2d, Time, 'dens_dt'             			, 'm/s'      ) 
!id_vvor_beta        		 	= register_diag_field(mod_name, 'vvor_beta' ,		axis_2d, Time, 'vvor_beta'           			, 'm/s'      ) 
id_vor    				= register_diag_field(mod_name, 'vor'   , 		axis_2d, Time, 'relative vorticity'  			, '1/s'      )
id_pv   				= register_diag_field(mod_name, 'pv'    , 		axis_2d, Time, 'absolute vorticity'  			, '1/s'      )
id_vq		   			= register_diag_field(mod_name, 'vq'    , 		axis_2d, Time, 'total meridional vorticity flux, vq'  	, 'm/s^2'    )
id_uq     				= register_diag_field(mod_name, 'uq'    , 		axis_2d, Time, 'total zonal vorticity flux, uq' 	, 'm/s^2'    )
id_stream		 		= register_diag_field(mod_name, 'stream', 		axis_2d, Time, 'streamfunction'     			, 'm^2/s'    )
id_trs    				= register_diag_field(mod_name, 'trs'   , 		axis_2d, Time, 'spectral tracer'     			, 'none'     )
id_tr 			    		= register_diag_field(mod_name, 'tr'    , 		axis_2d, Time, 'grid tracer'         			, 'none'     )

module_is_initialized = .true.

return
end subroutine barotropic_diagnostics_init

!--------------------------------------------------------------------------------------------

subroutine barotropic_diagnostics(Time, Grid, Phys, time_index)

type(time_type), 	intent(in) :: Time
type(phys_type), 	intent(in) :: Phys
type(grid_type), 	intent(in) :: Grid(1:12,1:12)
!type(testing_type), 	intent(in) :: Testing
integer,         	intent(in) :: time_index

logical :: used

if(.not.module_is_initialized) call error_mesg('barotropic_diagnostics', &
                              'barotropic_diagnostics_init has not been called', FATAL)

if(id_u      > 0) used = send_data(id_u      , Grid(1,1)%u     (:,:,time_index) , time)
if(id_v      > 0) used = send_data(id_v      , Grid(1,1)%v     (:,:,time_index) , time)
if(id_uprev  > 0) used = send_data(id_uprev      , Grid(1,1)%uprev     (:,:,time_index) , time)
if(id_ucurr  > 0) used = send_data(id_ucurr      , Grid(1,1)%ucurr     (:,:,time_index) , time)
if(id_vprev  > 0) used = send_data(id_vprev      , Grid(1,1)%vprev     (:,:,time_index) , time)
if(id_vcurr  > 0) used = send_data(id_vcurr      , Grid(1,1)%vcurr     (:,:,time_index) , time)
if(id_vorprev  > 0) used = send_data(id_vorprev      , Grid(1,1)%vorprev     (:,:,time_index) , time)
if(id_vorcurr  > 0) used = send_data(id_vorcurr      , Grid(1,1)%vorcurr     (:,:,time_index) , time)
if(id_utend  > 0) used = send_data(id_utend      , Grid(1,1)%utend     (:,:,time_index) , time)
if(id_vtend  > 0) used = send_data(id_vtend      , Grid(1,1)%vtend     (:,:,time_index) , time)
if(id_vortend      > 0) used = send_data(id_vortend      , Grid(1,1)%vortend     (:,:,time_index) , time)
if(id_vq     > 0) used = send_data(id_vq     , Grid(1,1)%vq    (:,:,time_index) , time)
if(id_uq     > 0) used = send_data(id_uq     , Grid(1,1)%uq    (:,:,time_index) , time)
if(id_dt_vor      > 0) used = send_data(id_dt_vor      , Grid(1,1)%dt_vor     (:,:,time_index) , time)
if(id_utrunc > 0) used = send_data(id_utrunc , Grid(1,1)%utrunc(:,:,time_index) , time)
if(id_vtrunc > 0) used = send_data(id_vtrunc , Grid(1,1)%vtrunc(:,:,time_index) , time)
if(id_vortrunc > 0) used = send_data(id_vortrunc , Grid(1,1)%vortrunc(:,:,time_index) , time)
if(id_ulindamp_m > 0) used = send_data(id_ulindamp_m , Grid(1,1)%ulindamp_m(:,:,time_index) , time)
if(id_ulindamp_e > 0) used = send_data(id_ulindamp_e , Grid(1,1)%ulindamp_e(:,:,time_index) , time)
if(id_vlindamp > 0) used = send_data(id_vlindamp , Grid(1,1)%vlindamp(:,:,time_index) , time)
if(id_vorlindamp_m > 0) used = send_data(id_vorlindamp_m , Grid(1,1)%vorlindamp_m(:,:,time_index) , time)
if(id_vorlindamp_e > 0) used = send_data(id_vorlindamp_e , Grid(1,1)%vorlindamp_e(:,:,time_index) , time)
if(id_dx_u2_v2 > 0) used = send_data(id_dx_u2_v2 , Grid(1,1)%dx_u2_v2(:,:,time_index) , time)
if(id_dy_u2_v2 > 0) used = send_data(id_dy_u2_v2 , Grid(1,1)%dy_u2_v2(:,:,time_index) , time)
if(id_ustir > 0) used = send_data(id_ustir , Grid(1,1)%ustir(:,:,time_index) , time)
if(id_vstir > 0) used = send_data(id_vstir , Grid(1,1)%vstir(:,:,time_index) , time)
if(id_vorstir > 0) used = send_data(id_vorstir , Grid(1,1)%vorstir(:,:,time_index) , time)
if(id_uspecdamp > 0) used = send_data(id_uspecdamp , Grid(1,1)%uspecdamp(:,:,time_index) , time)
if(id_vspecdamp > 0) used = send_data(id_vspecdamp , Grid(1,1)%vspecdamp(:,:,time_index) , time)
if(id_vorspecdamp > 0) used = send_data(id_vorspecdamp , Grid(1,1)%vorspecdamp(:,:,time_index) , time)
if(id_urobdamp > 0) used = send_data(id_urobdamp , Grid(1,1)%urobdamp(:,:,time_index) , time)
if(id_vrobdamp > 0) used = send_data(id_vrobdamp , Grid(1,1)%vrobdamp(:,:,time_index) , time)
if(id_vorrobdamp > 0) used = send_data(id_vorrobdamp , Grid(1,1)%vorrobdamp(:,:,time_index) , time)
if(id_pvor_advec > 0) used = send_data(id_pvor_advec , Grid(1,1)%pvor_advec(:,:,time_index) , time)
if(id_rvormean_advec > 0) used = send_data(id_rvormean_advec , Grid(1,1)%rvormean_advec(:,:,time_index) , time)
if(id_rvorprime_advec > 0) used = send_data(id_rvorprime_advec , Grid(1,1)%rvorprime_advec(:,:,time_index) , time)
if(id_rvor_advec > 0) used = send_data(id_rvor_advec , Grid(1,1)%rvor_advec(:,:,time_index) , time)
if(id_beta > 0) used = send_data(id_beta , Grid(1,1)%beta(:,:,time_index) , time)
if(id_beta_star > 0) used = send_data(id_beta_star , Grid(1,1)%beta_star(:,:,time_index) , time)
if(id_vormean_x > 0) used = send_data(id_vormean_x , Grid(1,1)%vormean_x(:,:,time_index) , time)
if(id_vormean_y > 0) used = send_data(id_vormean_y , Grid(1,1)%vormean_y(:,:,time_index) , time)
if(id_vorprime_x > 0) used = send_data(id_vorprime_x , Grid(1,1)%vorprime_x(:,:,time_index) , time)
if(id_vorprime_y > 0) used = send_data(id_vorprime_y , Grid(1,1)%vorprime_y(:,:,time_index) , time)
if(id_vormean_y_tend > 0) used = send_data(id_vormean_y_tend , Grid(1,1)%vormean_y_tend(:,:,time_index) , time)
if(id_vor_x > 0) used = send_data(id_vor_x , Grid(1,1)%vor_x(:,:,time_index) , time)
if(id_vor_y > 0) used = send_data(id_vor_y , Grid(1,1)%vor_y(:,:,time_index) , time)
if(id_vorprime > 0) used = send_data(id_vorprime , Grid(1,1)%vorprime(:,:,time_index) , time)
if(id_vorprime_prev > 0) used = send_data(id_vorprime_prev , Grid(1,1)%vorprime_prev(:,:,time_index) , time)
if(id_vorprime_curr > 0) used = send_data(id_vorprime_curr , Grid(1,1)%vorprime_curr(:,:,time_index) , time)
if(id_vork > 0) used = send_data(id_vork , Grid(1,1)%vork(:,:,time_index) , time)
if(id_vork_prev > 0) used = send_data(id_vork_prev , Grid(1,1)%vork_prev(:,:,time_index) , time)
if(id_vork_curr > 0) used = send_data(id_vork_curr , Grid(1,1)%vork_curr(:,:,time_index) , time)
if(id_vortendk > 0) used = send_data(id_vortendk , Grid(1,1)%vortendk(:,:,time_index) , time)
if(id_vormeantend > 0) used = send_data(id_vormeantend , Grid(1,1)%vormeantend(:,:,time_index) , time)
if(id_vorprimetend > 0) used = send_data(id_vorprimetend , Grid(1,1)%vorprimetend(:,:,time_index) , time)
if(id_enstrophy_rvormean_adevc_Matlab > 0) used = send_data(id_enstrophy_rvormean_adevc_Matlab , Grid(1,1)%enstrophy_rvormean_adevc_Matlab(:,:,time_index) , time)
if(id_enstrophy_rvorprime_adevc_Matlab > 0) used = send_data(id_enstrophy_rvorprime_adevc_Matlab , Grid(1,1)%enstrophy_rvorprime_adevc_Matlab(:,:,time_index) , time)
if(id_tester > 0) used = send_data(id_tester , Grid(1,1)%tester(:,:,time_index) , time)
if(id_tester2 > 0) used = send_data(id_tester2 , Grid(1,1)%tester2(:,:,time_index) , time)
if(id_tester3 > 0) used = send_data(id_tester3 , Grid(1,1)%tester3(:,:,time_index) , time)
if(id_tester4 > 0) used = send_data(id_tester4 , Grid(1,1)%tester4(:,:,time_index) , time)
if(id_tester5 > 0) used = send_data(id_tester5 , Grid(1,1)%tester5(:,:,time_index) , time)
if(id_tester6 > 0) used = send_data(id_tester6 , Grid(1,1)%tester6(:,:,time_index) , time)
if(id_v_y > 0) used = send_data(id_v_y , Grid(1,1)%v_y(:,:,time_index) , time)
if(id_rvor_adveck > 0) used = send_data(id_rvor_adveck , Grid(1,1)%rvor_adveck(:,:,time_index) , time)
if(id_u_dvordx > 0) used = send_data(id_u_dvordx , Grid(1,1)%u_dvordx(:,:,time_index) , time)
if(id_v_dvordy > 0) used = send_data(id_v_dvordy , Grid(1,1)%v_dvordy(:,:,time_index) , time)
if(id_vor_dudx > 0) used = send_data(id_u_dvordx , Grid(1,1)%vor_dudx(:,:,time_index) , time)
if(id_vor_dvdy > 0) used = send_data(id_u_dvordx , Grid(1,1)%vor_dvdy(:,:,time_index) , time)
! ENERGY TERMS
if(id_energy > 0) used = send_data(id_energy , Grid(1,1)%energy(:,:,time_index) , time)
if(id_energy_tend > 0) used = send_data(id_energy_tend , Grid(1,1)%energy_tend(:,:,time_index) , time)
if(id_energy_voradvec > 0) used = send_data(id_energy_voradvec , Grid(1,1)%energy_voradvec(:,:,time_index) , time)
if(id_energy_gradterm > 0) used = send_data(id_energy_gradterm , Grid(1,1)%energy_gradterm(:,:,time_index) , time)
if(id_energy_trunc > 0) used = send_data(id_energy_trunc , Grid(1,1)%energy_trunc(:,:,time_index) , time)
if(id_energy_lindamp_m > 0) used = send_data(id_energy_lindamp_m , Grid(1,1)%energy_lindamp_m(:,:,time_index) , time)
if(id_energy_lindamp_e > 0) used = send_data(id_energy_lindamp_e , Grid(1,1)%energy_lindamp_e(:,:,time_index) , time)
if(id_energy_specdamp > 0) used = send_data(id_energy_specdamp , Grid(1,1)%energy_specdamp(:,:,time_index) , time)
if(id_energy_stir > 0) used = send_data(id_energy_stir , Grid(1,1)%energy_stir(:,:,time_index) , time)
if(id_energy_robdamp > 0) used = send_data(id_energy_robdamp , Grid(1,1)%energy_robdamp(:,:,time_index) , time)
! MEAN ENERGY TERMS
if(id_energy_mean > 0) used = send_data(id_energy_mean , Grid(1,1)%energy_mean(:,:,time_index) , time)
if(id_energy_tend_mean > 0) used = send_data(id_energy_tend_mean , Grid(1,1)%energy_tend_mean(:,:,time_index) , time)
if(id_energy_voradvec_mean > 0) used = send_data(id_energy_voradvec_mean , Grid(1,1)%energy_voradvec_mean(:,:,time_index) , time)
if(id_energy_gradterm_mean > 0) used = send_data(id_energy_gradterm_mean , Grid(1,1)%energy_gradterm_mean(:,:,time_index) , time)
if(id_energy_trunc_mean > 0) used = send_data(id_energy_trunc_mean , Grid(1,1)%energy_trunc_mean(:,:,time_index) , time)
if(id_energy_lindamp_m_mean > 0) used = send_data(id_energy_lindamp_m_mean , Grid(1,1)%energy_lindamp_m_mean(:,:,time_index) , time)
if(id_energy_lindamp_e_mean > 0) used = send_data(id_energy_lindamp_e_mean , Grid(1,1)%energy_lindamp_e_mean(:,:,time_index) , time)
if(id_energy_specdamp_mean > 0) used = send_data(id_energy_specdamp_mean , Grid(1,1)%energy_specdamp_mean(:,:,time_index) , time)
if(id_energy_stir_mean > 0) used = send_data(id_energy_stir_mean , Grid(1,1)%energy_stir_mean(:,:,time_index) , time)
if(id_energy_robdamp_mean > 0) used = send_data(id_energy_robdamp_mean , Grid(1,1)%energy_robdamp_mean(:,:,time_index) , time)
! ENSTROPHY TERMS
if(id_enstrophy > 0) used = send_data(id_enstrophy , Grid(1,1)%enstrophy(:,:,time_index) , time)
if(id_enstrophy_tend > 0) used = send_data(id_enstrophy_tend , Grid(1,1)%enstrophy_tend(:,:,time_index) , time)
if(id_enstrophy_dt_vor > 0) used = send_data(id_enstrophy_dt_vor , Grid(1,1)%enstrophy_dt_vor(:,:,time_index) , time)
if(id_enstrophy_pvor_advec > 0) used = send_data(id_enstrophy_pvor_advec , Grid(1,1)%enstrophy_pvor_advec(:,:,time_index) , time)
if(id_enstrophy_rvormean_advec > 0) used = send_data(id_enstrophy_rvormean_advec , Grid(1,1)%enstrophy_rvormean_advec(:,:,time_index) , time)
if(id_enstrophy_rvorprime_advec > 0) used = send_data(id_enstrophy_rvorprime_advec , Grid(1,1)%enstrophy_rvorprime_advec(:,:,time_index) , time)
if(id_enstrophy_trunc > 0) used = send_data(id_enstrophy_trunc , Grid(1,1)%enstrophy_trunc(:,:,time_index) , time)
if(id_enstrophy_lindamp_m > 0) used = send_data(id_enstrophy_lindamp_m , Grid(1,1)%enstrophy_lindamp_m(:,:,time_index) , time)
if(id_enstrophy_lindamp_e > 0) used = send_data(id_enstrophy_lindamp_e , Grid(1,1)%enstrophy_lindamp_e(:,:,time_index) , time)
if(id_enstrophy_specdamp > 0) used = send_data(id_enstrophy_specdamp , Grid(1,1)%enstrophy_specdamp(:,:,time_index) , time)
if(id_enstrophy_stir > 0) used = send_data(id_enstrophy_stir , Grid(1,1)%enstrophy_stir(:,:,time_index) , time)
if(id_enstrophy_robdamp > 0) used = send_data(id_enstrophy_robdamp , Grid(1,1)%enstrophy_robdamp(:,:,time_index) , time)
if(id_enstrophy_tendk > 0) used = send_data(id_enstrophy_tendk , Grid(1,1)%enstrophy_tendk(:,:,time_index) , time)
if(id_enstrophy_pvor_adveck > 0) used = send_data(id_enstrophy_pvor_adveck , Grid(1,1)%enstrophy_pvor_adveck(:,:,time_index) , time)
if(id_enstrophy_rvormean_adveck > 0) used = send_data(id_enstrophy_rvormean_adveck , Grid(1,1)%enstrophy_rvormean_adveck(:,:,time_index) , time)
if(id_enstrophy_trunck > 0) used = send_data(id_enstrophy_trunck , Grid(1,1)%enstrophy_trunck(:,:,time_index) , time)
if(id_enstrophy_lindamp_ek > 0) used = send_data(id_enstrophy_lindamp_ek , Grid(1,1)%enstrophy_lindamp_ek(:,:,time_index) , time)
if(id_enstrophy_lindamp_mk > 0) used = send_data(id_enstrophy_lindamp_mk , Grid(1,1)%enstrophy_lindamp_mk(:,:,time_index) , time)
if(id_enstrophy_specdampk > 0) used = send_data(id_enstrophy_specdampk , Grid(1,1)%enstrophy_specdampk(:,:,time_index) , time)
if(id_enstrophy_robdampk > 0) used = send_data(id_enstrophy_robdampk , Grid(1,1)%enstrophy_robdampk(:,:,time_index) , time)
if(id_enstrophy_stirk > 0) used = send_data(id_enstrophy_stirk , Grid(1,1)%enstrophy_stirk(:,:,time_index) , time)
if(id_enstrophy_rvorprime_adveck_SUM > 0) used = send_data(id_enstrophy_rvorprime_adveck_SUM , Grid(1,1)%enstrophy_rvorprime_adveck_SUM(:,:,time_index) , time)
if(id_I_higher_freq > 0) used = send_data(id_I_higher_freq , Grid(1,1)%I_higher_freq(:,:,time_index) , time)
if(id_vorprime_yk > 0) used = send_data(id_vorprime_yk , Grid(1,1)%vorprime_yk(:,:,time_index) , time)
if(id_v_yk > 0) used = send_data(id_v_yk , Grid(1,1)%v_yk(:,:,time_index) , time)
if(id_vk > 0) used = send_data(id_vk , Grid(1,1)%vk(:,:,time_index) , time)
if(id_vorprimek_curr > 0) used = send_data(id_vorprimek_curr , Grid(1,1)%vorprimek_curr(:,:,time_index) , time)

! FLUX TERMS
if(id_vorprime_fluxk > 0) used = send_data(id_vorprime_fluxk , Grid(1,1)%vorprime_fluxk(:,:,time_index) , time)
if(id_I_pp > 0) used = send_data(id_I_pp , Grid(1,1)%I_pp(:,:,time_index) , time)
if(id_I_sp > 0) used = send_data(id_I_sp , Grid(1,1)%I_sp(:,:,time_index) , time)
if(id_I_ps > 0) used = send_data(id_I_ps , Grid(1,1)%I_ps(:,:,time_index) , time)
if(id_I_hp > 0) used = send_data(id_I_hp , Grid(1,1)%I_hp(:,:,time_index) , time)
if(id_I_ph > 0) used = send_data(id_I_ph , Grid(1,1)%I_ph(:,:,time_index) , time)
if(id_I_ss > 0) used = send_data(id_I_ss , Grid(1,1)%I_ss(:,:,time_index) , time)
if(id_I_hh > 0) used = send_data(id_I_hh , Grid(1,1)%I_hh(:,:,time_index) , time)
if(id_I_sh > 0) used = send_data(id_I_sh , Grid(1,1)%I_sh(:,:,time_index) , time)
if(id_I_hs > 0) used = send_data(id_I_hs , Grid(1,1)%I_hs(:,:,time_index) , time)
!if(id_I2_pp > 0) used = send_data(id_I_pp , Grid(1,1)%I_pp(:,:,time_index) , time)
!if(id_I2_sp > 0) used = send_data(id_I_sp , Grid(1,1)%I_sp(:,:,time_index) , time)
!if(id_I2_ps > 0) used = send_data(id_I_ps , Grid(1,1)%I_ps(:,:,time_index) , time)
!if(id_I2_hp > 0) used = send_data(id_I_hp , Grid(1,1)%I_hp(:,:,time_index) , time)
!if(id_I2_ph > 0) used = send_data(id_I_ph , Grid(1,1)%I_ph(:,:,time_index) , time)
!if(id_I2_ss > 0) used = send_data(id_I_ss , Grid(1,1)%I_ss(:,:,time_index) , time)
!if(id_I2_hh > 0) used = send_data(id_I_hh , Grid(1,1)%I_hh(:,:,time_index) , time)
!if(id_I2_sh > 0) used = send_data(id_I_sh , Grid(1,1)%I_sh(:,:,time_index) , time)
!if(id_I2_hs > 0) used = send_data(id_I_hs , Grid(1,1)%I_hs(:,:,time_index) , time)
if(id_I_0p > 0) used = send_data(id_I_0p , Grid(1,1)%I_0p(:,:,time_index) , time)
if(id_I_p0 > 0) used = send_data(id_I_p0 , Grid(1,1)%I_p0(:,:,time_index) , time)
if(id_I_0s > 0) used = send_data(id_I_0s , Grid(1,1)%I_0s(:,:,time_index) , time)
if(id_I_s0 > 0) used = send_data(id_I_s0 , Grid(1,1)%I_s0(:,:,time_index) , time)
if(id_I_0h > 0) used = send_data(id_I_0h , Grid(1,1)%I_0h(:,:,time_index) , time)
if(id_I_h0 > 0) used = send_data(id_I_h0 , Grid(1,1)%I_h0(:,:,time_index) , time)
if(id_I_00 > 0) used = send_data(id_I_00 , Grid(1,1)%I_00(:,:,time_index) , time)
if(id_I_XX > 0) used = send_data(id_I_XX , Grid(1,1)%I_XX(:,:,time_index) , time)

!do p=1,12
!   do q = 1,12
!       if(id_I(p,q) > 0) used = send_data(id_I(p,q) , Grid(p,q)%I(:,:,time_index) , time)
!   end do
!end do

!if(id_I > 0) used = send_data(id_I , Grid(1,1)%I(:,:,time_index) , time)
!if(id_I2 > 0) used = send_data(id_I2 , Grid(1,2)%I(:,:,time_index) , time)

!if(id_I(1,1) > 0) used = send_data(id_I(1,1) , Grid(1,1)%I(:,:,time_index) , time)
!if(id_I(1,2) > 0) used = send_data(id_I(1,2) , Grid(1,2)%I(:,:,time_index) , time)

!if(id_spec_deriv_vq_real(p,1) > 0) used = send_data(id_spec_deriv_vq_real(p,1) , Grid(p,1)%spec_deriv_vq_real(:,:,time_index) , time)
!if(id_spec_deriv_vq_imag(p,1) > 0) used = send_data(id_spec_deriv_vq_imag(p,1) , Grid(p,1)%spec_deriv_vq_imag(:,:,time_index) , time)

do p=1,128
       if(id_spec_deriv_vq_real(p,1) > 0) used = send_data(id_spec_deriv_vq_real(p,1) , Grid(p,1)%spec_deriv_vq_real(:,:,time_index) , time)
       if(id_spec_deriv_vq_imag(p,1) > 0) used = send_data(id_spec_deriv_vq_imag(p,1) , Grid(p,1)%spec_deriv_vq_imag(:,:,time_index) , time)
end do

do p=1,12
   do q = 1,12
       if(id_I1(p,q) > 0) used = send_data(id_I1(p,q) , Grid(p,q)%I1(:,:,time_index) , time)
   end do
end do

do p=1,12
   do q = 1,12
       if(id_I2(p,q) > 0) used = send_data(id_I2(p,q) , Grid(p,q)%I2(:,:,time_index) , time)
   end do
end do

do p=1,12
   do q = 1,12
       if(id_I3(p,q) > 0) used = send_data(id_I3(p,q) , Grid(p,q)%I3(:,:,time_index) , time)
   end do
end do

do p=1,12
   do q = 1,12
       if(id_I4(p,q) > 0) used = send_data(id_I4(p,q) , Grid(p,q)%I4(:,:,time_index) , time)
   end do
end do

!if(id_f_test > 0) used = send_data(id_f_test , Testing%f_test(:,:,:,time_index) , time)
!if(id_f_test  > 0) used = send_data(id_f_test,  f_test, Time)
!if(id_enstrophy_rvorprime_adveck > 0) used = send_data(id_enstrophy_rvorprime_adveck , Testing%enstrophy_rvorprime_adveck(:,:,:,time_index) , time)

!print *, "ENSTIRK", Grid(1,1)%enstrophy_stirk(22,22,time_index)
!print *, "F_TEST", Testing%f_test(22,22,22,time_index)

!do i = is, ie 
!   if(id_enstrophy_rvorprime_adveck > 0) used = send_data(id_enstrophy_rvorprime_adveck, Testing%enstrophy_rvorprime_adveck(:,:,time_index) , time)
!end do

! JOE TERMS
!if(id_sink                  > 0) used = send_data(id_sink      , Grid(1,1)%sink     (:,:,time_index) , time)
!if(id_sink_beta             > 0) used = send_data(id_sink_beta      , Grid(1,1)%sink_beta     (:,:,time_index) , time)
!if(id_robertssink           > 0) used = send_data(id_robertssink      , Grid(1,1)%robertssink     (:,:,time_index) , time)
!if(id_source                > 0) used = send_data(id_source      , Grid(1,1)%source     (:,:,time_index) , time)
!if(id_source_beta           > 0) used = send_data(id_source_beta      , Grid(1,1)%source_beta     (:,:,time_index) , time)
!if(id_ensflux               > 0) used = send_data(id_ensflux      , Grid(1,1)%ensflux     (:,:,time_index) , time)
!if(id_ensflux_div           > 0) used = send_data(id_ensflux_div      , Grid(1,1)%ensflux_div     (:,:,time_index) , time)
!if(id_ensflux_div_beta      > 0) used = send_data(id_ensflux_div_beta      , Grid(1,1)%ensflux_div_beta     (:,:,time_index) , time)
!if(id_dens_dt               > 0) used = send_data(id_dens_dt      , Grid(1,1)%dens_dt     (:,:,time_index) , time)
!if(id_vvor_beta      > 0) used = send_data(id_vvor_beta      , Grid(1,1)%vvor_beta     (:,:,time_index) , time)
if(id_vor    > 0) used = send_data(id_vor    , Grid(1,1)%vor   (:,:,time_index) , time)
if(id_pv     > 0) used = send_data(id_pv     , Grid(1,1)%pv    (:,:)            , time)
if(id_stream > 0) used = send_data(id_stream , Grid(1,1)%stream(:,:)            , time)
if(id_tr     > 0) used = send_data(id_tr     , Grid(1,1)%tr    (:,:,time_index) , time)
if(id_trs    > 0) used = send_data(id_trs    , Grid(1,1)%trs   (:,:,time_index) , time)

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
