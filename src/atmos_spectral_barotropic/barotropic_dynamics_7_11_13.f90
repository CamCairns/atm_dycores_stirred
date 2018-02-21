! ==============  BIG COMMENT ====================
!
! In this m,odified code I tried to apply the vor_div_from_uv_grid subroutine to the constituent fourier coefficent products in an effort to calculate (vq)_y accurately in spectral space. This did not work
! There is a check in the code (comment "COMMENT ON vor_div_from_uv_grid SUBROUTINE") where I confirm that there is a numerical error in the product rule, ie.
!				(AB)_y =/= A(B)_y + (A)_yB
! The reason for this is that if we assume some error "a" associated with any derivative (that is not the original derivative vor_div_from_uv_grid(-vq,uq) which we hold as "truth") then the product rule will incur an error A*a + a*B. In the example of calculating the vorticity flux using the product rule this error amounts to around 10% (if we note that we then multiply this value further by zeta the this is consistent with the residuals found in the Interaction coefficent calculated).
!
! Calculating the derivative by first taking an fft will not assuage this. This is because we are now taking (for a product of two fields) N^2 derivatives and incurring N^2 errors. Of course these errors will average out but there is till a significant problem (in excess of the error in the simple product rule). This is clear by comparing the calculations
!
!				FFT(vor_div_from_uv_grid(A))	&	vor_div_from_uv_grid(FFT(A))
!
! Where A is some field, such as u or v.
!
! BOTTOM LINE: There is no tricky wave of calculating the derivative accurately and being able to calculate in spectral space. Solution: Redefine what the derivative is so that it is the product rule version (RHS of 1st eq) that calculates how q develops in the model. I believe this is defensible and might work.

module barotropic_dynamics_mod

use               fms_mod, only: open_namelist_file,   &
                                 open_restart_file,    &
                                 file_exist,           &
                                 check_nml_error,      &
                                 error_mesg,           &
                                 FATAL, WARNING,       &
                                 write_version_number, &
                                 mpp_pe,               &
                                 mpp_root_pe,          &
                                 read_data,            &
                                 write_data,           &
                                 set_domain,           &
                                 close_file,           &
                                 stdlog

use    time_manager_mod,  only : time_type,      &
                                 get_time,       &
                                 operator(==),   &
                                 operator(-)

use         constants_mod,  only: radius, omega, pi

use         transforms_mod, only: transforms_init,         transforms_end,          &
                                  get_grid_boundaries,                              &
                                  trans_spherical_to_grid, trans_grid_to_spherical, &  
                                  compute_laplacian,                                &
                                  get_sin_lat,             get_cos_lat,             &
                                  get_deg_lon,             get_deg_lat,             &
                                  get_grid_domain,         get_spec_domain,         &
                                  spectral_domain,         grid_domain,             &
                                  vor_div_from_uv_grid,    uv_grid_from_vor_div,    &
                                  horizontal_advection

use   spectral_damping_mod, only: spectral_damping_init,  &
                                  compute_spectral_damping

use           leapfrog_mod, only: leapfrog

use       fv_advection_mod, only: fv_advection_init, &
                                  a_grid_horiz_advection

use           stirring_mod, only: stirring_init, stirring, stirring_end

use fft_mod,       only: fft_init, fft_end, fft_grid_to_fourier, fft_fourier_to_grid

!===============================================================================================
implicit none
private
!===============================================================================================
logical :: used

public :: barotropic_dynamics_init, &
          barotropic_dynamics,      &
          barotropic_dynamics_end,  &
          dynamics_type,            &
          grid_type,                &
          spectral_type,            &
          tendency_type


! version information 
!===================================================================
character(len=128) :: version = '$Id: barotropic_dynamics.f90,v 10.0 2003/10/24 22:00:57 fms Exp $'
character(len=128) :: tagname = '$Name: latest $'
!===================================================================

type grid_type
! Original model variables
   real, pointer, dimension(:,:,:) :: u=>NULL(), v=>NULL(), vor=>NULL(), trs=>NULL(), tr=>NULL(), &
! Momentum, Vorticity budget variables
utend=>NULL(), vtend=>NULL(), vortend=>NULL(), vq=>NULL(), uq=>NULL(), utrunc=>NULL(), vtrunc=>NULL(), & vortrunc=>NULL(), dt_vor=>NULL(), vlindamp=>NULL(), &
ulindamp_e=>NULL(), ulindamp_m=>NULL(), vorlindamp_e=>NULL(), vorlindamp_m=>NULL(), uspecdamp=>NULL(), vspecdamp=>NULL(), vorspecdamp=>NULL(),ustir=>NULL(), vstir=>NULL(), vorstir=>NULL(), &
urobdamp=>NULL(), vrobdamp=>NULL(), vorrobdamp=>NULL(), dx_u2_v2=>NULL(), dy_u2_v2=>NULL(), rvor_advec=>NULL(), rvormean_advec=>NULL(), rvorprime_advec=>NULL(), pvor_advec=>NULL(), beta=>NULL(), &
beta_star=>NULL(), tester=>NULL(), tester2=>NULL(), tester3=>NULL(), vormean_x=>NULL(), vormean_y=>NULL(), vorprime_x=>NULL(), vorprime_y=>NULL(), vormean_y_tend=>NULL(), vor_y=>NULL(), vor_x=>NULL(), & 
enstrophy_rvorprime_adevc_Matlab=>NULL(), enstrophy_rvormean_adevc_Matlab=>NULL(), vorprime_prev=>NULL(), vorprime_curr=>NULL(), vorprime=>NULL(), vork=>NULL(), vork_prev=>NULL(), vork_curr=>NULL(), &
vortendk=>NULL(), vormeantend=>NULL(), vorprimetend=>NULL(), I_XX=>NULL(), v_y=>NULL(), &
! Some more variables
uprev=>NULL(), ucurr=>NULL(), vprev=>NULL(), vcurr=>NULL(), vorprev=>NULL(), vorcurr=>NULL(), rvor_adveck=>NULL(), &
! Energy
energy=>NULL(), energy_tend=>NULL(), energy_voradvec=>NULL(), energy_gradterm=>NULL(), energy_trunc=>NULL(), energy_lindamp_m=>NULL(), energy_lindamp_e=>NULL(), energy_specdamp=>NULL(), & 
energy_stir=>NULL(), energy_robdamp=>NULL(), energy_tend_mean=>NULL(), energy_voradvec_mean=>NULL(), energy_gradterm_mean=>NULL(), energy_trunc_mean=>NULL(), energy_lindamp_m_mean=>NULL(), &
energy_lindamp_e_mean=>NULL(), energy_specdamp_mean=>NULL(), energy_stir_mean=>NULL(), energy_robdamp_mean=>NULL(), energy_mean=>NULL(), &
! Enstrophy
enstrophy=>NULL(), enstrophy_tend=>NULL(), enstrophy_dt_vor=>NULL(), enstrophy_pvor_advec=>NULL(), enstrophy_rvormean_advec=>NULL(), enstrophy_rvorprime_advec=>NULL(), enstrophy_trunc=>NULL(), enstrophy_lindamp_m=>NULL(), enstrophy_lindamp_e=>NULL(), enstrophy_specdamp=>NULL(), enstrophy_stir=>NULL(), enstrophy_robdamp=>NULL(), vorprime_yk=>NULL(), v_yk=>NULL(), vorprimek_curr=>NULL(), &
vk=>NULL(), &
enstrophy_tendk=>NULL(), enstrophy_pvor_adveck=>NULL(), enstrophy_rvormean_adveck=>NULL(), enstrophy_trunck=>NULL(), enstrophy_lindamp_ek=>NULL(), enstrophy_lindamp_mk=>NULL(), enstrophy_specdampk=>NULL(),&
enstrophy_robdampk=>NULL(), enstrophy_stirk=>NULL(), enstrophy_rvorprime_adveck_SUM=>NULL(), I_higher_freq=>NULL(), &
! Interaction
I_pp=>NULL(), I_sp=>NULL(), I_ps=>NULL(), I_hp=>NULL(), I_ph=>NULL(), I_ss=>NULL(), I_hh=>NULL(), I_sh=>NULL(), I_hs=>NULL(), I_0p=>NULL(), I_p0=>NULL(), I_0s=>NULL(), I_s0=>NULL(), I_0h=>NULL(), &
I_h0=>NULL(), I_00=>NULL(), &
I2_pp=>NULL(), I2_sp=>NULL(), I2_ps=>NULL(), I2_hp=>NULL(), I2_ph=>NULL(), I2_ss=>NULL(), I2_hh=>NULL(), I2_sh=>NULL(), I2_hs=>NULL(), &
! FLUX TERM
vorprime_fluxk=>NULL(), &
! Joe variables
vvor_beta=>NULL(), source=>NULL(), source_beta=>NULL(), sink=>NULL(), sink_beta=>NULL(), dens_dt=>NULL(), robertssink=>NULL(), ensflux=>NULL(), ensflux_div=>NULL(),ensflux_div_beta=>NULL()
!real, pointer :: I(3,3) 
real, pointer, dimension(:,:)   :: pv=>NULL(), stream=>NULL()
real, POINTER, DIMENSION(:,:,:) :: I1, I2, I3, I4, spec_deriv_vq_real=>NULL(), spec_deriv_vq_imag=>NULL()
!real, POINTER, DIMENSION(:,:,:) :: negI
end type


TYPE(grid_type), ALLOCATABLE, DIMENSION(:,:) :: Grid

type spectral_type
   complex, pointer, dimension(:,:,:) :: vor=>NULL(), trs=>NULL(), roberts=>NULL()
end type
type tendency_type
   real, pointer, dimension(:,:) :: u=>NULL(), v=>NULL(), trs=>NULL(), tr=>NULL()
end type

!type testing_type
!   real, dimension(:,:,:), allocatable :: I(3)!=>NULL()
!   allocatable, dimension(:) testing
!end type

!type testing_type
! dimension(3) :: testing
!Dyn%testing(1)%I => NULL()
!Dyn%testing(2)%I => NULL()
!Dyn%testing(3)%I => NULL()
!end type 

integer :: is, ie, js, je, ms, me, ns, ne, icnt, p, q

type dynamics_type
   type(grid_type)     :: grid(1:12,1:12)
   type(spectral_type) :: spec
   type(tendency_type) :: tend
!   type(testing_type)  :: testing(1:12,1:12)
   integer             :: num_lon, num_lat
   logical             :: grid_tracer, spec_tracer
end type

integer, parameter :: num_time_levels = 2

logical :: module_is_initialized = .false.

real,  allocatable, dimension(:)   :: sin_lat, cos_lat, rad_lat, rad_lon, &
                                      deg_lat, deg_lon, u_init, v_init, &
                                      coriolis, glon_bnd, glat_bnd

integer :: pe, npes

! namelist parameters with default values

logical  :: check_fourier_imag = .false.
logical  :: south_to_north     = .true.
logical  :: triang_trunc       = .true.

real     :: robert_coeff       = 0.04
real     :: longitude_origin   = 0.0

character(len=64) :: damping_option = 'resolution_dependent'
integer  :: damping_order      = 4
real     :: damping_coeff      = 1.e-04

real     :: zeta_0     = 8.e-05
real     :: tau        = 691200.
real     :: tau1       = 691200.
integer  :: m_0        = 4
integer  :: read_file  = 1
integer  :: damp       = 1.0
integer  :: spec_calc = 0
real  	 :: mult       = 0.0
real     :: linmult    = 0.0
real     :: eddy_width = 15.0
real     :: eddy_lat   = 45.0

logical  :: spec_tracer      = .true.
logical  :: grid_tracer      = .true.

integer  :: num_lat            = 128
integer  :: num_lon            = 256
integer  :: num_fourier        = 85
integer  :: num_spherical      = 86
integer  :: fourier_inc        = 1

real, dimension(2) :: valid_range_v = (/-1.e3,1.e3/)

namelist /barotropic_dynamics_nml/ check_fourier_imag, south_to_north, &
                                 triang_trunc,                         &
                                 num_lon, num_lat, num_fourier,        &
                                 num_spherical, fourier_inc,           &
                                 longitude_origin, damping_option,     &
                                 damping_order,     damping_coeff,     &
                                 robert_coeff,                         &
                                 spec_tracer, grid_tracer,             &
                                 eddy_lat, eddy_width, zeta_0, m_0, damp, mult, linmult, tau, tau1, spec_calc, read_file, &
                                 valid_range_v

contains

!===============================================================================================

subroutine barotropic_dynamics_init (Dyn,  Time, Time_init, dt_real, id_lon, id_lat, id_lonb, id_latb)

type(dynamics_type), intent(inout)  :: Dyn
type(time_type)    , intent(in)     :: Time, Time_init
real, intent(in) :: dt_real
integer, intent(out) :: id_lon, id_lat, id_lonb, id_latb

integer :: i, j, p, q, k

real,    allocatable, dimension(:)   :: glon_bnd, glat_bnd
complex, allocatable, dimension(:,:) :: div
real :: xx, yy, dd

integer  :: ierr, io, unit, pe
logical  :: root

! < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < >

call write_version_number (version, tagname)

pe   =  mpp_pe()
root = (pe == mpp_root_pe())

if (file_exist('input.nml')) then
  unit = open_namelist_file ()
  ierr=1
  do while (ierr /= 0)
    read  (unit, nml=barotropic_dynamics_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'barotropic_dynamics_nml')
  enddo
  10 call close_file (unit)
endif


if (root) write (stdlog(), nml=barotropic_dynamics_nml)

 call transforms_init(radius, num_lat, num_lon, num_fourier, fourier_inc, num_spherical, &
         south_to_north=south_to_north,   &
         triang_trunc=triang_trunc,       &
         longitude_origin=longitude_origin       )

 call get_grid_domain(is,ie,js,je)
 call get_spec_domain(ms,me,ns,ne)

Dyn%num_lon      = num_lon
Dyn%num_lat      = num_lat
Dyn%spec_tracer  = spec_tracer
Dyn%grid_tracer  = grid_tracer

allocate (sin_lat   (js:je))
allocate (cos_lat   (js:je))
allocate (deg_lat   (js:je))
allocate (deg_lon   (is:ie))
allocate (rad_lat   (js:je))
allocate (rad_lon   (is:ie))
allocate (coriolis  (js:je))

allocate (glon_bnd  (num_lon + 1))
allocate (glat_bnd  (num_lat + 1))

 call get_deg_lon (deg_lon)
 call get_deg_lat (deg_lat)
 call get_sin_lat (sin_lat)
 call get_cos_lat (cos_lat)
 call get_grid_boundaries (glon_bnd, glat_bnd, global=.true.)

coriolis = 2*omega*sin_lat

rad_lat = deg_lat*atan(1.0)/45.0
rad_lon = deg_lon*atan(1.0)/45.0

 call spectral_damping_init(damping_coeff, damping_order, damping_option, num_fourier, num_spherical, 1, 0., 0., 0.)

 call stirring_init(dt_real, Time, id_lon, id_lat, id_lonb, id_latb)

allocate (Dyn%spec%vor (ms:me, ns:ne, num_time_levels))
allocate (Dyn%grid(1,1)%u   (is:ie, js:je, num_time_levels))
allocate (Dyn%grid(1,1)%v   (is:ie, js:je, num_time_levels))
allocate (Dyn%grid(1,1)%vor (is:ie, js:je, num_time_levels))

allocate (Dyn%tend%u        (is:ie, js:je))
allocate (Dyn%tend%v        (is:ie, js:je))

allocate (Dyn%Grid(1,1)%stream   (is:ie, js:je))
allocate (Dyn%Grid(1,1)%pv       (is:ie, js:je))

!---------------------- Extra terms -------------------------------------------
! Momentum and Vorticity Budget Terms
allocate (Dyn%Grid(1,1)%ucurr		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%uprev		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vcurr		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vprev		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorcurr		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprev		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%utend		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vtend		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vortend		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vq			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%uq			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%dt_vor		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vlindamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%ulindamp_e		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%ulindamp_m		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorlindamp_e		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorlindamp_m		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%utrunc		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vtrunc		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vortrunc		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%dx_u2_v2		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%dy_u2_v2		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%uspecdamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vspecdamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorspecdamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%ustir		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vstir		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorstir		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%urobdamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vrobdamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorrobdamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%rvor_advec		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%rvormean_advec	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%rvorprime_advec	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%pvor_advec		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%beta			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%beta_star		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vormean_x		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vormean_y		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprime_x		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprime_y		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vormean_y_tend	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vor_y		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vor_x		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%v_y		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprime_prev	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprime_curr	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprime		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vork			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vork_prev		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vork_curr		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vortendk		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vormeantend		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprimetend		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_rvormean_adevc_Matlab(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_rvorprime_adevc_Matlab(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%tester			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%tester2			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%tester3			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_XX			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%rvor_adveck		(is:ie, js:je, num_time_levels))

allocate (Dyn%spec%roberts		(ms:me, ns:ne, num_time_levels))


! Energy Equation Terms
allocate (Dyn%Grid(1,1)%energy		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_tend		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_voradvec	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_gradterm	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_trunc		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_lindamp_m	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_lindamp_e	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_specdamp	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_stir		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_robdamp	(is:ie, js:je, num_time_levels))

allocate (Dyn%Grid(1,1)%energy_mean		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_tend_mean	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_voradvec_mean	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_gradterm_mean	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_trunc_mean	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_lindamp_m_mean(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_lindamp_e_mean(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_specdamp_mean	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_stir_mean	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%energy_robdamp_mean	(is:ie, js:je, num_time_levels))

! Enstrophy Equation Terms
allocate (Dyn%Grid(1,1)%enstrophy			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_tend			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_dt_vor		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_pvor_advec		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_rvormean_advec	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_rvorprime_advec	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_trunc			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_lindamp_e		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_lindamp_m		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_specdamp		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_stir			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_robdamp		(is:ie, js:je, num_time_levels))

! Spectral Enstrophy Equation
allocate (Dyn%Grid(1,1)%enstrophy_tendk			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_pvor_adveck		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_rvormean_adveck	(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_trunck		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_lindamp_ek		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_lindamp_mk		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_specdampk		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_robdampk		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_stirk			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprime_yk			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%v_yk				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vorprimek_curr			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%vk				(is:ie, js:je, num_time_levels))

! Flux Term 
allocate (Dyn%Grid(1,1)%vorprime_fluxk			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%enstrophy_rvorprime_adveck_SUM(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_higher_freq			(is:ie, js:je, num_time_levels))
!I
allocate (Dyn%Grid(1,1)%I_pp				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_sp				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_ps				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_hp				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_ph				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_ss				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_hh				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_sh				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_hs				(is:ie, js:je, num_time_levels))
!I2
allocate (Dyn%Grid(1,1)%I2_pp				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_sp				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_ps				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_hp				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_ph				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_ss				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_hh				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_sh				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I2_hs				(is:ie, js:je, num_time_levels))
! I OTHER
allocate (Dyn%Grid(1,1)%I_0p				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_p0				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_0s				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_s0				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_0h				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_h0				(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%I_00				(is:ie, js:je, num_time_levels))

allocate(Grid(ie/2,ie/2))

do p = 1,128
   allocate(Dyn%Grid(p,1)%spec_deriv_vq_real(is:ie, js:je, num_time_levels))
   allocate(Dyn%Grid(p,1)%spec_deriv_vq_imag(is:ie, js:je, num_time_levels))
end do

do p = 1,12! ie/2
 do q = 1,12! ie/2 
   allocate(Dyn%Grid(p,q)%I1(is:ie, js:je, num_time_levels))
   allocate(Dyn%Grid(p,q)%I2(is:ie, js:je, num_time_levels))
   allocate(Dyn%Grid(p,q)%I3(is:ie, js:je, num_time_levels))
   allocate(Dyn%Grid(p,q)%I4(is:ie, js:je, num_time_levels))
!   allocate (Dyn%Grid(1,q)%enstrophy_robdamp		(is:ie, js:je, num_time_levels))
 end do
end do

!allocate (Dyn%Testing%f_test			(is:ie, is:ie, js:je, num_time_levels))
!allocate (Dyn%Testing%enstrophy_rvorprime_adveck(is:ie, is:ie, js:je, num_time_levels))

! Joe Terms
allocate (Dyn%Grid(1,1)%vvor_beta		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%source		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%source_beta		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%sink			(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%sink_beta		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%robertssink		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%dens_dt		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%ensflux		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%ensflux_div		(is:ie, js:je, num_time_levels))
allocate (Dyn%Grid(1,1)%ensflux_div_beta	(is:ie, js:je, num_time_levels))

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allocate (div (ms:me, ns:ne))

 call fv_advection_init(num_lon, num_lat, glat_bnd, 360./float(fourier_inc))
if(Dyn%grid_tracer) then
  allocate(Dyn%Grid(1,1)%tr (is:ie, js:je, num_time_levels))
  allocate(Dyn%Tend%tr (is:ie, js:je))
endif

if(Dyn%spec_tracer) then
  allocate(Dyn%Grid(1,1)%trs (is:ie, js:je, num_time_levels))
  allocate(Dyn%Tend%trs (is:ie, js:je))
  allocate(Dyn%Spec%trs (ms:me, ns:ne, num_time_levels))
endif

if(Time == Time_init) then

if(read_file == 1) then
	open (unit=1, file='/scratch/z3410755/u.txt', status='old', action='read')
	read (1,*) icnt
	allocate(u_init(icnt))
		do i = 1, icnt
		read (1,*) u_init(i)
		end do
	  do j = js, je
	    Dyn%Grid(1,1)%u(:,j,1) =   mult*u_init(j)
	    Dyn%Grid(1,1)%v(:,j,1) = 0.0 ! Dyn%Grid(1,1)%v cannot be given any other value here (by mass continuity). If another value is given it is ignored. To load a zonally varying v profile (that still obeys 					    ! mass continuity) use read_file = 4 initialisation
	  end do
endif

if(read_file == 2) then
 call read_data('INPUT/u.nc', 'u',   Dyn%Grid(1,1)%u  (:,:,1), grid_domain, timelevel=2)
   do j = js, je
Dyn%Grid(1,1)%u  (:,j,1) = sum(Dyn%Grid(1,1)%u(:,j,1))/128
Dyn%Grid(1,1)%v(:,j,1) = 0.0
end do
endif


if(read_file == 3) then
 call read_data('INPUT/u.nc', 'ug',   Dyn%Grid(1,1)%u  (:,:,1), grid_domain, timelevel=2)
   do j = js, je
Dyn%Grid(1,1)%u  (:,j,1) = sum(Dyn%Grid(1,1)%u(:,j,1))/128
Dyn%Grid(1,1)%v(:,j,1) = 0.0
end do
endif


if(read_file == 4) then ! read full lon-lat u & v wind profile from .nc files
 call read_data('/srv/ccrc/data15/z3410755/init_data/u.nc', 'u',   Dyn%Grid(1,1)%u  (:,:,1), grid_domain, timelevel=1)
 call read_data('/srv/ccrc/data15/z3410755/init_data/v.nc', 'v',   Dyn%Grid(1,1)%v  (:,:,1), grid_domain, timelevel=1)
endif

  call vor_div_from_uv_grid(Dyn%Grid(1,1)%u(:,:,1), Dyn%Grid(1,1)%v(:,:,1), Dyn%Spec%vor(:,:,1), div) 

  call trans_spherical_to_grid(Dyn%Spec%vor(:,:,1), Dyn%Grid(1,1)%vor(:,:,1))

  do j = js, je
    do i = is, ie
      yy = (deg_lat(j)- eddy_lat)/eddy_width
      Dyn%Grid(1,1)%vor(i,j,1) = Dyn%Grid(1,1)%vor(i,j,1) + &
              0.5*zeta_0*cos_lat(j)*exp(-yy*yy)*cos(m_0*rad_lon(i))
    end do
  end do
 
 call trans_grid_to_spherical(Dyn%Grid(1,1)%vor(:,:,1), Dyn%Spec%vor(:,:,1))
  
  div = (0.,0.)
  call uv_grid_from_vor_div   (Dyn%Spec%vor(:,:,1), div,        &
                              Dyn%Grid(1,1)%u  (:,:,1), Dyn%Grid(1,1)%v  (:,:,1))

if(Dyn%grid_tracer) then
    Dyn%Grid(1,1)%tr = 0.0
    do j = js, je
      if(deg_lat(j) > 10.0 .and. deg_lat(j) < 20.0) Dyn%Grid(1,1)%tr(:,j,1) =  1.0
      if(deg_lat(j) > 70.0 )                        Dyn%Grid(1,1)%tr(:,j,1) = -1.0
    end do
  endif
  
  if(Dyn%spec_tracer) then
    Dyn%Grid(1,1)%trs = 0.0
    do j = js, je
      if(deg_lat(j) > 10.0 .and. deg_lat(j) < 20.0) Dyn%Grid(1,1)%trs(:,j,1) =  1.0
      if(deg_lat(j) > 70.0 )                        Dyn%Grid(1,1)%trs(:,j,1) = -1.0
    end do
    call trans_grid_to_spherical(Dyn%Grid(1,1)%trs(:,:,1), Dyn%Spec%trs(:,:,1))
  endif
  
else
!print *, "calling"
  call read_restart(Dyn)
endif

module_is_initialized = .true.

return
end subroutine barotropic_dynamics_init

!===============================================================================================

subroutine barotropic_dynamics(Time, Time_init, Dyn, previous, current, future, delta_t)

type(time_type)    , intent(in)     :: Time, Time_init
type(dynamics_type), intent(inout)  :: Dyn
integer, intent(in   )  :: previous, current, future
real,    intent(in   )  :: delta_t

! < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < >

complex, dimension(ms:me, ns:ne)  :: dt_vors, dt_vors_prev, damped, dt_divs, dt_divs_throwaway, stream, zeros, spec_diss
complex, dimension(ms:me, ns:ne)  :: vors_rub, divs_rub, dt_vors_rub, dt_divs_rub, rob_vors_diff
real,    dimension(is:ie, js:je)  :: dt_vorg, stir, uprime, vormean, vorprime, energyprev, energymean_prev
real,    dimension(is:ie, js:je)  :: uprev, vprev, vorprev, zerog, oneg, rob_v_diff, rob_u_diff, rob_vorg_diff, energytendcheck, energytendcheck_mean
real,    dimension(is:ie, js:je)  :: utendmean, vtendmean, vqmean, uqmean, dx_u2_v2mean, dy_u2_v2mean, utruncmean, vtruncmean, ulindamp_mmean, &
				     ulindamp_emean, vlindampmean, uspecdampmean, vspecdampmean, ustirmean, vstirmean, urobdampmean, vrobdampmean, dt_vorg_rub
real,    dimension(is:ie, js:je)  :: enstrophytendcheck, vortendmean, pvor_advecmean, rvor_advecmean, vortruncmean, vorlindamp_emean, vorlindamp_mmean, vorspecdampmean, vorstirmean, vorrobdampmean, &
				     vortendprime, pvor_advecprime, rvor_advecprime, vortruncprime, vorlindamp_eprime, vorlindamp_mprime, vorspecdampprime, vorstirprime, vorrobdampprime, vprime, &
				     enstrophytrunc_mean, enstrophytrunc_prime, enstrophystir_mean, enstrophystir_prime, vorprimetend, rvor_advec2, tester, tester2, tester3
real,    dimension(is:ie, js:je)  :: umean, vmean, umean_curr, vmean_curr,  umean_prev, vmean_prev, vormean_prev, vorprime_prev, enstrophyprev, vorprime_curr, vormean_curr, vprime_curr, pvmean_curr
real,    dimension(is:ie, js:je)  :: f_array, pvmean_advec, beta_star, vormean_y_prev
real,    dimension(js:je-2)	  :: dlat
complex, dimension(is:ie, js:je):: uk, vk, v_yk, vork, vork_curr, vork_prev, vorprimetendk, vorprimek, vorprimek_curr, vorprimek_prev, pvor_adveck, rvormean_adveck, enstrophy_pvor_adveck, & 
				     enstrophy_rvormean_adveck, vorprime_xk, vorprime_yk, vor_xk, vor_yk, vorprime_fluxk
complex, dimension(is:ie, is:ie, js:je):: vor_flux_spec
complex, dimension(is:ie, js:je)  :: I_pp, I_sp, I_ps, I_hp, I_ph, I_ss, I_hh, I_sh, I_hs, I_0p, I_p0, I_0s, I_s0, I_0h, I_h0, I_00
complex, dimension(is:ie, js:je)  :: I2_pp, I2_sp, I2_ps, I2_hp, I2_ph, I2_ss, I2_hh, I2_sh, I2_hs
complex, dimension(js:je+1, is:ie) :: vorprimel_curr, rubl
complex, dimension(is:ie+1, js:je):: vortendk, enstrophy_tendk, vbetak, vdy_vormeank, vortrunck, vorlindamp_ek, vorlindamp_mk, vorspecdampk, vorrobdampk, vorstirk, vor_xkII, dvork, vor_ykII, enstrophy_trunck, enstrophy_lindamp_ek, enstrophy_lindamp_mk, enstrophy_specdampk, enstrophy_robdampk, enstrophy_stirk
real, dimension(is:ie, is:ie, js:je) :: spec_deriv_vq_real, spec_deriv_vq_imag
complex, dimension(is:ie, is:ie, is:ie, js:je) :: enstrophy_rvorprime_adveck, enstrophy_rvorprime_adveck2
integer, dimension(is:ie+1) :: kmat
integer, dimension(js:je) :: lmat

integer :: j
integer :: p, q, k
 character(3) :: str1, str2
 character(6) :: str3, str4
! character(1024) :: filename

! < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < >

if(.not.module_is_initialized) then
  call error_mesg('barotropic_dynamics','dynamics has not been initialized ', FATAL)
endif

zeros = (0.,0.) ! Zero array with spectral dimensions
zerog = (0.,0.) ! Zero array with grid dimensions
oneg = (0.,0.) ! Ones array with grid dimensions

!Dyn%testing(1,1)%I = 1
!Dyn%testing(1,12)%I = 2
!print *, Dyn%testing(1,1)%I(1:5,1,:)
!print *, Dyn%testing(1,2)%I(1:5,1,:)
!print *, Dyn%testing(1,12)%I(1:5,1,:)
!print *, size(Dyn%testing(1,1)%I,1), size(Dyn%testing(1,1)%I,2), size(Dyn%testing(1,1)%I,3)

! ---------------------------------- PREVIOUS & CURRENT TERMS ---------------------------------------
! PREVIOUS
uprev = Dyn%Grid(1,1)%u(:,:,previous) ! Save these variable (u @ t = i-1) in order to calculate tendancies later on (need to specify like this because the 212 time scheme writes over it at step 6)
vprev = Dyn%Grid(1,1)%v(:,:,previous)
vorprev = Dyn%Grid(1,1)%vor(:,:,previous)
energyprev = Dyn%Grid(1,1)%energy(:,:,previous) ! Havn't written this to a restart (don't use it to calculate the tendency). Just here for check.
energymean_prev = Dyn%Grid(1,1)%energy_mean(:,:,previous) 
vormean_y_prev = Dyn%Grid(1,1)%vormean_y(:,:,previous) 

do j = js, je
 umean_prev(:,j)= sum(uprev(:,j))/count(uprev(:,j) > -10000) 
 vmean_prev(:,j)= sum(vprev(:,j))/count(vprev(:,j) > -10000)
 vormean_prev(:,j)= sum(Dyn%Grid(1,1)%vor(:,j,previous))/count(Dyn%Grid(1,1)%u(:,j,previous) > -10000)
end do
vorprime_prev = Dyn%Grid(1,1)%vor(:,:,previous)-vormean_prev(:,:)
enstrophyprev = vorprime_prev*vorprime_prev
vork_prev(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vor(:,:,previous))
vork_prev(ie/2+2:ie+1,:) = conjg(vork_prev(ie/2+1:2:-1,:))
! CURRENT
do j = js, je
  umean_curr(:,j)= sum(Dyn%Grid(1,1)%u(:,j,current))/count(Dyn%Grid(1,1)%u(:,j,current) > -10000) ! Calculate zonal mean, zonal flow (needed for linear damping)
  vmean_curr(:,j)= sum(Dyn%Grid(1,1)%v(:,j,current))/count(Dyn%Grid(1,1)%v(:,j,future) > -10000)
  vormean_curr(:,j)= sum(Dyn%Grid(1,1)%vor(:,j,current))/count(Dyn%Grid(1,1)%u(:,j,current) > -10000) ! Calculate zonal mean, vorticity (needed for calculating beta_star)
  pvmean_curr(:,j)= sum(Dyn%Grid(1,1)%vor(:,j,current))/count(Dyn%Grid(1,1)%u(:,j,current) > -10000) + coriolis(j) ! Calculate zonal mean, pv (needed for calculating beta_star)
end do
vprime_curr = Dyn%Grid(1,1)%v(:,:,current)-vmean_curr(:,:)
vorprime_prev = Dyn%Grid(1,1)%vor(:,:,previous)-vormean_prev(:,:)
vorprime_curr = Dyn%Grid(1,1)%vor(:,:,current)-vormean_curr(:,:)

vork_prev(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vor(:,:,previous))
vork_prev(ie/2+2:ie+1,:) = conjg(vork_prev(ie/2+1:2:-1,:))
vork_curr(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vor(:,:,current))
vork_curr(ie/2+2:ie+1,:) = conjg(vork_curr(ie/2+1:2:-1,:))

vorprimek_prev(1:ie/2+1,:) = fft_grid_to_fourier(vorprime_prev)
vorprimek_prev(ie/2+2:ie+1,:) = conjg(vorprimek_prev(ie/2+1:2:-1,:))
vorprimek_curr(1:ie/2+1,:) = fft_grid_to_fourier(vorprime_curr)
vorprimek_curr(ie/2+2:ie+1,:) = conjg(vorprimek_curr(ie/2+1:2:-1,:))

Dyn%Grid(1,1)%vorprime_prev(:,:,future) = vorprime_prev
Dyn%Grid(1,1)%vorprime_curr(:,:,future) = vorprime_curr
! ====================== FFT along LATITUDE??? ==================
!call fft_end
!call fft_init(size(vorprime_curr,2))
!
!print *, size(transpose(vorprime_curr(1:10,:)),1)!, size(transpose(vorprime_curr(1,:)),2)
!print *, size(fft_grid_to_fourier(transpose(vorprime_curr(1:10,:))),1)
!print *, size(transpose(vorprime_curr),1), size(transpose(vorprime_curr),2)

! call fft_end
! call fft_init(size(vorprime_curr,1))

!vorprimel_curr(:,1:je/2+1) = fft_grid_to_fourier(transpose(vorprime_curr)) ! FFT taken along lat (use for enstrophy_rvorprime_adveck)
!vorprimel_curr(:,je/2+2:je+1) = conjg(vorprimel_curr(:,je/2+1:2:-1))
! ===============================================================

Dyn%Grid(1,1)%uprev(:,:,future) = Dyn%Grid(1,1)%v(:,:,previous)
Dyn%Grid(1,1)%ucurr(:,:,future) = Dyn%Grid(1,1)%u(:,:,current)
Dyn%Grid(1,1)%vprev(:,:,future) = Dyn%Grid(1,1)%v(:,:,previous)
Dyn%Grid(1,1)%vcurr(:,:,future) = Dyn%Grid(1,1)%v(:,:,current)
Dyn%Grid(1,1)%vorprev(:,:,future) = Dyn%Grid(1,1)%vor(:,:,previous)
Dyn%Grid(1,1)%vorcurr(:,:,future) = Dyn%Grid(1,1)%vor(:,:,current)

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do j = js, je
  Dyn%Grid(1,1)%pv(:,j)  = Dyn%Grid(1,1)%vor(:,j,current)  + coriolis(j) 
  f_array(:,j) = coriolis(j) ! needed below in advection terms calc
end do

!Dyn%Testing%f_test(:,:,:,future) = 1.0
!print *, size(Dyn%Testing%f_test,1), size(Dyn%Testing%f_test,2), size(Dyn%Testing%f_test,3), size(Dyn%Testing%f_test,4)
!print *,  "F_TEST", Dyn%Testing%f_test(22,22,22,future)

Dyn%Tend%u = Dyn%Tend%u + Dyn%Grid(1,1)%pv*Dyn%Grid(1,1)%v(:,:,current)
Dyn%Tend%v = Dyn%Tend%v - Dyn%Grid(1,1)%pv*Dyn%Grid(1,1)%u(:,:,current)

! ---------------------------------- ADVECTION TERMS ---------------------------------------
Dyn%Grid(1,1)%vq(:,:,future) = Dyn%Grid(1,1)%pv*Dyn%Grid(1,1)%v(:,:,current) ! Full advection terms
Dyn%Grid(1,1)%uq(:,:,future) = Dyn%Grid(1,1)%pv*Dyn%Grid(1,1)%u(:,:,current)! Note the positive sign
!print *, "VQ", Dyn%Grid(1,1)%vq(10,10,future)
!print *, "-UQ", -Dyn%Grid(1,1)%uq(10,10,future)

! RELATIVE VORTICITY ADVECTION BREAK INTO ZETA_BAR and ZETA_PRIME COMPONENT
! ZETA_BAR
!print *, "VORMEAN_CURR", vormean_curr(22,47)
!print *, "VORPRIME_CURR", vorprime_curr(22,47)

 call vor_div_from_uv_grid (vormean_curr*Dyn%Grid(1,1)%v(:,:,current), -vormean_curr*Dyn%Grid(1,1)%u(:,:,current), dt_vors_rub, dt_divs_rub)  
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%rvormean_advec(:,:,future))
!print *, 'RVORMEAN_ADVEC', Dyn%Grid(1,1)%rvormean_advec(22,47,future)
! ZETA_PRIME
 call vor_div_from_uv_grid (vorprime_curr*Dyn%Grid(1,1)%v(:,:,current), -vorprime_curr*Dyn%Grid(1,1)%u(:,:,current), dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%rvorprime_advec(:,:,future))
!print *, 'RVORPRIME_ADVEC', Dyn%Grid(1,1)%rvorprime_advec(22,47,future)
 call vor_div_from_uv_grid (Dyn%Grid(1,1)%vor(:,:,current)*Dyn%Grid(1,1)%v(:,:,current), -Dyn%Grid(1,1)%vor(:,:,current)*Dyn%Grid(1,1)%u(:,:,current), dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%rvor_advec(:,:,future))
!print *, "RVOR_ADVEC", Dyn%Grid(1,1)%rvor_advec(22,47,future)
!print *, "TOTAL", Dyn%Grid(1,1)%rvormean_advec(22,47,future) + Dyn%Grid(1,1)%rvorprime_advec(22,47,future)
! PLANETARY VORTICITY ADVECTION
 call vor_div_from_uv_grid (f_array*Dyn%Grid(1,1)%v(:,:,current), -f_array*Dyn%Grid(1,1)%u(:,:,current), dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%pvor_advec(:,:,future))
!print *, 'PVOR_ADVEC' ,Dyn%Grid(1,1)%pvor_advec(22,47,future)
! BETA
Dyn%Grid(1,1)%beta(:,:,future) = -Dyn%Grid(1,1)%pvor_advec(:,:,future)/Dyn%Grid(1,1)%v(:,:,current)!print *, "V", Dyn%Grid(1,1)%v(1,:,current)
!print *, "Beta", Dyn%Grid(1,1)%beta(1,:,future)
! ZETA_BAR_X
 call vor_div_from_uv_grid (zerog, vormean_curr, dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vormean_x(:,:,future))
! ZETA_BAR_Y
 call vor_div_from_uv_grid (-vormean_curr, zerog, dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vormean_y(:,:,future))
! ZETA_PRIME_X
 call vor_div_from_uv_grid (zerog, vorprime_curr, dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vorprime_x(:,:,future))
! ZETA_PRIME_Y
 call vor_div_from_uv_grid (-vorprime_curr, zerog, dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vorprime_y(:,:,future))
! ZETA_X
 call vor_div_from_uv_grid (zerog, Dyn%Grid(1,1)%vor(:,:,current), dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vor_x(:,:,future))
! ZETA_Y
 call vor_div_from_uv_grid (-Dyn%Grid(1,1)%vor(:,:,current), zerog, dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vor_y(:,:,future))
! V_Y
 call vor_div_from_uv_grid (-Dyn%Grid(1,1)%u(:,:,current), zerog, dt_vors_rub, dt_divs_rub)  ! CHANGE BACK TO V WHEN FINISHED!!!!
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%v_y(:,:,future))
!print *,  'VORMEAN', vormean(22,:)
!print *,  'VORMEAN_Y', Dyn%Grid(1,1)%vormean_y(22,:,future)
! PVMEAN_ADVEC & BETA_STAR
Dyn%Grid(1,1)%beta_star(:,:,future) = Dyn%Grid(1,1)%beta(:,:,future) + Dyn%Grid(1,1)%vormean_y(:,:,future)
!print *, "BETA_STAR",  Dyn%Grid(1,1)%beta_star(22,47,future)

!========= COMMENT ON vor_div_from_uv_grid SUBROUTINE =============
! I have found that applying the product rule to this derivative operator does not give an exact result. That is to say
!			 (AB)_y =/= A(B)_y + (A)_yB
! Proof of this (using v*q as the fields of interest)
!LHS
 call vor_div_from_uv_grid (Dyn%Grid(1,1)%vor(:,:,current)*Dyn%Grid(1,1)%v(:,:,current), zerog, dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%tester(:,:,future))
! 1st TERM RHS
! call vor_div_from_uv_grid (Dyn%Grid(1,1)%vor(:,:,current), zerog, dt_vors_rub, dt_divs_rub) 
! call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%tester2(:,:,future))
! Dyn%Grid(1,1)%tester2(:,:,future) = Dyn%Grid(1,1)%v(:,:,current)*Dyn%Grid(1,1)%tester2(:,:,future)
! 2nd TERM RHS
! call vor_div_from_uv_grid (Dyn%Grid(1,1)%v(:,:,current), zerog, dt_vors_rub, dt_divs_rub) 
! call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%tester3(:,:,future))
! Dyn%Grid(1,1)%tester3(:,:,future) = Dyn%Grid(1,1)%vor(:,:,current)*Dyn%Grid(1,1)%tester3(:,:,future)
! Mathematically tester = tester2 + tester3
!
! Now I want to see if I can can calculate this derivative exactly by breaking the terms in spectral space and recombining.
! BUILD UP SPECTRAL ARRAYS; Has (almost) the same structure as Matlab code ie wavenumber takes the form [0 1 2 ... N/2 -N/2 -N/2+1 ... -1]
uk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%u(:,:,current))
uk(ie/2+1:ie,:) = conjg(uk(ie/2+1:2:-1,:))
vk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%v(:,:,current))
vk(ie/2+1:ie,:) = conjg(vk(ie/2+1:2:-1,:))
vork(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vor(:,:,current))
vork(ie/2+1:ie,:) = conjg(vork(ie/2+1:2:-1,:))

do p=1,size(kmat)
       do q=1,size(kmat)
	   vor_flux_spec(p,q,:) = vk(p,:)*vork(q,:)      
	end do
end do

!print *, size(kmat)
!print *, size(vk,1), size(vk,2)
!print *, size(vor_flux_spec(:,[1:ie/2 ie/2+2:ie],:),1), size(vor_flux_spec(:,[1:ie/2 ie/2+2:ie],:),2), size(vor_flux_spec(:,[1:ie/2 ie/2+2:ie],:),3)

do p=1,size(kmat) 
!	 print *, p
	 call vor_div_from_uv_grid (real(vor_flux_spec(p,1:ie,:)), zerog, dt_vors_rub, dt_divs_rub) 
         call trans_spherical_to_grid(dt_vors_rub, dt_vorg_rub)
!	 print *, dt_vorg_rub
	 Dyn%Grid(p,1)%spec_deriv_vq_real(:,:,future) = dt_vorg_rub
	 call vor_div_from_uv_grid (aimag(vor_flux_spec(p,1:ie,:)), zerog, dt_vors_rub, dt_divs_rub) 
         call trans_spherical_to_grid(dt_vors_rub, dt_vorg_rub)
	 Dyn%Grid(p,1)%spec_deriv_vq_imag(:,:,future) = dt_vorg_rub
! TEST ON U
	 call vor_div_from_uv_grid (real(uk), zerog, dt_vors_rub, dt_divs_rub) 
         call trans_spherical_to_grid(dt_vors_rub, dt_vorg_rub)
	 Dyn%Grid(1,1)%tester2(:,:,future) = dt_vorg_rub
	 call vor_div_from_uv_grid (aimag(uk), zerog, dt_vors_rub, dt_divs_rub) 
         call trans_spherical_to_grid(dt_vors_rub, dt_vorg_rub)
	 Dyn%Grid(1,1)%tester3(:,:,future) = dt_vorg_rub
end do

!print *, Dyn%Grid(13,1)%spec_deriv_vq_real(2,:,future)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!========= COMMENT ON horizontal_advection SUBROUTINE =============
! Both methods give (almost) the same result. I think the difference can be attributed to some numerical errors in the different ways they are calculated  (small => not important)
! METHOD 1
! call trans_grid_to_spherical(vorprime_curr, vors_rub)
! call horizontal_advection   (vors_rub, 0*Dyn%Grid(1,1)%u(:,:,current), Dyn%Grid(1,1)%v(:,:,current), Dyn%Grid(1,1)%tester(:,:,future))
! Dyn%Grid(1,1)%tester(:,:,future) = Dyn%Grid(1,1)%tester(:,:,future)/Dyn%Grid(1,1)%v(:,:,current)
!print *, Dyn%Grid(1,1)%tester(:,:,future)
! call trans_grid_to_spherical(f_array, vors_rub)
! call horizontal_advection   (vors_rub, 0*Dyn%Grid(1,1)%u(:,:,current), Dyn%Grid(1,1)%v(:,:,current), Dyn%Grid(1,1)%beta(:,:,future))
!Dyn%Grid(1,1)%beta(:,:,future) = Dyn%Grid(1,1)%beta(:,:,future)/Dyn%Grid(1,1)%v(:,:,current)
!print *, "METHOD 1: U.GRAD(f)", Dyn%Grid(1,1)%tester(10,10,future)
! METHOD 2
! call vor_div_from_uv_grid (Dyn%Grid(1,1)%vor(:,:,current)*Dyn%Grid(1,1)%v(:,:,current), -Dyn%Grid(1,1)%vor(:,:,current)*Dyn%Grid(1,1)%u(:,:,current), dt_vors_rub, dt_divs_rub) 
! call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%rvor_advec(:,:,future))
!print *, "METHOD 2: U.GRAD(VOR)", Dyn%Grid(1,1)%rvor_advec(10,10,future)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 call vor_div_from_uv_grid (Dyn%Tend%u, Dyn%Tend%v, dt_vors, dt_divs) ! Calculate DT_VOR before linear damping is applied
 call trans_spherical_to_grid(dt_vors, Dyn%Grid(1,1)%dt_vor(:,:,future))
!print *, "DT_VOR", Dyn%Grid(1,1)%dt_vor(22,47,future) 
!print *, "DT_VOR - XXX ...", Dyn%Grid(1,1)%dt_vor(22,47,future) - (Dyn%Grid(1,1)%rvormean_advec(22,47,future) + Dyn%Grid(1,1)%rvorprime_advec(22,47,future) + Dyn%Grid(1,1)%pvor_advec(22,47,future))
! Further subdivide the advection term of q=(zeta + f) into a vorticity advection component and a beta*v component
!vor_beta = beta*Dyn%Grid(1,1)%v(:,:,current)
!vor_advec = Dyn%Grid(1,1)%dt_vor(:,:,future) - vor_beta

!===================== COMMENT ON DT_VOR ===========================
! The sum of planetary vorticity advection (pvor_advec) and relative vorticity advection (rvor_advec) gives the absolute vorticity advection given in dt_vor
!print *, "ABSOLUTE ADVEC", Dyn%Grid(1,1)%pvor_advec(22,47,future) + Dyn%Grid(1,1)%rvor_advec(22,47,future)
!print *, "(PVOR + RVOR) - ABSOLUTE ADVEC", (Dyn%Grid(1,1)%pvor_advec(22,47,future) + Dyn%Grid(1,1)%rvor_advec(22,47,future)) - Dyn%Grid(1,1)%dt_vor(22,47,future)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ---------  JOE CODE (SNIPPET FOR QUASILINEAR and LINEAR EQUATIONS) ----------------------
! for QL integrations:
!Dyn%Tend%u=Dyn%Tend%u+linmult*Dyn%Grid(1,1)%dt_vor(:,:,future)
!make a linear model - product of eddies removed.
!Dyn%Tend%u = Dyn%Tend%u + Dyn%Grid(1,1)%pv*Dyn%Grid(1,1)%v(:,:,current)
!Dyn%Tend%v = Dyn%Tend%v - Dyn%Grid(1,1)%pv*(Dyn%Grid(1,1)%u(:,:,current)-umean)
!do j = js, je
!  Dyn%Grid(1,1)%pv(:,j)  = Dyn%Grid(1,1)%vor(:,j,current)  + coriolis(j) 
!end do
! add total vorticity times zonal mean u, so still missing products of eddies (u'zeta')
!Dyn%Tend%v = Dyn%Tend%v - Dyn%Grid(1,1)%pv*umean
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ----------------------------- LINEAR DAMPING ---------------------------------------------
if(damp > 1) then ! do linear damping
 Dyn%Grid(1,1)%vlindamp(:,:,future) =  Dyn%Grid(1,1)%v(:,:,current)/tau
 Dyn%Grid(1,1)%ulindamp_e(:,:,future) = (Dyn%Grid(1,1)%u(:,:,current)-umean_curr)/tau
 Dyn%Grid(1,1)%ulindamp_m(:,:,future) = umean_curr/tau1
 !print *, "V (EDDY) LINDAMP", Dyn%Grid(1,1)%vlindamp(10,10,future)
 !print *, "U EDDY LINDAMP", Dyn%Grid(1,1)%ulindamp_e(10,10,future)
 !print *, "U MEAN LINDAMP", Dyn%Grid(1,1)%ulindamp_m(10,10,future)
 Dyn%Tend%v = Dyn%Tend%v - Dyn%Grid(1,1)%vlindamp(:,:,future)  ! Damping on merdional (eddy) flow
 Dyn%Tend%u = Dyn%Tend%u - Dyn%Grid(1,1)%ulindamp_e(:,:,future) ! Damping on zonal eddy flow
 Dyn%Tend%u = Dyn%Tend%u - Dyn%Grid(1,1)%ulindamp_m(:,:,future) ! Damping on zonal mean flow
 
 call vor_div_from_uv_grid (Dyn%Grid(1,1)%ulindamp_e(:,:,future), Dyn%Grid(1,1)%vlindamp(:,:,future), dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vorlindamp_e(:,:,future)) ! Calculate the effect of eddy damping of u and v on vorticity
 call vor_div_from_uv_grid (Dyn%Grid(1,1)%ulindamp_m(:,:,future), zerog, dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vorlindamp_m(:,:,future)) ! Calculate the effect of mean flow damping of u on vorticity
 !print *, "VOR EDDY LINDAMP", Dyn%Grid(1,1)%vorlindamp_e(10,10,future) 
 !print *, "VOR MEAN LINDAMP", Dyn%Grid(1,1)%vorlindamp_m(10,10,future) 
endif
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! --------------------------------- TRUNCATION TERMS ---------------------------------------
! U & V
 call trans_grid_to_spherical(Dyn%Tend%u, dt_vors_rub)
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%utrunc(:,:,current))
 call trans_grid_to_spherical(Dyn%Tend%v, dt_vors_rub)
 call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vtrunc(:,:,current))
Dyn%Grid(1,1)%utrunc(:,:,future) = Dyn%Tend%u - Dyn%Grid(1,1)%utrunc(:,:,current)
Dyn%Grid(1,1)%vtrunc(:,:,future) = Dyn%Tend%v - Dyn%Grid(1,1)%vtrunc(:,:,current)
!print *, "UTRUNC", Dyn%Grid(1,1)%utrunc(10,10,future)
!print *, "VTRUNC", Dyn%Grid(1,1)%vtrunc(10,10,future)

!============= COMMENT ON U^2 TRUNCATION TERM =====================
! call trans_grid_to_spherical(Dyn%Tend%u**2, dt_vors_rub)
! call trans_spherical_to_grid(dt_vors_rub, u2trunc)
!print *, "U^2 TRUNC", u2trunc(40,40)
! This truncation TRUNC(uvq) does not work for the dt_u^2 budget. Instead u*TRUNC(u) works. On the one hand this makes sense, the dt_u^2 equation is just (up to a factor of 0.5) u*dt_u
!	u*dt_u = 0.5*dt_u^2 = u*vq - u*TRUNC
!On the other hand I would have thought TRUNC(uvq) should also have worked. It would be good to understand why it does not.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! VORTICITY
if(previous == current) then
 call trans_grid_to_spherical(Dyn%Grid(1,1)%vor(:,:,current), vors_rub)
 call trans_spherical_to_grid(vors_rub, Dyn%Grid(1,1)%vortrunc(:,:,future))
Dyn%Grid(1,1)%vortrunc(:,:,future) = (Dyn%Grid(1,1)%vor(:,:,current)-Dyn%Grid(1,1)%vortrunc(:,:,future))/delta_t
!print *, "VORTRUNC", Dyn%Grid(1,1)%vortrunc(10,10,future)
else
 call trans_grid_to_spherical(Dyn%Grid(1,1)%vor(:,:,previous), vors_rub)
 call trans_spherical_to_grid(vors_rub, Dyn%Grid(1,1)%vortrunc(:,:,future))
Dyn%Grid(1,1)%vortrunc(:,:,future) = (Dyn%Grid(1,1)%vor(:,:,previous)-Dyn%Grid(1,1)%vortrunc(:,:,future))/delta_t
!print *, "VORTRUNC", Dyn%Grid(1,1)%vortrunc(10,10,future)
endif
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -------------------------------------- GRAD(U^2 + V^2) -----------------------------------------
 call vor_div_from_uv_grid (Dyn%Tend%u, Dyn%Tend%v, dt_vors, dt_divs) ! 3. Moving into spectral space, recast equ. of motion (EoF) in vor-div; compute spectral divergence of {u_t,v_t}={(f+zeta)u,(f+zeta)v} 
 call uv_grid_from_vor_div(dt_vors,  zeros, Dyn%Grid(1,1)%dx_u2_v2(:,:,current), Dyn%Grid(1,1)%dy_u2_v2(:,:,current))
Dyn%Grid(1,1)%dx_u2_v2(:,:,future) = Dyn%Tend%u - Dyn%Grid(1,1)%dx_u2_v2(:,:,current) - Dyn%Grid(1,1)%utrunc(:,:,future)
Dyn%Grid(1,1)%dy_u2_v2(:,:,future) = Dyn%Tend%v - Dyn%Grid(1,1)%dy_u2_v2(:,:,current) - Dyn%Grid(1,1)%vtrunc(:,:,future)
!print *, "DX_U2_V2", Dyn%Grid(1,1)%dx_u2_v2(10,10,future)
!print *, "DY_U2_V2", Dyn%Grid(1,1)%dy_u2_v2(10,10,future)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -------------------  JOE CODE (SNIPPET FOR LINEARISED EQUATIONS)-------------------------
!pvmean, for linearized equations:
!do j = js, je
!  Dyn%Grid(1,1)%pv(:,j)  = sum(Dyn%Grid(1,1)%vor(:,j,current))/count(Dyn%Grid(1,1)%u(:,j,current) > -10000)   + coriolis(j) 
!end do
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -------------------------------- SPECTRAL DAMPING ----------------------------------------
if(damp > 0) then
  dt_vors_rub = dt_vors
  call compute_spectral_damping(Dyn%Spec%vor(:,:,previous), dt_vors, delta_t) ! 4. Compute spectral damping, the subroutine (located in src/atmos_spectral/spectral_damping.f90 outputs dt_vors)
  dt_vors_rub = dt_vors - dt_vors_rub
  call uv_grid_from_vor_div(dt_vors_rub,  zeros, Dyn%Grid(1,1)%uspecdamp(:,:,future),  Dyn%Grid(1,1)%vspecdamp(:,:,future))
  call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%vorspecdamp(:,:,future))
endif
!print *, "USPECDAMP", Dyn%Grid(1,1)%uspecdamp(10,10,future)
!print *, "VSPECDAMP", Dyn%Grid(1,1)%vspecdamp(10,10,future)
!print *, "VORSPECDAMP", Dyn%Grid(1,1)%vorspecdamp(10,10,future)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ------------------------------- STIRRING -------------------------------------------------
dt_vors_rub = dt_vors
 call stirring(Time, dt_vors)
dt_vors_rub = dt_vors - dt_vors_rub
 call uv_grid_from_vor_div(dt_vors_rub,  zeros, Dyn%Grid(1,1)%ustir(:,:,future),  Dyn%Grid(1,1)%vstir(:,:,future))
 call trans_spherical_to_grid(dt_vors_rub,Dyn%Grid(1,1)%vorstir(:,:,future))
!print *, "U STIR", Dyn%Grid(1,1)%ustir(10,10,future)
!print *, "V STIR", Dyn%Grid(1,1)%vstir(10,10,future)
!print *, "VOR STIR", Dyn%Grid(1,1)%vorstir(10,10,future)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -------------------  JOE CODE (SNIPPET FOR ENSTROPHY EQUATION)----------------------------
!calc the enstrophy flux in grid space, take the divergence later.
!Dyn%Grid(1,1)%vtend(:,:,current)=Dyn%Grid(1,1)%vor(:,:,current)*Dyn%Grid(1,1)%v(:,:,current)
!
!do j = js, je
!vormean(:,j)= sum(Dyn%Grid(1,1)%vor(:,j,current))/count(Dyn%Grid(1,1)%u(:,j,current) > -10000) 
!end do
!vorprime = Dyn%Grid(1,1)%vor(:,:,current)-vormean(:,:)
!Dyn%Grid(1,1)%source(:,:,current) = stir*vorprime
!Dyn%Grid(1,1)%sink(:,:,current) = Dyn%Grid(1,1)%vortend(:,:,current)*vorprime
!Dyn%Grid(1,1)%ensflux(:,:,current) = 0.5*vorprime*vorprime*Dyn%Grid(1,1)%v(:,:,current)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ---------------------------------- ROBERTS DAMPING ---------------------------------------
Dyn%spec%roberts(:,:,future) = Dyn%spec%vor(:,:,current)
 call leapfrog(Dyn%Spec%vor , dt_vors  , previous, current, future, delta_t, robert_coeff) ! 5. Apply Roberts damping/filter using leapfrog scheme, again in spectral space
Dyn%spec%roberts(:,:,future) = Dyn%spec%vor(:,:,current)-Dyn%spec%roberts(:,:,future) ! This is the roberts damping term: rob_coeff*(phi_{n-1} - 2*phi_n + phi_{n+1})

!print *, "ROBO CURRENT", Dyn%spec%roberts(10,10,current)
rob_vors_diff = Dyn%spec%vor(:,:,future)-Dyn%spec%roberts(:,:,current) ! Gives the time-stepped spec%vor that results from using an unfiltered spec%vor(previous)
 call uv_grid_from_vor_div(rob_vors_diff,  zeros, rob_u_diff,  rob_v_diff)
 call trans_spherical_to_grid(rob_vors_diff, rob_vorg_diff)

 call trans_spherical_to_grid(Dyn%Spec%vor(:,:,future), Dyn%Grid(1,1)%vor(:,:,future)) ! 6. Transforms from spherical harmonics to grid space and then from vor-div to u-v EoF
 call uv_grid_from_vor_div(Dyn%Spec%vor(:,:,future),  zeros, Dyn%Grid(1,1)%u(:,:,future),  Dyn%Grid(1,1)%v(:,:,future))

Dyn%Grid(1,1)%urobdamp(:,:,future) = (Dyn%Grid(1,1)%u(:,:,future)-rob_u_diff(:,:))/delta_t
Dyn%Grid(1,1)%vrobdamp(:,:,future) = (Dyn%Grid(1,1)%v(:,:,future)-rob_v_diff(:,:))/delta_t
Dyn%Grid(1,1)%vorrobdamp(:,:,future) = (Dyn%Grid(1,1)%vor(:,:,future)-rob_vorg_diff(:,:))/delta_t
!print *, "UROBDAMP", Dyn%Grid(1,1)%urobdamp(10,10,future)
!print *, "VROBDAMP", Dyn%Grid(1,1)%vrobdamp(10,10,future)
!print *, "VORROBDAMP", Dyn%Grid(1,1)%vorrobdamp(10,10,future)

!======= COMMENT ON CALCULATION OF ROBERTS FILTER TERM =============
! I have written an extensive comment on the roberts filter operation and calculation, see src/atmos_spectral/model/Roberts_damping.pdf
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! AFTER THIS POINT THERE IS NO MORE "ORGINAL" CODE (BESIDES A SMALL SNIPPET FOR THE CALCULATION OF THE STREAMFUNCTION AND TRACER ADVECTION). THE REST IS CALCULATIONS OF TENDENCIES, OTHER FIELDS THAT ARE 
! DERIVED FROM U,V AND VOR FIELDS (EG. ENERGY & ENSTROPHY) AND BUDGETS FOR ALL THESE VARIABLES.
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! -----------------  JOE CODE (EFFECT OF KILLING CERTAIN WAVENUMBERS)------------------------
! try killing the affect of waves 1:3
!dt_vors(1:3,:)=0;
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -----------------  JOE CODE (SNIPPET FOR ENSTOPHY EQUATION)-------------------------------
!Dyn%Grid(1,1)%robertssink(:,:,current)=roberts*vorprime
!Dyn%Grid(1,1)%dens_dt(:,:,current)=dt_vorg - Dyn%Grid(1,1)%vor(:,:,current)
!
! call trans_grid_to_spherical  (vorprime, dt_divs_throwaway)
! call horizontal_advection     (dt_divs_throwaway, Dyn%Grid(1,1)%u(:,:,current)-umean, Dyn%Grid(1,1)%v(:,:,current), Dyn%Grid(1,1)%ensflux_div(:,:,current))
!Dyn%Grid(1,1)%ensflux_div(:,:,current) = Dyn%Grid(1,1)%ensflux_div(:,:,current)*vorprime
! 
!Zonal mean absolute vorticity
!do j = js, je
!  vormean(:,j)  = sum(Dyn%Grid(1,1)%vor(:,j,current))/count(Dyn%Grid(1,1)%u(:,j,current) > -10000)   + coriolis(j) 
!end do
!Zonal mean absolute vorticity in spherical, just for advection later
!call trans_grid_to_spherical  (vormean, dt_divs_throwaway)
!get v'beta*, on the way to getting vvor'beta*
!call horizontal_advection     (dt_divs_throwaway, 0*Dyn%Grid(1,1)%u(:,:,current), Dyn%Grid(1,1)%v(:,:,current), Dyn%Grid(1,1)%vvor_beta(:,:,current))
!get beta, to use for calculating source and sink of vorticity flux budget
!beta = Dyn%Grid(1,1)%vvor_beta(:,:,current)/Dyn%Grid(1,1)%v(:,:,current)
!Dyn%Grid(1,1)%vvor_beta(:,:,current) = Dyn%Grid(1,1)%vvor_beta(:,:,current)*vorprime
!vvor_beta
!Dyn%Grid(1,1)%ensflux_div_beta(:,:,current) = Dyn%Grid(1,1)%ensflux_div(:,:,current)/beta
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ------------------------------ TEST VARIABLE ----------------------------
!Nx = 36
!do k=1,Nx
! test_length(k,:)=2*pi*k/Nx
!end do
!
!print *, size(test_length,1), size(test_length,2)
!print *, "f(x)", cos(test_length)
!print *, "F(k)", fft_grid_to_fourier(cos(test_length))
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ------------------------------ FUTURE VALUES -------------------------------------------
! Take the mean for momentum budget
do j = js, je 
 umean(:,j)= sum(Dyn%Grid(1,1)%u(:,j,future))/count(Dyn%Grid(1,1)%u(:,j,future) > -10000) 
 vmean(:,j)= sum(Dyn%Grid(1,1)%v(:,j,future))/count(Dyn%Grid(1,1)%v(:,j,future) > -10000)
 utendmean(:,j) = sum(Dyn%Grid(1,1)%utend(:,j,future))/count(Dyn%Grid(1,1)%utend(:,j,future) > -10000)
 vtendmean(:,j) = sum(Dyn%Grid(1,1)%vtend(:,j,future))/count(Dyn%Grid(1,1)%vtend(:,j,future) > -10000)
 vqmean(:,j) = sum(Dyn%Grid(1,1)%vq(:,j,future))/count(Dyn%Grid(1,1)%vq(:,j,future) > -10000)
 uqmean(:,j) = sum(Dyn%Grid(1,1)%uq(:,j,future))/count(Dyn%Grid(1,1)%uq(:,j,future) > -10000)
 dx_u2_v2mean(:,j) = sum(Dyn%Grid(1,1)%dx_u2_v2(:,j,future))/count(Dyn%Grid(1,1)%dx_u2_v2(:,j,future) > -10000) ! Term shoud be zero (d/dx)
 dy_u2_v2mean(:,j) = sum(Dyn%Grid(1,1)%dy_u2_v2(:,j,future))/count(Dyn%Grid(1,1)%dy_u2_v2(:,j,future) > -10000)
 utruncmean(:,j) = sum(Dyn%Grid(1,1)%utrunc(:,j,future))/count(Dyn%Grid(1,1)%utrunc(:,j,future) > -10000) ! Both truncation terms should have a zonal mean of ~ zero (composed entirely of high frequency "eddies")
 vtruncmean(:,j) = sum(Dyn%Grid(1,1)%vtrunc(:,j,future))/count(Dyn%Grid(1,1)%vtrunc(:,j,future) > -10000)
 ulindamp_mmean(:,j) = sum(Dyn%Grid(1,1)%ulindamp_m(:,j,future))/count(Dyn%Grid(1,1)%ulindamp_m(:,j,future) > -10000)
 ulindamp_emean(:,j) = sum(Dyn%Grid(1,1)%ulindamp_e(:,j,future))/count(Dyn%Grid(1,1)%ulindamp_e(:,j,future) > -10000) ! should be zero
 vlindampmean(:,j) = sum(Dyn%Grid(1,1)%vlindamp(:,j,future))/count(Dyn%Grid(1,1)%vlindamp(:,j,future) > -10000) ! should be zero
 uspecdampmean(:,j) = sum(Dyn%Grid(1,1)%uspecdamp(:,j,future))/count(Dyn%Grid(1,1)%uspecdamp(:,j,future) > -10000)
 vspecdampmean(:,j) = sum(Dyn%Grid(1,1)%vspecdamp(:,j,future))/count(Dyn%Grid(1,1)%vspecdamp(:,j,future) > -10000)
 ustirmean(:,j) = sum(Dyn%Grid(1,1)%ustir(:,j,future))/count(Dyn%Grid(1,1)%ustir(:,j,future) > -10000)
 vstirmean(:,j) = sum(Dyn%Grid(1,1)%vstir(:,j,future))/count(Dyn%Grid(1,1)%vstir(:,j,future) > -10000)
 urobdampmean(:,j) = sum(Dyn%Grid(1,1)%urobdamp(:,j,future))/count(Dyn%Grid(1,1)%urobdamp(:,j,future) > -10000)
 vrobdampmean(:,j) = sum(Dyn%Grid(1,1)%vrobdamp(:,j,future))/count(Dyn%Grid(1,1)%vrobdamp(:,j,future) > -10000)
end do
 vprime = Dyn%Grid(1,1)%v(:,:,future)-vmean(:,:)

! Find the eddy quantities to calculate vorticity budget
do j = js, je 
 vormean(:,j)= sum(Dyn%Grid(1,1)%vor(:,j,future))/count(Dyn%Grid(1,1)%u(:,j,future) > -10000) 
 vortendmean(:,j) = sum(Dyn%Grid(1,1)%vortend(:,j,future))/count(Dyn%Grid(1,1)%vortend(:,j,future) > -10000)
 pvor_advecmean(:,j) = sum(Dyn%Grid(1,1)%pvor_advec(:,j,future))/count(Dyn%Grid(1,1)%pvor_advec(:,j,future) > -10000)
 rvor_advecmean(:,j) = sum(Dyn%Grid(1,1)%rvor_advec(:,j,future))/count(Dyn%Grid(1,1)%rvor_advec(:,j,future) > -10000)
 vortruncmean(:,j) = sum(Dyn%Grid(1,1)%vortrunc(:,j,future))/count(Dyn%Grid(1,1)%vortrunc(:,j,future) > -10000) 
 vorlindamp_emean(:,j) = sum(Dyn%Grid(1,1)%vorlindamp_e(:,j,future))/count(Dyn%Grid(1,1)%vorlindamp_e(:,j,future) > -10000)
 vorlindamp_mmean(:,j) = sum(Dyn%Grid(1,1)%vorlindamp_m(:,j,future))/count(Dyn%Grid(1,1)%vorlindamp_m(:,j,future) > -10000)
 vorspecdampmean(:,j) = sum(Dyn%Grid(1,1)%vorspecdamp(:,j,future))/count(Dyn%Grid(1,1)%vorspecdamp(:,j,future) > -10000)
 vorstirmean(:,j) = sum(Dyn%Grid(1,1)%vorstir(:,j,future))/count(Dyn%Grid(1,1)%vorstir(:,j,future) > -10000)
 vorrobdampmean(:,j) = sum(Dyn%Grid(1,1)%vorrobdamp(:,j,future))/count(Dyn%Grid(1,1)%vorrobdamp(:,j,future) > -10000)
end do
 vorprime = Dyn%Grid(1,1)%vor(:,:,future)-vormean(:,:)
 vortendprime = Dyn%Grid(1,1)%vortend(:,:,future)-vortendmean(:,:)
 pvor_advecprime = Dyn%Grid(1,1)%pvor_advec(:,:,future)-pvor_advecmean(:,:)
 rvor_advecprime = Dyn%Grid(1,1)%rvor_advec(:,:,future)-rvor_advecmean(:,:)
 vortruncprime = Dyn%Grid(1,1)%vortrunc(:,:,future)-vortruncmean(:,:)
 vorlindamp_eprime = Dyn%Grid(1,1)%vorlindamp_e(:,:,future)-vorlindamp_emean(:,:)
 vorlindamp_mprime = Dyn%Grid(1,1)%vorlindamp_m(:,:,future)-vorlindamp_mmean(:,:)
 vorspecdampprime = Dyn%Grid(1,1)%vorspecdamp(:,:,future)-vorspecdampmean(:,:)
 vorstirprime = Dyn%Grid(1,1)%vorstir(:,:,future)-vorstirmean(:,:)
 vorrobdampprime = Dyn%Grid(1,1)%vorrobdamp(:,:,future)-vorrobdampmean(:,:)

! Spectral Terms
vorprimek(1:je/2+1,:) = fft_grid_to_fourier(vorprime)
vorprimek(je/2+2:je+1,:) = conjg(vorprimek(je/2+1:2:-1,:))
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ------------------------------ ENERGY MEAN_ENERGY & ENSTROPHY ----------------------------
Dyn%Grid(1,1)%energy(:,:,future) = (Dyn%Grid(1,1)%u(:,:,future)**2+Dyn%Grid(1,1)%v(:,:,future)**2)
!print *, "ENERGY", Dyn%Grid(1,1)%energy(10,10,future)
Dyn%Grid(1,1)%energy_mean(:,:,future) = (umean**2+vmean**2)
Dyn%Grid(1,1)%enstrophy(:,:,future) = vorprime*vorprime
!print *, "ENSTROPHY", Dyn%Grid(1,1)%enstrophy(22,47,future)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -------------------------- U,V,VOR, ENERGY & ENSTROPHY TENDENCIES ------------------------------
if(previous == current) then
Dyn%Grid(1,1)%utend(:,:,future) = (Dyn%Grid(1,1)%u(:,:,future)-Dyn%Grid(1,1)%u(:,:,current))/delta_t
Dyn%Grid(1,1)%vtend(:,:,future) = (Dyn%Grid(1,1)%v(:,:,future)-Dyn%Grid(1,1)%v(:,:,current))/delta_t
Dyn%Grid(1,1)%vortend(:,:,future) = (Dyn%Grid(1,1)%vor(:,:,future)-Dyn%Grid(1,1)%vor(:,:,current))/delta_t
Dyn%Grid(1,1)%vormeantend(:,:,future) = (vormean-vormean_curr)/delta_t
Dyn%Grid(1,1)%vorprimetend(:,:,future) = (vorprime-vorprime_curr)/delta_t
energytendcheck(:,:) = (Dyn%Grid(1,1)%energy(:,:,future)-Dyn%Grid(1,1)%energy(:,:,current))/delta_t
energytendcheck_mean(:,:) = (Dyn%Grid(1,1)%energy_mean(:,:,future)-Dyn%Grid(1,1)%energy_mean(:,:,current))/delta_t
enstrophytendcheck(:,:) = (Dyn%Grid(1,1)%enstrophy(:,:,future)-Dyn%Grid(1,1)%enstrophy(:,:,current))/delta_t
Dyn%Grid(1,1)%vormean_y_tend(:,:,future) = (Dyn%Grid(1,1)%vormean_y(:,:,future)-Dyn%Grid(1,1)%vormean_y(:,:,current))/delta_t
!vorprimetendk(:,:) = (vorprimek-vorprimek_curr)/delta_t
else
Dyn%Grid(1,1)%utend(:,:,future) = (Dyn%Grid(1,1)%u(:,:,future)-uprev)/delta_t
Dyn%Grid(1,1)%vtend(:,:,future) = (Dyn%Grid(1,1)%v(:,:,future)-vprev)/delta_t
Dyn%Grid(1,1)%vortend(:,:,future) = (Dyn%Grid(1,1)%vor(:,:,future)-vorprev)/delta_t
Dyn%Grid(1,1)%vormeantend(:,:,future) = (vormean-vormean_prev)/delta_t
Dyn%Grid(1,1)%vorprimetend(:,:,future) = (vorprime-vorprime_prev)/delta_t
energytendcheck(:,:) = (Dyn%Grid(1,1)%energy(:,:,future)-energyprev)/delta_t ! These two variables act as checks (confirm that the tendency calculated using the product rule is the same as central differences)
energytendcheck_mean(:,:) = (Dyn%Grid(1,1)%energy_mean(:,:,future)-energymean_prev)/delta_t ! These won't be correct for first two time steps (didn't put energyprev in restart file). Thats ok, just a check.
enstrophytendcheck(:,:) = (Dyn%Grid(1,1)%enstrophy(:,:,future)-enstrophyprev)/delta_t ! See above comment
Dyn%Grid(1,1)%vormean_y_tend(:,:,future) = (Dyn%Grid(1,1)%vormean_y(:,:,future)-vormean_y_prev)/delta_t
!vorprimetendk(:,:) = (vorprimek-vorprimek_prev)/delta_t
endif
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!print *, "HELLO", Dyn%Grid(1,1)%vortend(22,47,future) - (Dyn%Grid(1,1)%vormeantend(22,47,future) + Dyn%Grid(1,1)%vorprimetend(22,47,future))
!print *, "VOTMEAN TEND", Dyn%Grid(1,1)%vormeantend(22,47,future)
!print *, "VORPRIME TEND", Dyn%Grid(1,1)%vorprimetend(22,47,future)
! -------------------------------- ENERGY BUDGET TERMS --------------------------------------------
!============= COMMENT ON ENERGY CALCULATIONS ======================
! The U^2 and V^2 budget are calculated using the product rule (this way we can leverage all the work/code that went into creating the U and V budgets). 
! The central difference anlogy to the product rule is 
!	D(a*b)_n = D(a)_n*b_{n+1} + a_n*D(b)_n
! where D(a)_n = a_{n+1} - a_{n-1} is the central difference operator.
! For a = b this simplifies to 
!	D(a^2)_n = D(a)_n*(a_{n+1} + a_{n-1})
! Calculating u2tend in these two different ways confirms they are exact.
! Because the above product rule demands knowledge of u @ t=n+1 I calculate all the terms at the end of the subroutine, when this value is known.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Dyn%Grid(1,1)%energy_tend(:,:,future) = Dyn%Grid(1,1)%utend(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) + Dyn%Grid(1,1)%vtend(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_TEND", Dyn%Grid(1,1)%energy_tend(10,10,future)
!print *, "ENERGYTENDCHECK", energytendcheck(10,10)
Dyn%Grid(1,1)%energy_voradvec(:,:,future) = Dyn%Grid(1,1)%vq(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) - Dyn%Grid(1,1)%uq(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_VORADVEC", Dyn%Grid(1,1)%energy_voradvec(10,10,future)
Dyn%Grid(1,1)%energy_gradterm(:,:,future) = Dyn%Grid(1,1)%dx_u2_v2(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) + Dyn%Grid(1,1)%dy_u2_v2(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_GRADTERM", Dyn%Grid(1,1)%energy_gradterm(10,10,future)
Dyn%Grid(1,1)%energy_trunc(:,:,future) = Dyn%Grid(1,1)%utrunc(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) + Dyn%Grid(1,1)%vtrunc(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_TRUNC", Dyn%Grid(1,1)%energy_trunc(10,10,future)
Dyn%Grid(1,1)%energy_lindamp_m(:,:,future) = Dyn%Grid(1,1)%ulindamp_m(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future))
!print *, "ENERGY_LINDAMP_M", Dyn%Grid(1,1)%energy_lindamp_m(10,10,future)
Dyn%Grid(1,1)%energy_lindamp_e(:,:,future) = Dyn%Grid(1,1)%ulindamp_e(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) + Dyn%Grid(1,1)%vlindamp(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_LINDAMP_E", Dyn%Grid(1,1)%energy_lindamp_e(10,10,future)
Dyn%Grid(1,1)%energy_specdamp(:,:,future) = Dyn%Grid(1,1)%uspecdamp(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) + Dyn%Grid(1,1)%vspecdamp(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_SPECDAMP", Dyn%Grid(1,1)%energy_specdamp(10,10,future)
Dyn%Grid(1,1)%energy_stir(:,:,future) = Dyn%Grid(1,1)%ustir(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) + Dyn%Grid(1,1)%vstir(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_STIR", Dyn%Grid(1,1)%energy_stir(10,10,future)
Dyn%Grid(1,1)%energy_robdamp(:,:,future) = Dyn%Grid(1,1)%urobdamp(:,:,future)*(uprev(:,:)+Dyn%Grid(1,1)%u(:,:,future)) + Dyn%Grid(1,1)%vrobdamp(:,:,future)*(vprev(:,:)+Dyn%Grid(1,1)%v(:,:,future))
!print *, "ENERGY_ROBDAMP", Dyn%Grid(1,1)%energy_robdamp(10,10,future)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ----------------------------- MEAN ENERGY BUDGET TERMS ------------------------------------------
! ======== Compute energy_tend_mean in a different way =============
!energymean_prev = (umean_prev**2+vmean_prev**2)
!Dyn%Grid(1,1)%energy_mean(:,:,future) = (umean**2+vmean**2)
!print *, "ENERGY_MEAN", Dyn%Grid(1,1)%energy_mean(10,10,future)
!if(previous == current) then
!energytendcheck_mean(:,:) = (Dyn%Grid(1,1)%energy_mean(:,:,future)-Dyn%Grid(1,1)%energy_mean(:,:,current))/delta_t
!else
!energytendcheck_mean(:,:) = (Dyn%Grid(1,1)%energy_mean(:,:,future)-energymean_prev)/delta_t
!endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute \bar{E} = \bar{u}^2 + \bar{v}^2 budget terms using product rule
Dyn%Grid(1,1)%energy_tend_mean(:,:,future) = utendmean(:,:)*(umean_prev(:,:)+umean(:,:)) + vtendmean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYTEND_MEAN", Dyn%Grid(1,1)%energy_tend_mean(1,10,future)
!print *, "ENERGYTENDCHECK_MEAN", energytendcheck_mean(1,10)
Dyn%Grid(1,1)%energy_voradvec_mean(:,:,future) = vqmean(:,:)*(umean_prev(:,:)+umean(:,:)) - uqmean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYVORADVEC_MEAN", Dyn%Grid(1,1)%energy_voradvec_mean(1,10,future)
Dyn%Grid(1,1)%energy_gradterm_mean(:,:,future) = dx_u2_v2mean(:,:)*(umean_prev(:,:)+umean(:,:)) + dy_u2_v2mean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYGRADTERM_MEAN", Dyn%Grid(1,1)%energy_gradterm_mean(1,10,future)
Dyn%Grid(1,1)%energy_trunc_mean(:,:,future) = utruncmean(:,:)*(umean_prev(:,:)+umean(:,:)) + vtruncmean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYTRUNC_MEAN", Dyn%Grid(1,1)%energy_trunc_mean(1,10,future)
Dyn%Grid(1,1)%energy_lindamp_m_mean(:,:,future) = ulindamp_mmean(:,:)*(umean_prev(:,:)+umean(:,:))
!print *, "ENERGYLINDAMP_M_MEAN", Dyn%Grid(1,1)%energy_lindamp_m_mean(1,10,future)
Dyn%Grid(1,1)%energy_lindamp_e_mean(:,:,future) = ulindamp_emean(:,:)*(umean_prev(:,:)+umean(:,:)) + vlindampmean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYLINDAMP_E_MEAN", Dyn%Grid(1,1)%energy_lindamp_e_mean(1,10,future)
Dyn%Grid(1,1)%energy_specdamp_mean(:,:,future) = uspecdampmean(:,:)*(umean_prev(:,:)+umean(:,:)) + vspecdampmean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYSPECDAMP_MEAN", Dyn%Grid(1,1)%energy_specdamp_mean(1,10,future)
Dyn%Grid(1,1)%energy_stir_mean(:,:,future) = ustirmean(:,:)*(umean_prev(:,:)+umean(:,:)) + vstirmean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYSTIR_MEAN", Dyn%Grid(1,1)%energy_stir_mean(1,10,future)
Dyn%Grid(1,1)%energy_robdamp_mean(:,:,future) = urobdampmean(:,:)*(umean_prev(:,:)+umean(:,:)) + vrobdampmean(:,:)*(vmean_prev(:,:)+vmean(:,:))
!print *, "ENERGYROBDAMP_MEAN", Dyn%Grid(1,1)%energy_robdamp_mean(1,10,future)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -------------------------------- ENSTROPHY BUDGET  -----------------------------------------------
!============= COMMENT ON ENSTROPHY CALCULATIONS ======================
! The VOR'^2 budget are calculated using the product rule (this way we can leverage all the work/code that went into creating the VOR budget). 
! The central difference anlogy to the product rule is 
!	D(a*b)_n = D(a)_n*b_{n+1} + a_n*Dn
! where D(a)_n = a_{n+1} - a_{n-1} is the central difference operator.
! For a = b this simplifies to 
!	D(a^2)_n = D(a)_n*(a_{n+1} + a_{n-1})
! Calculating u2tend in these two different ways confirms they are exact.
! Because the above product rule demands knowledge of u @ t=n+1 I calculate all the terms at the end of the subroutine, when this value is known.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Dyn%Grid(1,1)%enstrophy_tend(:,:,future) = 0.5*Dyn%Grid(1,1)%vorprimetend(:,:,future)*(vorprime_prev+vorprime)
!print *, "ENSTROPHY_TEND", Dyn%Grid(1,1)%enstrophy_tend(22,47,future)
!print *, "ENSTROPHYTENDCHECK", enstrophytendcheck(22,47)
!print *, "ENSTROPHYTENDCHECK- ENSTROPHY_TEND", enstrophytendcheck(22,47) - Dyn%Grid(1,1)%enstrophy_tend(22,47,future)
Dyn%Grid(1,1)%enstrophy_dt_vor(:,:,future) = Dyn%Grid(1,1)%dt_vor(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_DTVOR", Dyn%Grid(1,1)%enstrophy_dt_vor(22,47,future)
Dyn%Grid(1,1)%enstrophy_pvor_advec(:,:,future) = Dyn%Grid(1,1)%pvor_advec(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_PVORADVEC", Dyn%Grid(1,1)%enstrophy_pvor_advec(22,47,future)
Dyn%Grid(1,1)%enstrophy_rvormean_advec(:,:,future) = Dyn%Grid(1,1)%rvormean_advec(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_RVORMEAN_ADVEC", Dyn%Grid(1,1)%enstrophy_rvormean_advec(22,47,future)
Dyn%Grid(1,1)%enstrophy_rvorprime_advec(:,:,future) = Dyn%Grid(1,1)%rvorprime_advec(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_RVORPRIME_ADVEC", Dyn%Grid(1,1)%enstrophy_rvorprime_advec(22,47,future)
Dyn%Grid(1,1)%enstrophy_trunc(:,:,future) = Dyn%Grid(1,1)%vortrunc(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_TRUNC", Dyn%Grid(1,1)%enstrophy_trunc(22,47,future)
Dyn%Grid(1,1)%enstrophy_lindamp_m(:,:,future) = Dyn%Grid(1,1)%vorlindamp_m(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_LINDAMP_M", Dyn%Grid(1,1)%enstrophy_lindamp_m(22,47,future)
Dyn%Grid(1,1)%enstrophy_lindamp_e(:,:,future) = Dyn%Grid(1,1)%vorlindamp_e(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_LINDAMP_E", Dyn%Grid(1,1)%enstrophy_lindamp_e(22,47,future)
Dyn%Grid(1,1)%enstrophy_specdamp(:,:,future) = Dyn%Grid(1,1)%vorspecdamp(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_SPECDAMP", Dyn%Grid(1,1)%enstrophy_specdamp(22,47,future)
Dyn%Grid(1,1)%enstrophy_stir(:,:,future) = Dyn%Grid(1,1)%vorstir(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_STIR", Dyn%Grid(1,1)%enstrophy_stir(22,47,future)
Dyn%Grid(1,1)%enstrophy_robdamp(:,:,future) = Dyn%Grid(1,1)%vorrobdamp(:,:,future)*vorprime_curr
!print *, "ENSTROPHY_ROBDAMP", Dyn%Grid(1,1)%enstrophy_robdamp(22,47,future)

Dyn%Grid(1,1)%enstrophy_rvormean_adevc_Matlab(:,:,future) = (Dyn%Grid(1,1)%u(:,:,current)*Dyn%Grid(1,1)%vormean_x(:,:,future) + Dyn%Grid(1,1)%v(:,:,current)*Dyn%Grid(1,1)%vormean_y(:,:,future))*vorprime_curr
Dyn%Grid(1,1)%enstrophy_rvorprime_adevc_Matlab(:,:,future) = (Dyn%Grid(1,1)%u(:,:,current)*Dyn%Grid(1,1)%vorprime_x(:,:,future) + Dyn%Grid(1,1)%v(:,:,current)*Dyn%Grid(1,1)%vorprime_y(:,:,future))*vorprime_curr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ----------------------------- ENERGY MEAN_ENERGY & ENSTROPHY FUTURE VALUES -----------------------
! Need these for the product rule (see documentation) and for calculating energy_mean and enstrophy
do j = js, je 
 enstrophytrunc_mean(:,j)= sum(Dyn%Grid(1,1)%enstrophy_trunc(:,j,future))/count(Dyn%Grid(1,1)%u(:,j,future) > -10000)
 enstrophystir_mean(:,j)= sum(Dyn%Grid(1,1)%enstrophy_stir(:,j,future))/count(Dyn%Grid(1,1)%u(:,j,future) > -10000)
end do
enstrophytrunc_prime = Dyn%Grid(1,1)%enstrophy_trunc(:,:,future)-enstrophytrunc_mean(:,:)
enstrophystir_prime = Dyn%Grid(1,1)%enstrophy_stir(:,:,future)-enstrophystir_mean(:,:)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ------------------------------ ENSTROPHY SPECTRAL BUDGET ----------------------------

!============= COMMENT ON FORTRAN90 FFT =====================
! fft routine can be found in /src/shared/fft.f90
! fft_grid_to_fourier takes an fft along the first dimension of a 2D or 3D array
! The complex fourier components are passed in the following format.
!
!        fourier (1)     = cmplx ( a(0), b(0) )
!        fourier (2)     = cmplx ( a(1), b(1) )
!            :              :
!            :              :
!        fourier (n/2+1) = cmplx ( a(n/2), b(n/2) )
!
!   where n = length of each real transform
!
! There are two differences between the FFT subroutine in Fortran and the equivalent in Matlab:
!
! 1) The fortran fft has a factor of 1/n outside the forward fft while with the matlab fft this is associated with a backward (inverse) fft.
! Multiply the fortran output by n to get agreement. 
!
! 2) The Fortran transformed vector is length N+1 while Matlabs transformed vector is N. The extra element in Fortran is because it keeps the positive N/2th harmonic, Matlab does not ie.
!
! 		Matlab:  [0 1 ... N/2-1     -N/2 -N/2-1 ... -1] = N elements
!		Fortran: [0 1 ... N/2-1 N/2 -N/2 -N/2-1 ... -1] = N+1 elements
!
! There are a couple of ways to handle this; we could just modify the kmat vector of harmonics, or remove the extra elements from the Fortran output, removing the extra 64th harmonic shouldn't make a big
! difference (tiny higher frequency).
!
! Because barotropic_diagnostics only handles real Dyn%grid variable I save the resultant fourier transform as:
! 	[real(fft_grid_to_fourier) imag(fft_grid_to_fourier)]
! The first n/2+1 are the real components, the last n/2-1 components are the 2:n/2-1 complex components. 
! The 1st and n/2th complex components are always zero for the fft of a real variable. These zeros need to be added in when the .nc file is loaded into matlab (see .m code) 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

kmat = [0:ie/2, -ie/2:-1]
!lmat = [0:je/2-1, -je/2:-1]

! BUILD UP SPECTRAL ARRAYS; Has (almost) the same structure as Matlab code ie wavenumber takes the form [0 1 2 ... N/2 -N/2 -N/2+1 ... -1]
uk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%ucurr(:,:,future))
uk(ie/2+1:ie,:) = conjg(uk(ie/2+1:2:-1,:))
vk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vcurr(:,:,future))
vk(ie/2+1:ie,:) = conjg(vk(ie/2+1:2:-1,:))
v_yk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%v_y(:,:,future))
v_yk(ie/2+1:ie,:) = conjg(v_yk(ie/2+1:2:-1,:))
!vbetak(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%v(:,:,future)*Dyn%Grid(1,1)%beta(:,:,future));
!vbetak(ie/2+1:ie,:) = conjg(vbetak(ie/2+1:2:-1,:))
pvor_adveck(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%pvor_advec(:,:,future));
pvor_adveck(ie/2+1:ie,:) = conjg(pvor_adveck(ie/2+1:2:-1,:))
rvormean_adveck(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%rvormean_advec(:,:,future));
rvormean_adveck(ie/2+1:ie,:) = conjg(rvormean_adveck(ie/2+1:2:-1,:))
vorprime_xk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorprime_x(:,:,future))
vorprime_xk(ie/2+1:ie,:) = conjg(vorprime_xk(ie/2+1:2:-1,:))
vorprime_yk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorprime_y(:,:,future))
vorprime_yk(ie/2+1:ie,:) = conjg(vorprime_yk(ie/2+1:2:-1,:))
vorprimetendk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorprimetend(:,:,future))
vorprimetendk(ie/2+1:ie,:) = conjg(vorprimetendk(ie/2+1:2:-1,:))
vor_xk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vor_x(:,:,future))
vor_xk(ie/2+1:ie,:) = conjg(vor_xk(ie/2+1:2:-1,:))
vor_yk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vor_y(:,:,future))
vor_yk(ie/2+1:ie,:) = conjg(vor_yk(ie/2+1:2:-1,:))
vdy_vormeank(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%v(:,:,future)*Dyn%Grid(1,1)%vormean_y(:,:,future));
vdy_vormeank(ie/2+1:ie,:) = conjg(vdy_vormeank(ie/2+1:2:-1,:))
vortrunck(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vortrunc(:,:,future))
vortrunck(ie/2+1:ie,:) = conjg(vortrunck(ie/2+1:2:-1,:))
!print *, "VORTRUNCK", vortrunck(6,47)
vorlindamp_ek(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorlindamp_e(:,:,future))
vorlindamp_ek(ie/2+1:ie,:) = conjg(vorlindamp_ek(ie/2+1:2:-1,:))
vorlindamp_mk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorlindamp_m(:,:,future))
vorlindamp_mk(ie/2+1:ie,:) = conjg(vorlindamp_mk(ie/2+1:2:-1,:))
vorrobdampk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorrobdamp(:,:,future))
vorrobdampk(ie/2+1:ie,:) = conjg(vorrobdampk(ie/2+1:2:-1,:))
vorspecdampk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorspecdamp(:,:,future))
vorspecdampk(ie/2+1:ie,:) = conjg(vorspecdampk(ie/2+1:2:-1,:))
vorstirk(1:ie/2+1,:) = fft_grid_to_fourier(Dyn%Grid(1,1)%vorstir(:,:,future))
vorstirk(ie/2+1:ie,:) = conjg(vorstirk(ie/2+1:2:-1,:))

!do p=1,size(kmat)
!    vor_xkII(p,:) = ((0,1)*kmat(p)/(radius*cos(deg_lat(:)*pi/180)))*vork(p,:)
!end do
!
!dlat = (cshift(deg_lat,1) - cshift(deg_lat,-1))*pi/180
!
!print *, dlat
!
!do j = is, ie 
!  dvork(j,:) = cshift(vork(j,:),1) - cshift(vork(j,:),-1)
!  vor_ykII(j,:) = (radius**-1)*dvork(j,:))/dlat
!end do

Dyn%Grid(1,1)%vork(1:ie/2+1,:,future) = real(vork(1:ie/2+1,:)) ! First n/2+1 positions belong to real values
Dyn%Grid(1,1)%vork(ie/2+2:ie,:,future)= aimag(vork(2:ie/2,:)) ! last n/2-1 positions belong to imaginary values.

Dyn%Grid(1,1)%vork_curr(1:ie/2+1,:,future) = real(vork_curr(1:ie/2+1,:)) ! First n/2+1 positions belong to real values
Dyn%Grid(1,1)%vork_curr(ie/2+2:ie,:,future)= aimag(vork_curr(2:ie/2,:)) ! last n/2-1 positions belong to imaginary values.

Dyn%Grid(1,1)%vork_prev(1:ie/2+1,:,future) = real(vork_prev(1:ie/2+1,:)) ! First n/2+1 positions belong to real values
Dyn%Grid(1,1)%vork_prev(ie/2+2:ie,:,future)= aimag(vork_prev(2:ie/2,:)) ! last n/2-1 positions belong to imaginary values.

! --------------------------------------------------
vorprime_fluxk(1:ie/2+1,:) = ie**2*vk*conjg(vorprimek_curr)
enstrophy_tendk = ie**2*0.5*vorprimetendk*conjg(vorprimek_prev + vorprimek) ! Factors of ie**2 in so that the Fortan fft subroutine is equivalent to the Matlab fft 
enstrophy_pvor_adveck = ie**2*pvor_adveck*conjg(vorprimek_curr)
enstrophy_rvormean_adveck = ie**2*rvormean_adveck*conjg(vorprimek_curr)
enstrophy_trunck = ie**2*vortrunck*conjg(vorprimek_curr)
enstrophy_lindamp_ek = ie**2*vorlindamp_ek*conjg(vorprimek_curr)
enstrophy_lindamp_mk = ie**2*vorlindamp_mk*conjg(vorprimek_curr)
enstrophy_specdampk = ie**2*vorspecdampk*conjg(vorprimek_curr)
enstrophy_robdampk = ie**2*vorrobdampk*conjg(vorprimek_curr)
enstrophy_stirk = ie**2*vorstirk*conjg(vorprimek_curr)

Dyn%Grid(1,1)%vortendk(1:ie/2+1,:,future) = real(vorprimetendk(1:ie/2+1,:))
Dyn%Grid(1,1)%vortendk(ie/2+2:ie,:,future)= aimag(vorprimetendk(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_tendk(1:ie/2+1,:,future) = real(enstrophy_tendk(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_tendk(ie/2+2:ie,:,future)= aimag(enstrophy_tendk(2:ie/2,:))
Dyn%Grid(1,1)%vorprime_fluxk(1:ie/2+1,:,future) = real(vorprime_fluxk(1:ie/2+1,:))
Dyn%Grid(1,1)%vorprime_fluxk(ie/2+2:ie,:,future)= aimag(vorprime_fluxk(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_pvor_adveck(1:ie/2+1,:,future) = real(enstrophy_pvor_adveck(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_pvor_adveck(ie/2+2:ie,:,future)= aimag(enstrophy_pvor_adveck(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_rvormean_adveck(1:ie/2+1,:,future) = real(enstrophy_rvormean_adveck(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_rvormean_adveck(ie/2+2:ie,:,future)= aimag(enstrophy_rvormean_adveck(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_trunck(1:ie/2+1,:,future) = real(enstrophy_trunck(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_trunck(ie/2+2:ie,:,future)= aimag(enstrophy_trunck(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_lindamp_ek(1:ie/2+1,:,future) = real(enstrophy_lindamp_ek(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_lindamp_ek(ie/2+2:ie,:,future)= aimag(enstrophy_lindamp_ek(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_lindamp_mk(1:ie/2+1,:,future) = real(enstrophy_lindamp_mk(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_lindamp_mk(ie/2+2:ie,:,future)= aimag(enstrophy_lindamp_mk(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_specdampk(1:ie/2+1,:,future) = real(enstrophy_specdampk(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_specdampk(ie/2+2:ie,:,future)= aimag(enstrophy_specdampk(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_robdampk(1:ie/2+1,:,future) = real(enstrophy_robdampk(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_robdampk(ie/2+2:ie,:,future)= aimag(enstrophy_robdampk(2:ie/2,:))
Dyn%Grid(1,1)%enstrophy_stirk(1:ie/2+1,:,future) = real(enstrophy_stirk(1:ie/2+1,:))
Dyn%Grid(1,1)%enstrophy_stirk(ie/2+2:ie,:,future)= aimag(enstrophy_stirk(2:ie/2,:))

Dyn%Grid(1,1)%vorprime_yk(1:ie/2+1,:,future) = real(vorprime_yk(1:ie/2+1,:))
Dyn%Grid(1,1)%vorprime_yk(ie/2+2:ie,:,future)= aimag(vorprime_yk(2:ie/2,:))
Dyn%Grid(1,1)%v_yk(1:ie/2+1,:,future) = real(v_yk(1:ie/2+1,:))
Dyn%Grid(1,1)%v_yk(ie/2+2:ie,:,future)= aimag(v_yk(2:ie/2,:))
Dyn%Grid(1,1)%vorprimek_curr(1:ie/2+1,:,future) = real(vorprimek_curr(1:ie/2+1,:))
Dyn%Grid(1,1)%vorprimek_curr(ie/2+2:ie,:,future)= aimag(vorprimek_curr(2:ie/2,:))
Dyn%Grid(1,1)%vk(1:ie/2+1,:,future) = real(vk(1:ie/2+1,:))
Dyn%Grid(1,1)%vk(ie/2+2:ie,:,future)= aimag(vk(2:ie/2,:))

if(spec_calc > 0) then
! ------------ EXPERIMENT CALCULATING Y-DERIVATIVE USING FFT'S IN SUBROUTINE --------------
!
! call vor_div_from_uv_grid(Dyn%Grid(1,1)%vork(:,:,current)*Dyn%Grid(1,1)%vk(:,:,current), zerog, dt_vors_rub, dt_divs_rub) 
! call trans_spherical_to_grid(dt_vors_rub, Dyn%Grid(1,1)%rvor_adveck(:,:,future))
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!do p=1,size(kmat)
!    dx_vorprimek(p,:) = ((0,1)*kmat(p)/(radius*cos(deg_lat(:)*pi/180)))*vorprimek_curr(p,:)
!end do
!
!do p=1,size(lmat)
!    rubl(p,:) = ((0,1)*lmat(p)/(radius))*vorprimel_curr(p,:)
!end do
!
!dy_vorprimek =  2*fft_grid_to_fourier(transpose(fft_fourier_to_grid(rubl)))

!print *, "VORPRIME_YK", Dyn%Grid(1,1)%tester(1:5,:,future)
!print *, size(kmat)

!do p=1,size(kmat)
!       do q=1,size(kmat)
!             do k=1,size(kmat)
!                print *, [p, q, k]
!                  if (kmat(p)+kmat(q)==kmat(k)) then
!                       enstrophy_rvorprime_adveck(p,q,k,:) = (uk(p,:)*vorprime_xk(q,:) + vk(p,:)*vorprime_yk(q,:))*conjg(vorprimek_curr(k,:))
!                      tester(p,q,k,:) = vk(p,:)*vorprime_yk(q,:)*conjg(vorprimek_curr(k,:)) + 0.5*v_yk(p,:)*vorprimek_curr(q,:)*conjg(vorprimek_curr(k,:))
!		       tester2(p,q,k,:) =  vk(p,:)*vorprime_yk(q,:)*conjg(vorprimek_curr(k,:))
!                      tester3(p,q,k,:) =  0.5*v_yk(p,:)*vorprimek_curr(q,:)*conjg(vorprimek_curr(k,:))
!            	  else
!                      enstrophy_rvorprime_adveck(p,q,k,:) = (0,0)
!                      tester(p,q,k,:) = (0,0)
!                      tester2(p,q,k,:) = (0,0)
!                      tester3(p,q,k,:) = (0,0)
!                  end if
!            end do
!      end do
!end do

! ---------------- SPECTRAL SPACE vor_divs_from_uv_grid ------------------------------------

do p=1,size(kmat)
       do q=1,size(kmat)
	   vor_flux_spec(p,q,:) = vk(p,:)*vorprimek_curr(q,:)      
	end do
end do

do p=1,size(kmat)
	 call vor_div_from_uv_grid (real(vor_flux_spec(p,1:ie,:)), zerog, dt_vors_rub, dt_divs_rub) 
         call trans_spherical_to_grid(dt_vors_rub, spec_deriv_vq_real(p,:,:))
	 call vor_div_from_uv_grid (aimag(vor_flux_spec(p,1:ie,:)), zerog, dt_vors_rub, dt_divs_rub) 
         call trans_spherical_to_grid(dt_vors_rub, spec_deriv_vq_imag(p,:,:))
end do

do p=1,size(kmat)
    do q=1,size(kmat)
	  do k=1,size(kmat)
		if (kmat(p)+kmat(q)==kmat(k)) then
	  	   enstrophy_rvorprime_adveck2(p,q,k,:) = cmplx(spec_deriv_vq_real(p,q,:),spec_deriv_vq_imag(p,q,:))*conjg(vorprimek_curr(k,:))
		else
		   enstrophy_rvorprime_adveck2(p,q,k,:) = (0,0)
		end if
	  end do
    end do
end do

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!print *, "enstrophy_rvorprime_adveck", tester2(1:5,1,1,48)

!I_pp = sum(sum(enstrophy_rvorprime_adveck(2:4,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie-2:ie,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(2:4,ie-2:ie,:,:),1),1) +  & 
!	sum(sum(enstrophy_rvorprime_adveck(ie-2:ie,ie-2:ie,:,:),1),1)

!I_ps = sum(sum(enstrophy_rvorprime_adveck(2:4,5:13,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(2:4,ie-12:ie-3,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie-2:ie,ie-12:ie-3,:,:),1),1) + &
!	sum(sum(enstrophy_rvorprime_adveck(ie-2:ie,5:13,:,:),1),1)

!I_ph = sum(sum(enstrophy_rvorprime_adveck(2:4,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie-2:ie,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(2:4,ie/2+1:ie-13,:,:),1),1) + &
!	sum(sum(enstrophy_rvorprime_adveck(ie-2:ie,ie/2+1:ie-13,:,:),1),1)

!I_sp = sum(sum(enstrophy_rvorprime_adveck(5:13,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie-12:ie-3,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie-12:ie-3,ie-2:ie,:,:),1),1) + &
!	sum(sum(enstrophy_rvorprime_adveck(5:13,ie-2:ie,:,:),1),1)

!I_ss = sum(sum(enstrophy_rvorprime_adveck(5:13,5:13,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(5:13,ie-12:ie-3,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie-12:ie-3,5:13,:,:),1),1) + &!
!	sum(sum(enstrophy_rvorprime_adveck(ie-12:ie-3,ie-12:ie-3,:,:),1),1)

!I_sh = sum(sum(enstrophy_rvorprime_adveck(5:13,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(5:13,ie/2+1:ie-13,:,:),1),1) + &
!	 sum(sum(enstrophy_rvorprime_adveck(ie-12:ie-3,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie-12:ie-3,ie/2+1:ie-13,:,:),1),1)

!I_hp = sum(sum(enstrophy_rvorprime_adveck(14:ie/2,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(14:ie/2,ie-2:ie,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie/2+1:ie-13,2:4,:,:),1),1) + &
!	sum(sum(enstrophy_rvorprime_adveck(ie/2+1:ie-13,ie-2:ie,:,:),1),1)

!I_hs = sum(sum(enstrophy_rvorprime_adveck(14:ie/2,5:13,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie/2+1:ie-13,5:13,:,:),1),1) + &
!	 sum(sum(enstrophy_rvorprime_adveck(14:ie/2,ie-12:ie-3,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie/2+1:ie-13,ie-12:ie-3,:,:),1),1)

!I_hh = sum(sum(enstrophy_rvorprime_adveck(14:ie/2,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(14:ie/2,ie/2+1:ie-13,:,:),1),1) + &
!	 sum(sum(enstrophy_rvorprime_adveck(ie/2+1:ie-13,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck(ie/2+1:ie-13,ie/2+1:ie-13,:,:),1),1)

I_pp = sum(sum(enstrophy_rvorprime_adveck2(2:4,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie-2:ie,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(2:4,ie-2:ie,:,:),1),1) +  & 
	sum(sum(enstrophy_rvorprime_adveck2(ie-2:ie,ie-2:ie,:,:),1),1)

I_ps = sum(sum(enstrophy_rvorprime_adveck2(2:4,5:13,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(2:4,ie-12:ie-3,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie-2:ie,ie-12:ie-3,:,:),1),1) + &
	sum(sum(enstrophy_rvorprime_adveck2(ie-2:ie,5:13,:,:),1),1)

I_ph = sum(sum(enstrophy_rvorprime_adveck2(2:4,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie-2:ie,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(2:4,ie/2+1:ie-13,:,:),1),1) + &
	sum(sum(enstrophy_rvorprime_adveck2(ie-2:ie,ie/2+1:ie-13,:,:),1),1)

I_sp = sum(sum(enstrophy_rvorprime_adveck2(5:13,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie-12:ie-3,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie-12:ie-3,ie-2:ie,:,:),1),1) + &
	sum(sum(enstrophy_rvorprime_adveck2(5:13,ie-2:ie,:,:),1),1)

I_ss = sum(sum(enstrophy_rvorprime_adveck2(5:13,5:13,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(5:13,ie-12:ie-3,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie-12:ie-3,5:13,:,:),1),1) + &!
	sum(sum(enstrophy_rvorprime_adveck2(ie-12:ie-3,ie-12:ie-3,:,:),1),1)

I_sh = sum(sum(enstrophy_rvorprime_adveck2(5:13,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(5:13,ie/2+1:ie-13,:,:),1),1) + &
	 sum(sum(enstrophy_rvorprime_adveck2(ie-12:ie-3,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie-12:ie-3,ie/2+1:ie-13,:,:),1),1)

I_hp = sum(sum(enstrophy_rvorprime_adveck2(14:ie/2,2:4,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(14:ie/2,ie-2:ie,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie/2+1:ie-13,2:4,:,:),1),1) + &
	sum(sum(enstrophy_rvorprime_adveck2(ie/2+1:ie-13,ie-2:ie,:,:),1),1)

I_hs = sum(sum(enstrophy_rvorprime_adveck2(14:ie/2,5:13,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie/2+1:ie-13,5:13,:,:),1),1) + &
	 sum(sum(enstrophy_rvorprime_adveck2(14:ie/2,ie-12:ie-3,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie/2+1:ie-13,ie-12:ie-3,:,:),1),1)

I_hh = sum(sum(enstrophy_rvorprime_adveck2(14:ie/2,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(14:ie/2,ie/2+1:ie-13,:,:),1),1) + &
	 sum(sum(enstrophy_rvorprime_adveck2(ie/2+1:ie-13,14:ie/2,:,:),1),1) + sum(sum(enstrophy_rvorprime_adveck2(ie/2+1:ie-13,ie/2+1:ie-13,:,:),1),1)

I_0p = sum(enstrophy_rvorprime_adveck(1,2:4,:,:),1) + sum(enstrophy_rvorprime_adveck(1,ie-3:ie,:,:),1)
I_p0 = sum(enstrophy_rvorprime_adveck(2:4,1,:,:),1) + sum(enstrophy_rvorprime_adveck(ie-3:ie,1,:,:),1)
I_0s = sum(enstrophy_rvorprime_adveck(1,5:13,:,:),1) + sum(enstrophy_rvorprime_adveck(1,ie-12:ie-4,:,:),1)
I_s0 = sum(enstrophy_rvorprime_adveck(5:13,1,:,:),1) + sum(enstrophy_rvorprime_adveck(ie-12:ie-4,1,:,:),1)
I_0h = sum(enstrophy_rvorprime_adveck(1,14:ie/2,:,:),1) + sum(enstrophy_rvorprime_adveck(1,ie/2+1:ie-13,:,:),1)
I_h0 = sum(enstrophy_rvorprime_adveck(14:ie/2,1,:,:),1) + sum(enstrophy_rvorprime_adveck(ie/2+1:ie-13,1,:,:),1)
I_00 = enstrophy_rvorprime_adveck(1,1,:,:)

!I
Dyn%Grid(1,1)%I_pp(1:ie/2+1,:,future) = ie**2*real(I_pp(1:ie/2+1,:))
Dyn%Grid(1,1)%I_pp(ie/2+2:ie,:,future)= ie**2*aimag(I_pp(2:ie/2,:))
Dyn%Grid(1,1)%I_sp(1:ie/2+1,:,future) = ie**2*real(I_sp(1:ie/2+1,:))
Dyn%Grid(1,1)%I_sp(ie/2+2:ie,:,future)= ie**2*aimag(I_sp(2:ie/2,:))
Dyn%Grid(1,1)%I_ps(1:ie/2+1,:,future) = ie**2*real(I_ps(1:ie/2+1,:))
Dyn%Grid(1,1)%I_ps(ie/2+2:ie,:,future)= ie**2*aimag(I_ps(2:ie/2,:))
Dyn%Grid(1,1)%I_hp(1:ie/2+1,:,future) = ie**2*real(I_hp(1:ie/2+1,:))
Dyn%Grid(1,1)%I_hp(ie/2+2:ie,:,future)= ie**2*aimag(I_hp(2:ie/2,:))
Dyn%Grid(1,1)%I_ph(1:ie/2+1,:,future) = ie**2*real(I_ph(1:ie/2+1,:))
Dyn%Grid(1,1)%I_ph(ie/2+2:ie,:,future)= ie**2*aimag(I_ph(2:ie/2,:))
Dyn%Grid(1,1)%I_ss(1:ie/2+1,:,future) = ie**2*real(I_ss(1:ie/2+1,:))
Dyn%Grid(1,1)%I_ss(ie/2+2:ie,:,future)= ie**2*aimag(I_ss(2:ie/2,:))
Dyn%Grid(1,1)%I_hh(1:ie/2+1,:,future) = ie**2*real(I_hh(1:ie/2+1,:))
Dyn%Grid(1,1)%I_hh(ie/2+2:ie,:,future)= ie**2*aimag(I_hh(2:ie/2,:))
Dyn%Grid(1,1)%I_sh(1:ie/2+1,:,future) = ie**2*real(I_sh(1:ie/2+1,:))
Dyn%Grid(1,1)%I_sh(ie/2+2:ie,:,future)= ie**2*aimag(I_sh(2:ie/2,:))
Dyn%Grid(1,1)%I_hs(1:ie/2+1,:,future) = ie**2*real(I_hs(1:ie/2+1,:))
Dyn%Grid(1,1)%I_hs(ie/2+2:ie,:,future)= ie**2*aimag(I_hs(2:ie/2,:))
!I2
Dyn%Grid(1,1)%I2_pp(1:ie/2+1,:,future) = ie**2*real(I2_pp(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_pp(ie/2+2:ie,:,future)= ie**2*aimag(I2_pp(2:ie/2,:))
Dyn%Grid(1,1)%I2_sp(1:ie/2+1,:,future) = ie**2*real(I2_sp(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_sp(ie/2+2:ie,:,future)= ie**2*aimag(I2_sp(2:ie/2,:))
Dyn%Grid(1,1)%I2_ps(1:ie/2+1,:,future) = ie**2*real(I2_ps(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_ps(ie/2+2:ie,:,future)= ie**2*aimag(I2_ps(2:ie/2,:))
Dyn%Grid(1,1)%I2_hp(1:ie/2+1,:,future) = ie**2*real(I2_hp(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_hp(ie/2+2:ie,:,future)= ie**2*aimag(I2_hp(2:ie/2,:))
Dyn%Grid(1,1)%I2_ph(1:ie/2+1,:,future) = ie**2*real(I2_ph(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_ph(ie/2+2:ie,:,future)= ie**2*aimag(I2_ph(2:ie/2,:))
Dyn%Grid(1,1)%I2_ss(1:ie/2+1,:,future) = ie**2*real(I2_ss(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_ss(ie/2+2:ie,:,future)= ie**2*aimag(I2_ss(2:ie/2,:))
Dyn%Grid(1,1)%I2_hh(1:ie/2+1,:,future) = ie**2*real(I2_hh(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_hh(ie/2+2:ie,:,future)= ie**2*aimag(I2_hh(2:ie/2,:))
Dyn%Grid(1,1)%I2_sh(1:ie/2+1,:,future) = ie**2*real(I2_sh(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_sh(ie/2+2:ie,:,future)= ie**2*aimag(I2_sh(2:ie/2,:))
Dyn%Grid(1,1)%I2_hs(1:ie/2+1,:,future) = ie**2*real(I2_hs(1:ie/2+1,:))
Dyn%Grid(1,1)%I2_hs(ie/2+2:ie,:,future)= ie**2*aimag(I2_hs(2:ie/2,:))
! Other
Dyn%Grid(1,1)%I_0p(1:ie/2+1,:,future) = ie**2*real(I_0p(1:ie/2+1,:))
Dyn%Grid(1,1)%I_0p(ie/2+2:ie,:,future)= ie**2*aimag(I_0p(2:ie/2,:))
Dyn%Grid(1,1)%I_p0(1:ie/2+1,:,future) = ie**2*real(I_p0(1:ie/2+1,:))
Dyn%Grid(1,1)%I_p0(ie/2+2:ie,:,future)= ie**2*aimag(I_p0(2:ie/2,:))
Dyn%Grid(1,1)%I_0s(1:ie/2+1,:,future) = ie**2*real(I_0s(1:ie/2+1,:))
Dyn%Grid(1,1)%I_0s(ie/2+2:ie,:,future)= ie**2*aimag(I_0s(2:ie/2,:))
Dyn%Grid(1,1)%I_s0(1:ie/2+1,:,future) = ie**2*real(I_s0(1:ie/2+1,:))
Dyn%Grid(1,1)%I_s0(ie/2+2:ie,:,future)= ie**2*aimag(I_s0(2:ie/2,:))
Dyn%Grid(1,1)%I_0h(1:ie/2+1,:,future) = ie**2*real(I_0h(1:ie/2+1,:))
Dyn%Grid(1,1)%I_0h(ie/2+2:ie,:,future)= ie**2*aimag(I_0h(2:ie/2,:))
Dyn%Grid(1,1)%I_h0(1:ie/2+1,:,future) = ie**2*real(I_h0(1:ie/2+1,:))
Dyn%Grid(1,1)%I_h0(ie/2+2:ie,:,future)= ie**2*aimag(I_h0(2:ie/2,:))
Dyn%Grid(1,1)%I_00(1:ie/2+1,:,future) = ie**2*real(I_00(1:ie/2+1,:))
Dyn%Grid(1,1)%I_00(ie/2+2:ie,:,future)= ie**2*aimag(I_00(2:ie/2,:))

!I_ps = sum(sum(enstrophy_rvorprime_adveck(2:4,5:13,:,:),1),1) ! A
!I_hp = sum(sum(enstrophy_rvorprime_adveck(14:ie/2-1,2:4,:,:),1),1) ! B
!I_ph = sum(sum(enstrophy_rvorprime_adveck(2:4,14:ie/2-1,:,:),1),1) ! B
!I_ss = sum(sum(enstrophy_rvorprime_adveck(5:13,5:13,:,:),1),1)
!I_sh = sum(sum(enstrophy_rvorprime_adveck(5:13,14:ie/2-1,:,:),1),1) ! C
!I_hs = sum(sum(enstrophy_rvorprime_adveck(14:ie/2-1,5:13,:,:),1),1) ! C
!I_hh = sum(sum(enstrophy_rvorprime_adveck(14:ie/2-1,14:ie/2-1,:,:),1),1)

!print *, "SUM OVER ALL", size(sum(sum(enstrophy_rvorprime_adveck(:,:,1:5,32:64),1),1),1), size(sum(sum(enstrophy_rvorprime_adveck(:,:,1:5,32:64),1),1),2)
!print *, enstrophy_rvorprime_adveck(1:2,1:2,1:2,1:2)
!print *, "SUBDIVISION", (I_pp(1:5,32:64) + I_sp(1:5,32:64) + I_ps(1:5,32:64) + I_ph(1:5,32:64) + I_hp(1:5,32:64) + I_ss(1:5,32:64) + I_hh(1:5,32:64) + I_sh(1:5,32:64) + I_hs(1:5,32:64))
!print *, "SUM OVER ALL - SUBDIVISION", sum(sum(enstrophy_rvorprime_adveck,1),1) - (I_pp(1:5,32:64) + I_sp(1:5,32:64) + I_ps(1:5,32:64) + I_ph(1:5,32:64) + I_hp(1:5,32:64) + I_ss(1:5,32:64) & 
!			+ I_hh(1:5,32:64) + I_sh(1:5,32:64) + I_hs(1:5,32:64))

Dyn%Grid(1,1)%enstrophy_rvorprime_adveck_SUM(1:ie/2+1,:,future) = ie**2*real(sum(sum(enstrophy_rvorprime_adveck(2:ie,2:ie,1:ie/2+1,:),1),1))
Dyn%Grid(1,1)%enstrophy_rvorprime_adveck_SUM(ie/2+2:ie,:,future)= ie**2*aimag(sum(sum(enstrophy_rvorprime_adveck(2:ie,2:ie,2:ie/2,:),1),1))

!Dyn%Grid(1,1)%tester(1:ie/2+1,:,future) = ie**2*real(sum(sum(tester(2:ie,2:ie,1:ie/2+1,:),1),1))
!Dyn%Grid(1,1)%tester(ie/2+2:ie,:,future)= ie**2*aimag(sum(sum(tester(2:ie,2:ie,2:ie/2,:),1),1))

!Dyn%Grid(1,1)%tester2(1:ie/2+1,:,future) = ie**2*real(sum(sum(tester2(2:ie,2:ie,1:ie/2+1,:),1),1))
!Dyn%Grid(1,1)%tester2(ie/2+2:ie,:,future)= ie**2*aimag(sum(sum(tester2(2:ie,2:ie,2:ie/2,:),1),1))

!Dyn%Grid(1,1)%tester3(1:ie/2+1,:,future) = ie**2*real(sum(sum(tester3(2:ie,2:ie,1:ie/2+1,:),1),1))
!Dyn%Grid(1,1)%tester3(ie/2+2:ie,:,future)= ie**2*aimag(sum(sum(tester3(2:ie,2:ie,2:ie/2,:),1),1))

!print *, Dyn%Grid(1,1)%tester3(1:5,48,:)

!print *, "SUM", Dyn%Grid(1,1)%enstrophy_rvorprime_adveck_SUM(8,40:45,future)
!print *, "SUM2", Dyn%Grid(1,1)%tester(8,40:45,future)
!print *, "I_pp + I_ps + I_ph + I_sp + I_ss + I_sh + I_hp + I_hs + I_hh", Dyn%Grid(1,1)%I_pp(8,40:45,future) + Dyn%Grid(1,1)%I_ps(8,40:45,future) + Dyn%Grid(1,1)%I_ph(8,40:45,future) + &
!						    			 Dyn%Grid(1,1)%I_sp(8,40:45,future) + Dyn%Grid(1,1)%I_ss(8,40:45,future) + Dyn%Grid(1,1)%I_sh(8,40:45,future) + &
!						    			 Dyn%Grid(1,1)%I_hp(8,40:45,future) + Dyn%Grid(1,1)%I_hs(8,40:45,future) + Dyn%Grid(1,1)%I_hh(8,40:45,future)

Dyn%Grid(1,1)%I_XX(:,:,future) = Dyn%Grid(1,1)%I_pp(:,:,future) + Dyn%Grid(1,1)%I_ps(:,:,future) + Dyn%Grid(1,1)%I_ph(:,:,future) + Dyn%Grid(1,1)%I_sp(:,:,future) + &
				 Dyn%Grid(1,1)%I_ss(:,:,future) + Dyn%Grid(1,1)%I_sh(:,:,future) + Dyn%Grid(1,1)%I_hp(:,:,future) + Dyn%Grid(1,1)%I_hs(:,:,future) + Dyn%Grid(1,1)%I_hh(:,:,future)

!print *, "SUM IT ALL", sum(sum(Dyn%Grid(1,1)%I_XX(:,:,future),2),1)
!print *, "I_ph", Dyn%Grid(1,1)%I_ph(1:20,45,future)
!print *, "I_tester", Dyn%Grid(1,1)%tester(1:20,45,future)
!print *, "SUM 2", Dyn%Grid(1,1)%enstrophy_rvorprime_adveck_SUM(2,45,future)

do p=1,12!size(kmat)-1
       do q=1,12!size(kmat)-1
	    Dyn%Grid(p,q)%I1(1:ie/2+1,:,future) = ie**2*real(enstrophy_rvorprime_adveck(p+1,q+1,1:ie/2+1,:))
	    Dyn%Grid(p,q)%I1(ie/2+2:ie,:,future) = ie**2*aimag(enstrophy_rvorprime_adveck(p+1,q+1,2:ie/2,:))
       end do
end do

do p=1,1!2!size(kmat)-1
       do q=1,12!size(kmat)-1
	    Dyn%Grid(p,q)%I2(1:ie/2+1,:,future) = ie**2*real(enstrophy_rvorprime_adveck(p+1,ie-q+1,1:ie/2+1,:))
	    Dyn%Grid(p,q)%I2(ie/2+2:ie,:,future) = ie**2*aimag(enstrophy_rvorprime_adveck(p+1,ie-q+1,2:ie/2,:))
       end do
end do

do p=1,1!2!size(kmat)-1
       do q=1,12!size(kmat)-1
	    Dyn%Grid(p,q)%I3(1:ie/2+1,:,future) = ie**2*real(enstrophy_rvorprime_adveck(ie-p+1,q+1,1:ie/2+1,:))
	    Dyn%Grid(p,q)%I3(ie/2+2:ie,:,future) = ie**2*aimag(enstrophy_rvorprime_adveck(ie-p+1,q+1,2:ie/2,:))
       end do
end do

do p=1,1!2!size(kmat)
       do q=1,12!size(kmat)
	    Dyn%Grid(p,q)%I4(1:ie/2+1,:,future) = ie**2*real(enstrophy_rvorprime_adveck(ie-p+1,ie-q+1,1:ie/2+1,:))
	    Dyn%Grid(p,q)%I4(ie/2+2:ie,:,future) = ie**2*aimag(enstrophy_rvorprime_adveck(ie-p+1,ie-q+1,2:ie/2,:))
       end do
end do

Dyn%Grid(1,1)%I_higher_freq(1:ie/2+1,:,future) = ie**2*real(sum(sum(enstrophy_rvorprime_adveck(14:128-12,14:128-12,1:ie/2+1,:),1),1))
Dyn%Grid(1,1)%I_higher_freq(ie/2+2:ie,:,future) = ie**2*aimag(sum(sum(enstrophy_rvorprime_adveck(14:128-12,14:128-12,2:ie/2,:),1),1))

endif

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(minval(Dyn%Grid(1,1)%v) < valid_range_v(1) .or. maxval(Dyn%Grid(1,1)%v) > valid_range_v(2)) then
  call error_mesg('barotropic_dynamics','meridional wind out of valid range balls', FATAL)
endif

if(Dyn%spec_tracer) call update_spec_tracer(Dyn%Spec%trs, Dyn%Grid(1,1)%trs, Dyn%Tend%trs, &
                         Dyn%Grid(1,1)%u, Dyn%Grid(1,1)%v, previous, current, future, delta_t)

if(Dyn%grid_tracer) call update_grid_tracer(Dyn%Grid(1,1)%tr, Dyn%Tend%tr, &
                         Dyn%Grid(1,1)%u, Dyn%Grid(1,1)%v, previous, current, future, delta_t)

stream = compute_laplacian(Dyn%Spec%vor(:,:,current), -1) 
 call trans_spherical_to_grid(stream, Dyn%Grid(1,1)%stream)

! -------------------------------- BUDGET EQUATIONS ----------------------------------------
! CHECK TO MAKE SURE BUDGET IS BEING BALANCED CORRECTLY (RESIDUAL SHOULD BE ~<= 10e-19)
! ZONAL MOMENTUM
!print *, "BUDGET: UTEND_AFTER = VQ - UTRUNC", Dyn%Grid(1,1)%utend(40,40,future) - (Dyn%Grid(1,1)%vq(40,40,future) - Dyn%Grid(1,1)%utrunc(40,40,future)) ! BUDGET WITH NO DAMPING (SET COEFFICENTS TO ZERO IN NAMELIST)
!print *, "BUDGET: UTEND_AFTER = VQ - UTRUNC - LIND + SPECD + STIR + ROBD", Dyn%Grid(1,1)%utend(10,10,future) - (Dyn%Grid(1,1)%vq(10,10,future) - Dyn%Grid(1,1)%utrunc(10,10,future) - Dyn%Grid(1,1)%ulindamp_m(10,10,future) &
!- Dyn%Grid(1,1)%ulindamp_e(10,10,future) + Dyn%Grid(1,1)%uspecdamp(10,10,future) + Dyn%Grid(1,1)%ustir(10,10,future) + Dyn%Grid(1,1)%urobdamp(10,10,future))
! MERIDIONAL MOMENTUM
!print *, "BUDGET: VTEND = -UQ - VTRUNC", Dyn%Grid(1,1)%vtend(10,10,future) - (-Dyn%Grid(1,1)%uq(10,10,future) - Dyn%Grid(1,1)%vtrunc(10,10,future)) 
!print *, "BUDGET: VTEND = -UQ - VTRUNC + VSPECDAMP + VROBDAMP", Dyn%Grid(1,1)%vtend(27,50,future) - (-Dyn%Grid(1,1)%uq(27,50,future) - Dyn%Grid(1,1)%vtrunc(27,50,future) + Dyn%Grid(1,1)%vspecdamp(27,50,future) &
! + Dyn%Grid(1,1)%vrobdamp(27,50,future)) ! BUDGET WITH SPECTRAL DAMPING & ROBERTS DAMPING
! VORTICTY
!print *, "BUDGET: VORTEND = DT_VOR - VORTRUNC", Dyn%Grid(1,1)%vortend(22,47,future) - (Dyn%Grid(1,1)%dt_vor(22,47,future) - Dyn%Grid(1,1)%vortrunc(22,47,future))
!print *, "BUDGET: VORTEND = DT_VOR - VORTRUNC - LIND + SPECD + STIR + ROBD", Dyn%Grid(1,1)%vortend(22,47,future) - (Dyn%Grid(1,1)%dt_vor(22,47,future) - Dyn%Grid(1,1)%vortrunc(22,47,future) - & 
!Dyn%Grid(1,1)%vorlindamp_e(22,47,future) - Dyn%Grid(1,1)%vorlindamp_m(22,47,future) + Dyn%Grid(1,1)%vorspecdamp(22,47,future) + Dyn%Grid(1,1)%vorstir(22,47,future) + Dyn%Grid(1,1)%vorrobdamp(22,47,future))
!print *, "BUDGET: VORTEND = PVOR_ADV + RVOR_ADV - VORTRUNC - LIND + SPECD + STIR + ROBD", Dyn%Grid(1,1)%vortend(22,47,future) - (Dyn%Grid(1,1)%pvor_advec(22,47,future) + Dyn%Grid(1,1)%rvormean_advec(22,47,future) + &
! Dyn%Grid(1,1)%rvorprime_advec(22,47,future)- Dyn%Grid(1,1)%vortrunc(22,47,future) - & 
!Dyn%Grid(1,1)%vorlindamp_e(22,47,future) - Dyn%Grid(1,1)%vorlindamp_m(22,47,future) + Dyn%Grid(1,1)%vorspecdamp(22,47,future) + Dyn%Grid(1,1)%vorstir(22,47,future) + Dyn%Grid(1,1)%vorrobdamp(22,47,future))
! ENERGY
!print *, "ENERGY_TEND = VORADVEC - GRADTERM - TRUNC - LINDAMP + SPECDAMP + STIR + ROBDAMP", Dyn%Grid(1,1)%energy_tend(10,10,future) - (Dyn%Grid(1,1)%energy_voradvec(10,10,future) & 
! - Dyn%Grid(1,1)%energy_gradterm(10,10,future) - Dyn%Grid(1,1)%energy_lindamp_m(10,10,future) - Dyn%Grid(1,1)%energy_lindamp_e(10,10,future) - Dyn%Grid(1,1)%energy_trunc(10,10,future) + Dyn%Grid(1,1)%energy_specdamp(10,10,future)&
! + Dyn%Grid(1,1)%energy_stir(10,10,future) + Dyn%Grid(1,1)%energy_robdamp(10,10,future))
! MEAN ENERGY
!print *, "ENERGY_TEND_MEAN = VORADVEC - GRADTERM - TRUNC - LINDAMP + SPECDAMP + STIR + ROBDAMP", Dyn%Grid(1,1)%energy_tend_mean(10,10,future) - (Dyn%Grid(1,1)%energy_voradvec_mean(10,10,future) & 
! - Dyn%Grid(1,1)%energy_gradterm_mean(10,10,future) - Dyn%Grid(1,1)%energy_lindamp_m_mean(10,10,future) - Dyn%Grid(1,1)%energy_lindamp_e_mean(10,10,future) - Dyn%Grid(1,1)%energy_trunc_mean(10,10,future) & 
! + Dyn%Grid(1,1)%energy_specdamp_mean(10,10,future) + Dyn%Grid(1,1)%energy_stir_mean(10,10,future) + Dyn%Grid(1,1)%energy_robdamp_mean(10,10,future))
! ENSTROPHY
!print *, "BUDGET: ENST_TEND = PVOR_ADVEC + RVOR_ADVEC - VORTRUNC", Dyn%Grid(1,1)%enstrophy_tend(22,47,future) -(Dyn%Grid(1,1)%enstrophy_pvor_advec(22,47,future) + Dyn%Grid(1,1)%enstrophy_rvormean_advec(22,47,future) + & 
!+ Dyn%Grid(1,1)%enstrophy_rvorprime_advec(22,47,future) - Dyn%Grid(1,1)%enstrophy_trunc(22,47,future))
!print *, "BUDGET: ENS_TEND = PVORAD + RVORAD - ENSTRUNC - LIND + SPECD + STIR + ROBD", Dyn%Grid(1,1)%enstrophy_tend(22,47,future) - (Dyn%Grid(1,1)%enstrophy_pvor_advec(22,47,future) + & 
! Dyn%Grid(1,1)%enstrophy_rvormean_advec(22,47,future) + Dyn%Grid(1,1)%enstrophy_rvorprime_advec(22,47,future) - Dyn%Grid(1,1)%enstrophy_trunc(22,47,future) - Dyn%Grid(1,1)%enstrophy_lindamp_m(22,47,future) - & 
!Dyn%Grid(1,1)%enstrophy_lindamp_e(22,47,future) + Dyn%Grid(1,1)%enstrophy_specdamp(22,47,future) + Dyn%Grid(1,1)%enstrophy_stir(22,47,future) + Dyn%Grid(1,1)%enstrophy_robdamp(22,47,future) )

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
return
end subroutine barotropic_dynamics

!===================================================================================

subroutine update_spec_tracer(tr_spec, tr_grid, dt_tr, ug, vg, &
                              previous, current, future, delta_t)

complex, intent(inout), dimension(ms:me, ns:ne, num_time_levels) :: tr_spec
real   , intent(inout), dimension(is:ie, js:je, num_time_levels) :: tr_grid
real   , intent(inout), dimension(is:ie, js:je                 ) :: dt_tr
real   , intent(in   ), dimension(is:ie, js:je, num_time_levels) :: ug, vg  
real   , intent(in   )  :: delta_t
integer, intent(in   )  :: previous, current, future

complex, dimension(ms:me, ns:ne) :: dt_trs

call horizontal_advection     (tr_spec(:,:,current), ug(:,:,current), vg(:,:,current), dt_tr)
call trans_grid_to_spherical  (dt_tr, dt_trs)
call compute_spectral_damping (tr_spec(:,:,previous), dt_trs, delta_t)
call leapfrog                 (tr_spec, dt_trs, previous, current, future, delta_t, robert_coeff)
call trans_spherical_to_grid  (tr_spec(:,:,future), tr_grid(:,:,future))

return
end subroutine update_spec_tracer
!==========================================================================

subroutine update_grid_tracer(tr_grid, dt_tr_grid, ug, vg, &
                              previous, current, future, delta_t)

real   , intent(inout), dimension(is:ie, js:je, num_time_levels) :: tr_grid
real   , intent(inout), dimension(is:ie, js:je                 ) :: dt_tr_grid
real   , intent(in   ), dimension(is:ie, js:je, num_time_levels) :: ug, vg

real   , intent(in   )  :: delta_t
integer, intent(in   )  :: previous, current, future

real, dimension(size(tr_grid,1),size(tr_grid,2)) :: tr_current, tr_future

tr_future = tr_grid(:,:,previous) + delta_t*dt_tr_grid
dt_tr_grid = 0.0
call a_grid_horiz_advection (ug(:,:,current), vg(:,:,current), tr_future, delta_t, dt_tr_grid)
tr_future = tr_future + delta_t*dt_tr_grid
tr_current = tr_grid(:,:,current) + &
    robert_coeff*(tr_grid(:,:,previous) + tr_future - 2.0*tr_grid(:,:,current))
tr_grid(:,:,current) = tr_current
tr_grid(:,:,future)  = tr_future

return
end subroutine update_grid_tracer

!==========================================================================

subroutine read_restart(Dyn)

type(dynamics_type), intent(inout)  :: Dyn

integer :: unit, m, n, nt 
real, dimension(ms:me, ns:ne) :: real_part, imag_part
if(file_exist('INPUT/barotropic_dynamics.res.nc')) then
  do nt = 1, 2
    call read_data('INPUT/barotropic_dynamics.res.nc', 'vors_real', real_part, spectral_domain, timelevel=nt)
    call read_data('INPUT/barotropic_dynamics.res.nc', 'vors_imag', imag_part, spectral_domain, timelevel=nt)
    do n=ns,ne
      do m=ms,me
        Dyn%Spec%vor(m,n,nt) = cmplx(real_part(m,n),imag_part(m,n))
      end do
    end do
    call read_data('INPUT/barotropic_dynamics.res.nc', 'roberts_real', real_part, spectral_domain, timelevel=nt) ! Read restart for Dyn%spec%roberts variable
    call read_data('INPUT/barotropic_dynamics.res.nc', 'roberts_imag', imag_part, spectral_domain, timelevel=nt)
    do n=ns,ne
      do m=ms,me
        Dyn%Spec%roberts(m,n,nt) = cmplx(real_part(m,n),imag_part(m,n))
      end do
    end do
    if(Dyn%spec_tracer) then
      call read_data('INPUT/barotropic_dynamics.res.nc', 'trs_real', real_part, spectral_domain, timelevel=nt)
      call read_data('INPUT/barotropic_dynamics.res.nc', 'trs_imag', imag_part, spectral_domain, timelevel=nt)
      do n=ns,ne
        do m=ms,me
           Dyn%Spec%trs(m,n,nt) = cmplx(real_part(m,n),imag_part(m,n))
        end do
      end do
    endif
    call read_data('INPUT/barotropic_dynamics.res.nc', 'u',   Dyn%Grid(1,1)%u  (:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/barotropic_dynamics.res.nc', 'v',   Dyn%Grid(1,1)%v  (:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/barotropic_dynamics.res.nc', 'vor', Dyn%Grid(1,1)%vor(:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/barotropic_dynamics.res.nc', 'energy', Dyn%Grid(1,1)%energy(:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/barotropic_dynamics.res.nc', 'vormean_y_tend', Dyn%Grid(1,1)%vormean_y_tend(:,:,nt), grid_domain, timelevel=nt)
    if(Dyn%spec_tracer) then
      call read_data('INPUT/barotropic_dynamics.res.nc', 'trs', Dyn%Grid(1,1)%trs(:,:,nt), grid_domain, timelevel=nt)
    endif
    if(Dyn%grid_tracer) then
      call read_data('INPUT/barotropic_dynamics.res.nc', 'tr', Dyn%Grid(1,1)%tr(:,:,nt), grid_domain, timelevel=nt)
    endif
  end do


else if(file_exist('INPUT/barotropic_dynamics.res')) then
  unit = open_restart_file(file='INPUT/barotropic_dynamics.res',action='read')

  do nt = 1, 2
    call set_domain(spectral_domain)
    call read_data(unit,Dyn%Spec%vor(:,:, nt))
    if(Dyn%spec_tracer) call read_data(unit,Dyn%Spec%trs(:,:, nt))

    call set_domain(grid_domain)
    call read_data(unit,Dyn%Grid(1,1)%u   (:,:, nt))
    call read_data(unit,Dyn%Grid(1,1)%v   (:,:, nt))
    call read_data(unit,Dyn%Grid(1,1)%vor (:,:, nt))
    if(Dyn%spec_tracer) call read_data(unit,Dyn%Grid(1,1)%trs(:,:, nt))
    if(Dyn%grid_tracer) call read_data(unit,Dyn%Grid(1,1)%tr (:,:, nt))
    
  end do
  call close_file(unit)


else
  call error_mesg('read_restart', 'restart does not exist', FATAL)
endif


return
end subroutine read_restart

!====================================================================

subroutine write_restart(Dyn, previous, current)

type(dynamics_type), intent(in)  :: Dyn
integer, intent(in) :: previous, current

integer :: unit, nt, nn

do nt = 1, 2
  if(nt == 1) nn = previous
  if(nt == 2) nn = current
  call write_data('RESTART/barotropic_dynamics.res.nc', 'vors_real',  real(Dyn%Spec%vor(:,:,nn)), spectral_domain)
  call write_data('RESTART/barotropic_dynamics.res.nc', 'vors_imag', aimag(Dyn%Spec%vor(:,:,nn)), spectral_domain)
  call write_data('RESTART/barotropic_dynamics.res.nc', 'roberts_real',  real(Dyn%Spec%roberts(:,:,nn)), spectral_domain) ! Write roberts into restart file
  call write_data('RESTART/barotropic_dynamics.res.nc', 'roberts_imag', aimag(Dyn%Spec%roberts(:,:,nn)), spectral_domain)
  if(Dyn%spec_tracer) then
    call write_data('RESTART/barotropic_dynamics.res.nc', 'trs_real',  real(Dyn%Spec%trs(:,:,nn)), spectral_domain)
    call write_data('RESTART/barotropic_dynamics.res.nc', 'trs_imag', aimag(Dyn%Spec%trs(:,:,nn)), spectral_domain)
  endif
  call write_data('RESTART/barotropic_dynamics.res.nc', 'u',   Dyn%Grid(1,1)%u  (:,:,nn), grid_domain)
  call write_data('RESTART/barotropic_dynamics.res.nc', 'v',   Dyn%Grid(1,1)%v  (:,:,nn), grid_domain)
  call write_data('RESTART/barotropic_dynamics.res.nc', 'vor', Dyn%Grid(1,1)%vor(:,:,nn), grid_domain)
  call write_data('RESTART/barotropic_dynamics.res.nc', 'energy', Dyn%Grid(1,1)%energy(:,:,nn), grid_domain)
  call write_data('RESTART/barotropic_dynamics.res.nc', 'vormean_y_tend', Dyn%Grid(1,1)%vormean_y_tend(:,:,nn), grid_domain)

!  call write_data('RESTART/enstrophy_rvorprime_adveck_SUM.res.nc', 'enstrophy_rvorprime_adveck_SUM', Dyn%Grid(1,1)%enstrophy_rvorprime_adveck_SUM(:,:,nn), grid_domain)
!   call write_data('RESTART/enstrophy_rvorprime_adveck_SUM2.nc', 'enstrophy_rvorprime_adveck_SUM2', Dyn%Grid(1,1)%enstrophy_rvorprime_adveck_SUM(:,:,nn), grid_domain)

  if(Dyn%spec_tracer) then
    call write_data('RESTART/barotropic_dynamics.res.nc', 'trs', Dyn%Grid(1,1)%trs(:,:,nn), grid_domain)
  endif
  if(Dyn%grid_tracer) then
    call write_data('RESTART/barotropic_dynamics.res.nc', 'tr', Dyn%Grid(1,1)%tr(:,:,nn), grid_domain)
  endif
enddo

! unit = open_restart_file(file='RESTART/barotropic_dynamics.res', action='write')
!
!do nt = 1, 2
!  if(nt == 1) nn = previous
! if(nt == 2) nn = current
! 
!  call set_domain(spectral_domain)
!  call write_data(unit,Dyn%Spec%vor(:,:, nn))
! if(Dyn%spec_tracer) call write_data(unit,Dyn%Spec%trs(:,:, nn))
!
!  call set_domain(grid_domain)
!  call write_data(unit,Dyn%Grid(1,1)%u   (:,:, nn))
!  call write_data(unit,Dyn%Grid(1,1)%v   (:,:, nn))
!  call write_data(unit,Dyn%Grid(1,1)%vor (:,:, nn))
!  if(Dyn%spec_tracer) call write_data(unit,Dyn%Grid(1,1)%trs(:,:, nn))
!  if(Dyn%grid_tracer) call write_data(unit,Dyn%Grid(1,1)%tr (:,:, nn))
!end do

! call close_file(unit)

end subroutine write_restart

!====================================================================

subroutine barotropic_dynamics_end (Dyn, previous, current)

type(dynamics_type), intent(inout)  :: Dyn
integer, intent(in) :: previous, current

if(.not.module_is_initialized) then
  call error_mesg('barotropic_dynamics','dynamics has not been initialized ', FATAL)
endif

call transforms_end()
call stirring_end()

call write_restart (Dyn, previous, current)

module_is_initialized = .false.

return
end subroutine barotropic_dynamics_end
!===================================================================================

end module barotropic_dynamics_mod
