# Makefile created by mkmf $Id: mkmf,v 13.0 2006/04/10 21:20:01 fms Exp $ 

include /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../bin/mkmf.template.ia64

SRCROOT = /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/

CPPDEFS = -Duse_libMPI -Duse_netCDF

.DEFAULT:
	-touch $@
all: fms.x
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 fms.o constants.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/time_manager/time_manager.F90
bgrid_cold_start.o: $(SRCROOT)atmos_bgrid/tools/bgrid_cold_start.f90 bgrid_horiz.o bgrid_halo.o gaussian_topog.o fms.o constants.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_cold_start.f90
bgrid_conserve_energy.o: $(SRCROOT)atmos_bgrid/model/bgrid_conserve_energy.f90 bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o bgrid_change_grid.o time_manager.o diag_manager.o fms.o constants.o mpp_domains.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_conserve_energy.f90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_sma.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
vert_advection.o: $(SRCROOT)atmos_shared/vert_advection/vert_advection.f90 fms.o mpp.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_shared/vert_advection/vert_advection.f90
bgrid_masks.o: $(SRCROOT)atmos_bgrid/tools/bgrid_masks.f90 bgrid_horiz.o bgrid_halo.o bgrid_vert.o fms.o mpp.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_masks.f90
atmos_model.o: $(SRCROOT)atmos_solo/atmos_model.f90 atmosphere.o time_manager.o fms.o fms_io.o mpp_io.o diag_manager.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_solo/atmos_model.f90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_util.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_misc.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_connect.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_read.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_read_2Ddecomp.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_write.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write_2Ddecomp.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
bgrid_horiz.o: $(SRCROOT)atmos_bgrid/tools/bgrid_horiz.f90 mpp_domains.o constants.o fms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_horiz.f90
atmosphere.o: $(SRCROOT)atmos_bgrid/driver/shallow/atmosphere.f90 bgrid_core_driver.o bgrid_prog_var.o bgrid_horiz.o time_manager.o fms.o bgrid_change_grid.o bgrid_vert.o bgrid_halo.o shallow_physics.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/driver/shallow/atmosphere.f90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_2d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_3d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_4d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms	$(SRCROOT)shared/fms/fms_io.F90
horiz_interp_spherical.o: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.f90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_spherical.f90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
bgrid_advection.o: $(SRCROOT)atmos_bgrid/model/bgrid_advection.f90 bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o bgrid_polar_filter.o bgrid_halo.o bgrid_change_grid.o vert_advection.o fms.o mpp.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_advection.f90
topography.o: $(SRCROOT)shared/topography/topography.f90 gaussian_topog.o horiz_interp.o fms.o fms_io.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/topography/topography.f90
bgrid_sponge.o: $(SRCROOT)atmos_bgrid/model/bgrid_sponge.f90 bgrid_horiz.o bgrid_masks.o bgrid_prog_var.o bgrid_change_grid.o bgrid_halo.o fms.o mpp.o mpp_domains.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_sponge.f90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 diag_data.o diag_axis.o diag_output.o fms.o mpp_domains.o time_manager.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_util.F90
bgrid_horiz_adjust.o: $(SRCROOT)atmos_bgrid/model/bgrid_horiz_adjust.f90 bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_halo.o bgrid_change_grid.o constants.o fms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_horiz_adjust.f90
horiz_interp.o: $(SRCROOT)shared/horiz_interp/horiz_interp.f90 fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp.f90
bgrid_vert.o: $(SRCROOT)atmos_bgrid/tools/bgrid_vert.f90 constants.o fms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_vert.f90
fft.o: $(SRCROOT)shared/fft/fft.F90 platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/fft/fft.F90
shallow_physics.o: $(SRCROOT)atmos_param/shallow_physics/shallow_physics.f90 fms.o time_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_param/shallow_physics/shallow_physics.f90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
bgrid_prog_var.o: $(SRCROOT)atmos_bgrid/tools/bgrid_prog_var.f90 bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_halo.o bgrid_cold_start.o fms.o field_manager.o tracer_manager.o platform.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_prog_var.f90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
horiz_interp_bicubic.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.f90 mpp.o fms.o horiz_interp_type.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.f90
bgrid_core_driver.o: $(SRCROOT)atmos_bgrid/model/bgrid_core_driver.f90 bgrid_core.o bgrid_horiz.o bgrid_vert.o bgrid_prog_var.o bgrid_diagnostics.o bgrid_integrals.o bgrid_conserve_energy.o field_manager.o tracer_manager.o time_manager.o fms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_core_driver.f90
horiz_interp_bilinear.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.f90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.f90
bgrid_horiz_diff.o: $(SRCROOT)atmos_bgrid/model/bgrid_horiz_diff.f90 bgrid_horiz.o bgrid_masks.o bgrid_prog_var.o bgrid_halo.o bgrid_change_grid.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_horiz_diff.f90
gaussian_topog.o: $(SRCROOT)shared/topography/gaussian_topog.f90 fms.o constants.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/topography/gaussian_topog.f90
bgrid_change_grid.o: $(SRCROOT)atmos_bgrid/tools/bgrid_change_grid.f90 bgrid_horiz.o fms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_change_grid.f90
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 time_manager.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.f90 mpp_domains.o fms.o diag_data.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_axis.f90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.f90 fms.o time_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/time_manager/get_cal_time.f90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.f90 time_manager.o mpp_io.o fms.o mpp.o diag_axis.o diag_data.o diag_util.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_manager.f90
bgrid_vert_adjust.o: $(SRCROOT)atmos_bgrid/model/bgrid_vert_adjust.f90 bgrid_vert.o fms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_vert_adjust.f90
bgrid_integrals.o: $(SRCROOT)atmos_bgrid/tools/bgrid_integrals.f90 bgrid_change_grid.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o time_manager.o fms.o constants.o mpp.o mpp_io.o mpp_domains.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_integrals.f90
constants.o: $(SRCROOT)shared/constants/constants.f90
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/constants/constants.f90
bgrid_diagnostics.o: $(SRCROOT)atmos_bgrid/tools/bgrid_diagnostics.f90 bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o bgrid_change_grid.o diag_manager.o time_manager.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_diagnostics.f90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
fms.o: $(SRCROOT)shared/fms/fms.f90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/fms/fms.f90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.f90 mpp_io.o mpp_domains.o mpp.o diag_axis.o time_manager.o fms.o platform.o diag_data.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_output.f90
horiz_interp_conserve.o: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.f90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_conserve.f90
bgrid_polar_filter.o: $(SRCROOT)atmos_bgrid/tools/bgrid_polar_filter.f90 bgrid_horiz.o fft.o fms.o mpp_domains.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_polar_filter.f90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
platform.o: $(SRCROOT)shared/platform/platform.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include	$(SRCROOT)shared/platform/platform.F90
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/system_clock.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_sma.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_nocomm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_sma.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_mpi.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_mpi.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_mpi.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_nocomm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_nocomm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_nocomm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_nocomm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum_int.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
fft99.o: $(SRCROOT)shared/fft/fft99.f90 constants.o mpp.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/fft/fft99.f90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_datatype.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_util.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_define.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_misc.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_update_domains2D.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_new.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_new.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_new.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_reduce.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_reduce.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_sum.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_field.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_new.h mpp_parameter.o mpp_data.o mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include -I/home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
horiz_interp_type.o: $(SRCROOT)shared/horiz_interp/horiz_interp_type.f90 mpp.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_type.f90
bgrid_core.o: $(SRCROOT)atmos_bgrid/model/bgrid_core.f90 bgrid_prog_var.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_advection.o bgrid_horiz_diff.o bgrid_horiz_adjust.o bgrid_vert_adjust.o bgrid_polar_filter.o bgrid_halo.o bgrid_sponge.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/model/bgrid_core.f90
bgrid_halo.o: $(SRCROOT)atmos_bgrid/tools/bgrid_halo.f90 bgrid_horiz.o fms.o mpp_domains.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_bgrid/tools/bgrid_halo.f90
./bgrid_prog_var.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_prog_var.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_prog_var.f90 .
./bgrid_polar_filter.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_polar_filter.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_polar_filter.f90 .
./shallow_physics.f90: $(SRCROOT)atmos_param/shallow_physics/shallow_physics.f90
	cp $(SRCROOT)atmos_param/shallow_physics/shallow_physics.f90 .
./mpp_global_sum.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_sum.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_sum.h .
./mpp_do_update_gsm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_gsm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_gsm.h .
./parse.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/field_manager/parse.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/field_manager/parse.inc .
./mpp_io_write.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_write.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_write.inc .
./mpp_do_redistribute_new.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_new.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_new.h .
./fft99.f90: $(SRCROOT)shared/fft/fft99.f90
	cp $(SRCROOT)shared/fft/fft99.f90 .
./mpp_transmit_sma.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_sma.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_sma.h .
./mpp_comm.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm.inc .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./horiz_interp_bilinear.f90: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.f90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.f90 .
./diag_manager.f90: $(SRCROOT)shared/diag_manager/diag_manager.f90
	cp $(SRCROOT)shared/diag_manager/diag_manager.f90 .
./mpp_write.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write.h .
./mpp_reduce_caf.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_caf.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_caf.h .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./mpp_do_update_new.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_new.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_new.h .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./mpp_global_field.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_field.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_field.h .
./mpp_do_updateV_new.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_new.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_new.h .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./mpp_util_nocomm.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_nocomm.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_nocomm.inc .
./bgrid_core.f90: $(SRCROOT)atmos_bgrid/model/bgrid_core.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_core.f90 .
./fms.f90: $(SRCROOT)shared/fms/fms.f90
	cp $(SRCROOT)shared/fms/fms.f90 .
./topography.f90: $(SRCROOT)shared/topography/topography.f90
	cp $(SRCROOT)shared/topography/topography.f90 .
./bgrid_masks.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_masks.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_masks.f90 .
./mpp_io_connect.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_connect.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_connect.inc .
./bgrid_horiz_adjust.f90: $(SRCROOT)atmos_bgrid/model/bgrid_horiz_adjust.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_horiz_adjust.f90 .
./mpp_comm_mpi.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_mpi.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_mpi.inc .
./mpp_do_updateV_caf.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_caf.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_caf.h .
./mpp_do_redistribute_caf.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_caf.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_caf.h .
./mpp_io_misc.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_misc.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_misc.inc .
./mpp_transmit_mpi.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_mpi.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_mpi.h .
./bgrid_advection.f90: $(SRCROOT)atmos_bgrid/model/bgrid_advection.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_advection.f90 .
./mpp_reduce_mpi.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_mpi.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_mpi.h .
./atmos_model.f90: $(SRCROOT)atmos_solo/atmos_model.f90
	cp $(SRCROOT)atmos_solo/atmos_model.f90 .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./system_clock.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/system_clock.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/system_clock.h .
./mpp_sum_mpi.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_mpi.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_mpi.h .
./horiz_interp_type.f90: $(SRCROOT)shared/horiz_interp/horiz_interp_type.f90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_type.f90 .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./mpp_sum_caf.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_caf.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_caf.h .
./mpp_domains_comm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.h .
./bgrid_diagnostics.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_diagnostics.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_diagnostics.f90 .
./mpp_domains_util.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_util.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_util.inc .
./fft.F90: $(SRCROOT)shared/fft/fft.F90
	cp $(SRCROOT)shared/fft/fft.F90 .
./mpp_domains_reduce.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_reduce.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_reduce.inc .
./mpp_transmit.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit.inc .
./bgrid_integrals.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_integrals.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_integrals.f90 .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./mpp_comm_sma.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_sma.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_sma.inc .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./bgrid_horiz.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_horiz.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_horiz.f90 .
./get_cal_time.f90: $(SRCROOT)shared/time_manager/get_cal_time.f90
	cp $(SRCROOT)shared/time_manager/get_cal_time.f90 .
./mpp_domains_comm.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.inc .
./mpp_do_update_old.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_old.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_old.h .
./mpp_do_redistribute_old.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_old.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_old.h .
./mpp_data_sma.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_sma.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_sma.inc .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./horiz_interp.f90: $(SRCROOT)shared/horiz_interp/horiz_interp.f90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp.f90 .
./mpp_reduce_sma.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_sma.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_sma.h .
./mpp_io_util.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_util.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_util.inc .
./atmosphere.f90: $(SRCROOT)atmos_bgrid/driver/shallow/atmosphere.f90
	cp $(SRCROOT)atmos_bgrid/driver/shallow/atmosphere.f90 .
./bgrid_sponge.f90: $(SRCROOT)atmos_bgrid/model/bgrid_sponge.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_sponge.f90 .
./bgrid_vert.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_vert.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_vert.f90 .
./vert_advection.f90: $(SRCROOT)atmos_shared/vert_advection/vert_advection.f90
	cp $(SRCROOT)atmos_shared/vert_advection/vert_advection.f90 .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./diag_output.f90: $(SRCROOT)shared/diag_manager/diag_output.f90
	cp $(SRCROOT)shared/diag_manager/diag_output.f90 .
./horiz_interp_conserve.f90: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.f90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.f90 .
./bgrid_vert_adjust.f90: $(SRCROOT)atmos_bgrid/model/bgrid_vert_adjust.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_vert_adjust.f90 .
./mpp_data_mpi.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_mpi.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_mpi.inc .
./mpp_comm_nocomm.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_nocomm.inc .
./bgrid_halo.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_halo.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_halo.f90 .
./constants.f90: $(SRCROOT)shared/constants/constants.f90
	cp $(SRCROOT)shared/constants/constants.f90 .
./fms_platform.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h .
./read_data_3d.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_3d.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_3d.inc .
./mpp_domains_misc.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_misc.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_misc.inc .
./mpp_sum_sma.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_sma.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_sma.h .
./read_data_4d.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_4d.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_4d.inc .
./bgrid_conserve_energy.f90: $(SRCROOT)atmos_bgrid/model/bgrid_conserve_energy.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_conserve_energy.f90 .
./mpp_do_global_field_caf.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_caf.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_caf.h .
./mpp_do_global_field_old.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_old.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_old.h .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./mpp_do_global_field_gsm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_gsm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_gsm.h .
./mpp_sum.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum.inc .
./horiz_interp_spherical.f90: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.f90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.f90 .
./mpp_do_updateV_gsm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_gsm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_gsm.h .
./mpp_read_2Ddecomp.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_read_2Ddecomp.h .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./mpp_sum_gsm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_gsm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_gsm.h .
./mpp_transmit_nocomm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_nocomm.h .
./bgrid_horiz_diff.f90: $(SRCROOT)atmos_bgrid/model/bgrid_horiz_diff.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_horiz_diff.f90 .
./diag_axis.f90: $(SRCROOT)shared/diag_manager/diag_axis.f90
	cp $(SRCROOT)shared/diag_manager/diag_axis.f90 .
./mpp_do_updateV_old.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_old.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_old.h .
./mpp_chksum.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum.h .
./gaussian_topog.f90: $(SRCROOT)shared/topography/gaussian_topog.f90
	cp $(SRCROOT)shared/topography/gaussian_topog.f90 .
./bgrid_change_grid.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_change_grid.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_change_grid.f90 .
./write_data.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/write_data.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/write_data.inc .
./mpp_io_read.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_read.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_read.inc .
./mpp_update_domains2D.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_update_domains2D.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_update_domains2D.h .
./mpp_reduce_gsm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_gsm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_gsm.h .
./mpp_sum_nocomm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_nocomm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_nocomm.h .
./read_data_2d.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_2d.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_2d.inc .
./mpp_global_reduce.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_reduce.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_reduce.h .
./mpp_write_2Ddecomp.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write_2Ddecomp.h .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./mpp_util.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util.inc .
./bgrid_core_driver.f90: $(SRCROOT)atmos_bgrid/model/bgrid_core_driver.f90
	cp $(SRCROOT)atmos_bgrid/model/bgrid_core_driver.f90 .
./mpp_do_update_caf.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_caf.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_caf.h .
./mpp_domains_define.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_define.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_define.inc .
./mpp_chksum_int.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum_int.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum_int.h .
./horiz_interp_bicubic.f90: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.f90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.f90 .
./mpp_do_global_field_new.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_new.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_new.h .
./bgrid_cold_start.f90: $(SRCROOT)atmos_bgrid/tools/bgrid_cold_start.f90
	cp $(SRCROOT)atmos_bgrid/tools/bgrid_cold_start.f90 .
./mpp_util_mpi.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_mpi.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_mpi.inc .
./mpp_reduce_nocomm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_nocomm.h .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./mpp_do_redistribute_gsm.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_gsm.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_gsm.h .
./mpp_datatype.h: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_datatype.h
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_datatype.h .
./mpp_util_sma.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_sma.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_sma.inc .
./mpp_data_nocomm.inc: /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_nocomm.inc
	cp /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_nocomm.inc .
SRC = $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_bgrid/tools/bgrid_cold_start.f90 $(SRCROOT)atmos_bgrid/model/bgrid_conserve_energy.f90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)atmos_shared/vert_advection/vert_advection.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_masks.f90 $(SRCROOT)atmos_solo/atmos_model.f90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)atmos_bgrid/tools/bgrid_horiz.f90 $(SRCROOT)atmos_bgrid/driver/shallow/atmosphere.f90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.f90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)atmos_bgrid/model/bgrid_advection.f90 $(SRCROOT)shared/topography/topography.f90 $(SRCROOT)atmos_bgrid/model/bgrid_sponge.f90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)atmos_bgrid/model/bgrid_horiz_adjust.f90 $(SRCROOT)shared/horiz_interp/horiz_interp.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_vert.f90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)atmos_param/shallow_physics/shallow_physics.f90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)atmos_bgrid/tools/bgrid_prog_var.f90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.f90 $(SRCROOT)atmos_bgrid/model/bgrid_core_driver.f90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.f90 $(SRCROOT)atmos_bgrid/model/bgrid_horiz_diff.f90 $(SRCROOT)shared/topography/gaussian_topog.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_change_grid.f90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/diag_manager/diag_axis.f90 $(SRCROOT)shared/time_manager/get_cal_time.f90 $(SRCROOT)shared/diag_manager/diag_manager.f90 $(SRCROOT)atmos_bgrid/model/bgrid_vert_adjust.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_integrals.f90 $(SRCROOT)shared/constants/constants.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_diagnostics.f90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/fms/fms.f90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/diag_manager/diag_output.f90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_polar_filter.f90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/fft/fft99.f90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.f90 $(SRCROOT)atmos_bgrid/model/bgrid_core.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_halo.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_sum.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/field_manager/parse.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_write.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_new.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_new.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_field.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_new.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_nocomm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_connect.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_misc.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_mpi.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_mpi.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/system_clock.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_mpi.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_util.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_reduce.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_sma.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_sma.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_util.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_nocomm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_3d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_misc.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_4d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_read_2Ddecomp.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_nocomm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/write_data.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_read.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_update_domains2D.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_nocomm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_2d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write_2Ddecomp.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_reduce.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum_int.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_define.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_new.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_nocomm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_sma.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_nocomm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_datatype.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_gsm.h
OBJ = time_manager.o bgrid_cold_start.o bgrid_conserve_energy.o mpp_data.o vert_advection.o bgrid_masks.o atmos_model.o mpp_io.o bgrid_horiz.o atmosphere.o fms_io.o horiz_interp_spherical.o tracer_manager.o bgrid_advection.o topography.o bgrid_sponge.o diag_util.o bgrid_horiz_adjust.o horiz_interp.o bgrid_vert.o fft.o shallow_physics.o memuse.o memutils.o bgrid_prog_var.o threadloc.o horiz_interp_bicubic.o bgrid_core_driver.o horiz_interp_bilinear.o bgrid_horiz_diff.o gaussian_topog.o bgrid_change_grid.o diag_data.o diag_axis.o get_cal_time.o diag_manager.o bgrid_vert_adjust.o bgrid_integrals.o constants.o bgrid_diagnostics.o mpp_pset.o fms.o field_manager.o diag_output.o horiz_interp_conserve.o bgrid_polar_filter.o nsclock.o mpp_parameter.o platform.o mpp.o fft99.o mpp_domains.o horiz_interp_type.o bgrid_core.o bgrid_halo.o
OFF = $(SRCROOT)atmos_bgrid/tools/bgrid_prog_var.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_polar_filter.f90 $(SRCROOT)atmos_param/shallow_physics/shallow_physics.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_sum.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/field_manager/parse.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_write.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_new.h $(SRCROOT)shared/fft/fft99.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm.inc $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.f90 $(SRCROOT)shared/diag_manager/diag_manager.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_caf.h $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)shared/mpp/mpp_pset.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_new.h $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/fms/fms_io.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_field.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_new.h $(SRCROOT)shared/memutils/memuse.c /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_nocomm.inc $(SRCROOT)atmos_bgrid/model/bgrid_core.f90 $(SRCROOT)shared/fms/fms.f90 $(SRCROOT)shared/topography/topography.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_masks.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_connect.inc $(SRCROOT)atmos_bgrid/model/bgrid_horiz_adjust.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_misc.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_mpi.h $(SRCROOT)atmos_bgrid/model/bgrid_advection.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)atmos_solo/atmos_model.f90 $(SRCROOT)shared/platform/platform.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/system_clock.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)shared/horiz_interp/horiz_interp_type.f90 $(SRCROOT)shared/mpp/mpp_domains.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.h $(SRCROOT)atmos_bgrid/tools/bgrid_diagnostics.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_util.inc $(SRCROOT)shared/fft/fft.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_reduce.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit.inc $(SRCROOT)atmos_bgrid/tools/bgrid_integrals.f90 $(SRCROOT)shared/diag_manager/diag_data.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_sma.inc $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)atmos_bgrid/tools/bgrid_horiz.f90 $(SRCROOT)shared/time_manager/get_cal_time.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_comm.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_sma.inc $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_util.inc $(SRCROOT)atmos_bgrid/driver/shallow/atmosphere.f90 $(SRCROOT)atmos_bgrid/model/bgrid_sponge.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_vert.f90 $(SRCROOT)atmos_shared/vert_advection/vert_advection.f90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/diag_manager/diag_output.f90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.f90 $(SRCROOT)atmos_bgrid/model/bgrid_vert_adjust.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)atmos_bgrid/tools/bgrid_halo.f90 $(SRCROOT)shared/constants/constants.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/include/fms_platform.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_3d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_misc.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_sma.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_4d.inc $(SRCROOT)atmos_bgrid/model/bgrid_conserve_energy.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_old.h $(SRCROOT)shared/mpp/threadloc.c /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum.inc $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)shared/time_manager/time_manager.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)atmos_bgrid/model/bgrid_horiz_diff.f90 $(SRCROOT)shared/diag_manager/diag_axis.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_updateV_old.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum.h $(SRCROOT)shared/topography/gaussian_topog.f90 $(SRCROOT)atmos_bgrid/tools/bgrid_change_grid.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/write_data.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_io_read.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_update_domains2D.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_sum_nocomm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/fms/read_data_2d.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_global_reduce.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)shared/mpp/nsclock.c /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util.inc $(SRCROOT)atmos_bgrid/model/bgrid_core_driver.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_update_caf.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_domains_define.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_chksum_int.h $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_global_field_new.h $(SRCROOT)atmos_bgrid/tools/bgrid_cold_start.f90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_mpi.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)shared/diag_manager/diag_util.F90 /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_do_redistribute_gsm.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_datatype.h /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_util_sma.inc /home/z3309205/atm_dycores/atm_dycores/exp/bgrid_shallow/../../src/shared/mpp/include/mpp_data_nocomm.inc
clean: neat
	-rm -f .cppdefs $(OBJ) fms.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
fms.x: $(OBJ)
	$(LD) $(OBJ) -o fms.x $(LDFLAGS)
