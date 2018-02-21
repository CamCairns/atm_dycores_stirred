# Makefile created by mkmf $Id: mkmf,v 13.0 2006/04/10 21:20:01 fms Exp $ 

include /home/cwc210/FMS_models/atm_dycores_stirred//bin/mkmf.template.ia64

SRCROOT = /home/cwc210/FMS_models/atm_dycores_stirred//src/

CPPDEFS = -Duse_libMPI -Duse_netCDF

.DEFAULT:
	-touch $@
all: fms.x
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 fms.o constants.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/time_manager/time_manager.F90
spherical_fourier.o: $(SRCROOT)atmos_spectral/tools/spherical_fourier.f90 fms.o mpp.o mpp_domains.o constants.o spec_mpp.o spherical.o gauss_and_legendre.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical_fourier.f90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_mpi.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
atmos_model.o: $(SRCROOT)atmos_solo/atmos_model.f90 atmosphere.o time_manager.o fms.o fms_io.o mpp_io.o diag_manager.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_solo/atmos_model.f90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h /sw/netcdf-4.3.0/include/netcdf.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_misc.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_connect.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_read.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_read_2Ddecomp.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_write.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write_2Ddecomp.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include -I/sw/netcdf-4.3.0/include -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
atmosphere.o: $(SRCROOT)atmos_spectral_barotropic/atmosphere.f90 fms.o transforms.o time_manager.o barotropic_dynamics.o barotropic_physics.o diag_manager.o barotropic_diagnostics.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral_barotropic/atmosphere.f90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_2d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_3d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_4d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms	$(SRCROOT)shared/fms/fms_io.F90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 diag_data.o diag_axis.o diag_output.o fms.o mpp_domains.o time_manager.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_util.F90
barotropic_physics.o: $(SRCROOT)atmos_spectral_barotropic/barotropic_physics.f90 fms.o transforms.o time_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral_barotropic/barotropic_physics.f90
spectral_damping.o: $(SRCROOT)atmos_spectral/model/spectral_damping.f90 fms.o transforms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/model/spectral_damping.f90
barotropic_diagnostics.o: $(SRCROOT)atmos_spectral_barotropic/barotropic_diagnostics.f90 fms.o transforms.o diag_manager.o time_manager.o barotropic_physics.o barotropic_dynamics.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral_barotropic/barotropic_diagnostics.f90
fft.o: $(SRCROOT)shared/fft/fft.F90 platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/fft/fft.F90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 /sw/netcdf-4.3.0/include/netcdf.inc time_manager.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/sw/netcdf-4.3.0/include	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.f90 mpp_domains.o fms.o diag_data.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_axis.f90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.f90 fms.o time_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/time_manager/get_cal_time.f90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.f90 time_manager.o mpp_io.o fms.o mpp.o diag_axis.o diag_data.o diag_util.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_manager.f90
transforms.o: $(SRCROOT)atmos_spectral/tools/transforms.f90 fms.o mpp.o mpp_domains.o spec_mpp.o constants.o spherical_fourier.o grid_fourier.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/tools/transforms.f90
spherical.o: $(SRCROOT)atmos_spectral/tools/spherical.f90 fms.o spec_mpp.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical.f90
constants.o: $(SRCROOT)shared/constants/constants.f90
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/constants/constants.f90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
fms.o: $(SRCROOT)shared/fms/fms.f90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/fms/fms.f90
spec_mpp.o: $(SRCROOT)atmos_spectral/tools/spec_mpp.f90 fms.o mpp_domains.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spec_mpp.f90
grid_fourier.o: $(SRCROOT)atmos_spectral/tools/grid_fourier.f90 fms.o constants.o fft.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/tools/grid_fourier.f90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
barotropic_dynamics.o: $(SRCROOT)atmos_spectral_barotropic/barotropic_dynamics.f90 fms.o time_manager.o constants.o transforms.o spectral_damping.o leapfrog.o fv_advection.o stirring.o fft.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral_barotropic/barotropic_dynamics.f90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.f90 mpp_io.o mpp_domains.o mpp.o diag_axis.o time_manager.o fms.o platform.o diag_data.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_output.f90
gauss_and_legendre.o: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.f90 fms.o constants.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/tools/gauss_and_legendre.f90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
leapfrog.o: $(SRCROOT)atmos_spectral/model/leapfrog.f90 fms.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/model/leapfrog.f90
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
stirring.o: $(SRCROOT)atmos_spectral/model/stirring.f90 constants.o time_manager.o fms.o transforms.o diag_manager.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/model/stirring.f90
platform.o: $(SRCROOT)shared/platform/platform.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include	$(SRCROOT)shared/platform/platform.F90
fv_advection.o: $(SRCROOT)atmos_spectral/model/fv_advection.f90 fms.o constants.o mpp_domains.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)atmos_spectral/model/fv_advection.f90
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/system_clock.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_mpi.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_nocomm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_mpi.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_nocomm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_nocomm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_nocomm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_nocomm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum_int.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
fft99.o: $(SRCROOT)shared/fft/fft99.f90 constants.o mpp.o
	$(FC) $(FFLAGS) -c	$(SRCROOT)shared/fft/fft99.f90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_datatype.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_define.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_misc.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_update_domains2D.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_reduce.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_reduce.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_sum.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_field.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_new.h mpp_parameter.o mpp_data.o mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include -I/home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
./spectral_damping.f90: $(SRCROOT)atmos_spectral/model/spectral_damping.f90
	cp $(SRCROOT)atmos_spectral/model/spectral_damping.f90 .
./barotropic_dynamics.f90: $(SRCROOT)atmos_spectral_barotropic/barotropic_dynamics.f90
	cp $(SRCROOT)atmos_spectral_barotropic/barotropic_dynamics.f90 .
./mpp_datatype.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_datatype.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_datatype.h .
./mpp_do_update_new.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_new.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_new.h .
./mpp_do_redistribute_caf.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_caf.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_caf.h .
./mpp_do_redistribute_old.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_old.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_old.h .
./fft99.f90: $(SRCROOT)shared/fft/fft99.f90
	cp $(SRCROOT)shared/fft/fft99.f90 .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./diag_manager.f90: $(SRCROOT)shared/diag_manager/diag_manager.f90
	cp $(SRCROOT)shared/diag_manager/diag_manager.f90 .
./mpp_do_updateV_old.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_old.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_old.h .
./parse.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/field_manager/parse.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/field_manager/parse.inc .
./mpp_global_sum.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_sum.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_sum.h .
./mpp_update_domains2D.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_update_domains2D.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_update_domains2D.h .
./grid_fourier.f90: $(SRCROOT)atmos_spectral/tools/grid_fourier.f90
	cp $(SRCROOT)atmos_spectral/tools/grid_fourier.f90 .
./spherical_fourier.f90: $(SRCROOT)atmos_spectral/tools/spherical_fourier.f90
	cp $(SRCROOT)atmos_spectral/tools/spherical_fourier.f90 .
./fms_platform.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./mpp_chksum.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum.h .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./mpp_comm_mpi.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_mpi.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_mpi.inc .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./mpp_global_field.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_field.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_field.h .
./mpp_do_update_gsm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_gsm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_gsm.h .
./mpp_reduce_mpi.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_mpi.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_mpi.h .
./mpp_comm.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm.inc .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./mpp_sum_nocomm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_nocomm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_nocomm.h .
./spherical.f90: $(SRCROOT)atmos_spectral/tools/spherical.f90
	cp $(SRCROOT)atmos_spectral/tools/spherical.f90 .
./mpp_do_global_field_caf.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_caf.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_caf.h .
./fv_advection.f90: $(SRCROOT)atmos_spectral/model/fv_advection.f90
	cp $(SRCROOT)atmos_spectral/model/fv_advection.f90 .
./stirring.f90: $(SRCROOT)atmos_spectral/model/stirring.f90
	cp $(SRCROOT)atmos_spectral/model/stirring.f90 .
./fms.f90: $(SRCROOT)shared/fms/fms.f90
	cp $(SRCROOT)shared/fms/fms.f90 .
./mpp_io_misc.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_misc.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_misc.inc .
./mpp_io_read.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_read.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_read.inc .
./mpp_do_updateV_gsm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_gsm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_gsm.h .
./atmosphere.f90: $(SRCROOT)atmos_spectral_barotropic/atmosphere.f90
	cp $(SRCROOT)atmos_spectral_barotropic/atmosphere.f90 .
./read_data_3d.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_3d.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_3d.inc .
./mpp_do_redistribute_new.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_new.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_new.h .
./mpp_util_nocomm.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_nocomm.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_nocomm.inc .
./atmos_model.f90: $(SRCROOT)atmos_solo/atmos_model.f90
	cp $(SRCROOT)atmos_solo/atmos_model.f90 .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./mpp_util.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util.inc .
./mpp_sum.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum.inc .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./mpp_transmit_mpi.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_mpi.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_mpi.h .
./mpp_io_util.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_util.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_util.inc .
./fft.F90: $(SRCROOT)shared/fft/fft.F90
	cp $(SRCROOT)shared/fft/fft.F90 .
./read_data_4d.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_4d.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_4d.inc .
./mpp_sum_caf.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_caf.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_caf.h .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./mpp_domains_comm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.h .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./get_cal_time.f90: $(SRCROOT)shared/time_manager/get_cal_time.f90
	cp $(SRCROOT)shared/time_manager/get_cal_time.f90 .
./mpp_write.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write.h .
./mpp_domains_comm.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.inc .
./mpp_write_2Ddecomp.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write_2Ddecomp.h .
./mpp_data_nocomm.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_nocomm.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_nocomm.inc .
./mpp_util_mpi.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_mpi.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_mpi.inc .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./mpp_do_global_field_old.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_old.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_old.h .
./spec_mpp.f90: $(SRCROOT)atmos_spectral/tools/spec_mpp.f90
	cp $(SRCROOT)atmos_spectral/tools/spec_mpp.f90 .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./diag_output.f90: $(SRCROOT)shared/diag_manager/diag_output.f90
	cp $(SRCROOT)shared/diag_manager/diag_output.f90 .
./mpp_data_sma.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_sma.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_sma.inc .
./mpp_io_write.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_write.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_write.inc .
./transforms.f90: $(SRCROOT)atmos_spectral/tools/transforms.f90
	cp $(SRCROOT)atmos_spectral/tools/transforms.f90 .
./mpp_io_connect.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_connect.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_connect.inc .
./mpp_sum_gsm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_gsm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_gsm.h .
./mpp_util_sma.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_sma.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_sma.inc .
./mpp_do_updateV_new.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_new.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_new.h .
./mpp_sum_mpi.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_mpi.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_mpi.h .
./mpp_do_updateV_caf.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_caf.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_caf.h .
./constants.f90: $(SRCROOT)shared/constants/constants.f90
	cp $(SRCROOT)shared/constants/constants.f90 .
./mpp_do_redistribute_gsm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_gsm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_gsm.h .
./mpp_reduce_gsm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_gsm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_gsm.h .
./read_data_2d.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_2d.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_2d.inc .
./mpp_reduce_nocomm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_nocomm.h .
./netcdf.inc: /sw/netcdf-4.3.0/include/netcdf.inc
	cp /sw/netcdf-4.3.0/include/netcdf.inc .
./mpp_domains_reduce.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_reduce.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_reduce.inc .
./mpp_do_global_field_gsm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_gsm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_gsm.h .
./mpp_domains_define.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_define.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_define.inc .
./system_clock.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/system_clock.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/system_clock.h .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./leapfrog.f90: $(SRCROOT)atmos_spectral/model/leapfrog.f90
	cp $(SRCROOT)atmos_spectral/model/leapfrog.f90 .
./mpp_transmit.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit.inc .
./mpp_do_global_field_new.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_new.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_new.h .
./mpp_data_mpi.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_mpi.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_mpi.inc .
./mpp_global_reduce.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_reduce.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_reduce.h .
./diag_axis.f90: $(SRCROOT)shared/diag_manager/diag_axis.f90
	cp $(SRCROOT)shared/diag_manager/diag_axis.f90 .
./mpp_do_update_caf.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_caf.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_caf.h .
./mpp_transmit_nocomm.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_nocomm.h .
./gauss_and_legendre.f90: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.f90
	cp $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.f90 .
./mpp_domains_misc.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_misc.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_misc.inc .
./barotropic_diagnostics.f90: $(SRCROOT)atmos_spectral_barotropic/barotropic_diagnostics.f90
	cp $(SRCROOT)atmos_spectral_barotropic/barotropic_diagnostics.f90 .
./write_data.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/write_data.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/write_data.inc .
./mpp_read_2Ddecomp.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_read_2Ddecomp.h .
./barotropic_physics.f90: $(SRCROOT)atmos_spectral_barotropic/barotropic_physics.f90
	cp $(SRCROOT)atmos_spectral_barotropic/barotropic_physics.f90 .
./mpp_domains_util.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_util.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_util.inc .
./mpp_chksum_int.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum_int.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum_int.h .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./mpp_comm_nocomm.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_nocomm.inc .
./mpp_reduce_sma.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_sma.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_sma.h .
./mpp_do_update_old.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_old.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_old.h .
./mpp_sum_sma.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_sma.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_sma.h .
./mpp_transmit_sma.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_sma.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_sma.h .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./mpp_comm_sma.inc: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_sma.inc
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_sma.inc .
./mpp_reduce_caf.h: /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_caf.h
	cp /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_caf.h .
SRC = $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.f90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)atmos_solo/atmos_model.f90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)atmos_spectral_barotropic/atmosphere.f90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)atmos_spectral_barotropic/barotropic_physics.f90 $(SRCROOT)atmos_spectral/model/spectral_damping.f90 $(SRCROOT)atmos_spectral_barotropic/barotropic_diagnostics.f90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/diag_manager/diag_axis.f90 $(SRCROOT)shared/time_manager/get_cal_time.f90 $(SRCROOT)shared/diag_manager/diag_manager.f90 $(SRCROOT)atmos_spectral/tools/transforms.f90 $(SRCROOT)atmos_spectral/tools/spherical.f90 $(SRCROOT)shared/constants/constants.f90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/fms/fms.f90 $(SRCROOT)atmos_spectral/tools/spec_mpp.f90 $(SRCROOT)atmos_spectral/tools/grid_fourier.f90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)atmos_spectral_barotropic/barotropic_dynamics.f90 $(SRCROOT)shared/diag_manager/diag_output.f90 $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.f90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)atmos_spectral/model/leapfrog.f90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)atmos_spectral/model/stirring.f90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)atmos_spectral/model/fv_advection.f90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/fft/fft99.f90 $(SRCROOT)shared/mpp/mpp_domains.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_datatype.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/field_manager/parse.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_sum.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_update_domains2D.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_mpi.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_field.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_nocomm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_misc.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_read.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_3d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_nocomm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_4d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write_2Ddecomp.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_nocomm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_mpi.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_write.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_connect.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_2d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_nocomm.h /sw/netcdf-4.3.0/include/netcdf.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_reduce.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_define.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/system_clock.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_mpi.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_reduce.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_nocomm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/write_data.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_misc.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_read_2Ddecomp.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum_int.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_nocomm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_caf.h
OBJ = time_manager.o spherical_fourier.o mpp_data.o atmos_model.o mpp_io.o atmosphere.o fms_io.o tracer_manager.o diag_util.o barotropic_physics.o spectral_damping.o barotropic_diagnostics.o fft.o memuse.o memutils.o threadloc.o diag_data.o diag_axis.o get_cal_time.o diag_manager.o transforms.o spherical.o constants.o mpp_pset.o fms.o spec_mpp.o grid_fourier.o field_manager.o barotropic_dynamics.o diag_output.o gauss_and_legendre.o nsclock.o leapfrog.o mpp_parameter.o stirring.o platform.o fv_advection.o mpp.o fft99.o mpp_domains.o
OFF = $(SRCROOT)atmos_spectral/model/spectral_damping.f90 $(SRCROOT)atmos_spectral_barotropic/barotropic_dynamics.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_datatype.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_old.h $(SRCROOT)shared/fft/fft99.f90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/diag_manager/diag_manager.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/field_manager/parse.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_sum.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_update_domains2D.h $(SRCROOT)atmos_spectral/tools/grid_fourier.f90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/include/fms_platform.h $(SRCROOT)shared/mpp/mpp_io.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum.h $(SRCROOT)shared/mpp/mpp_pset.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_mpi.inc $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/fms/fms_io.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_field.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm.inc $(SRCROOT)shared/memutils/memuse.c /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)atmos_spectral/tools/spherical.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_caf.h $(SRCROOT)atmos_spectral/model/fv_advection.f90 $(SRCROOT)atmos_spectral/model/stirring.f90 $(SRCROOT)shared/fms/fms.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_misc.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_read.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_gsm.h $(SRCROOT)atmos_spectral_barotropic/atmosphere.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_3d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_nocomm.inc $(SRCROOT)atmos_solo/atmos_model.f90 $(SRCROOT)shared/platform/platform.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum.inc $(SRCROOT)shared/mpp/mpp_domains.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_util.inc $(SRCROOT)shared/fft/fft.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_4d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_caf.h $(SRCROOT)shared/diag_manager/diag_data.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.h $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/time_manager/get_cal_time.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_comm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_write_2Ddecomp.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_nocomm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_mpi.inc $(SRCROOT)shared/memutils/memutils.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_old.h $(SRCROOT)atmos_spectral/tools/spec_mpp.f90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/diag_manager/diag_output.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_write.inc $(SRCROOT)atmos_spectral/tools/transforms.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_io_connect.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_util_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_mpi.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_updateV_caf.h $(SRCROOT)shared/constants/constants.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_redistribute_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/read_data_2d.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_nocomm.h /sw/netcdf-4.3.0/include/netcdf.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_reduce.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_gsm.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_define.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/system_clock.h $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_spectral/model/leapfrog.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_global_field_new.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_data_mpi.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_global_reduce.h $(SRCROOT)shared/diag_manager/diag_axis.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_caf.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)atmos_spectral_barotropic/barotropic_diagnostics.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/fms/write_data.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)atmos_spectral_barotropic/barotropic_physics.f90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_domains_util.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_chksum_int.h $(SRCROOT)shared/mpp/nsclock.c /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_nocomm.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_do_update_old.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_sum_sma.h /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_transmit_sma.h $(SRCROOT)shared/diag_manager/diag_util.F90 /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_comm_sma.inc /home/cwc210/FMS_models/atm_dycores_stirred//src/shared/mpp/include/mpp_reduce_caf.h
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
