OPTIONS= /NOLOGO
COPTIONS= /03

OBJFILES=energy_2D.obj recover_list_2D.obj \
	ini_divide_2D.obj keep_list_2D.obj \
	SPHYSICS_2D.obj getdata_2D.obj \
	check_limits_2D.obj \
	divide_2D.obj \
	movingObjects_2D.obj movingGate_2D.obj \
	movingPaddle_2D.obj movingWedge_2D.obj \
	updateNormals_2D.obj vorticity_2D.obj\
	periodicityCorrection_2D.obj \
	ac_NONE_2D.obj \
	kernel_correction_NC_2D.obj \
	ac_2D.obj \
	poute_2D.obj \
	gradients_calc_basic_2D.obj \
	self_BC_Dalrymple_2D.obj \
	celij_BC_Dalrymple_2D.obj \
	rigid_body_motion_2D.obj \
	variable_time_step_2D.obj \
	viscosity_artificial_2D.obj \
	correct_2D.obj \
	kernel_wendland5_2D.obj \
	EoS_Tait_2D.obj \
	densityFilter_MLS_2D.obj \
	ac_MLS_2D.obj \
	LU_decomposition_2D.obj \
	pre_celij_MLS_2D.obj \
	pre_self_MLS_2D.obj \
	step_predictor_corrector_2D.obj 
.f.obj:
	ifort $(OPTIONS) $(COPTIONS) /O3 /c $<

SPHYSICS_2D.exe: $(OBJFILES)
	xilink /OUT:$@ $(OPTIONS) $(OBJFILES)

clean:
	del *.mod *.obj
