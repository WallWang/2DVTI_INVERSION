#include "mpi.h"
#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include "iofile.h"
#include "RaySub.h"
#include "bsplines.h"
#include "sparse_matrix.h"
#include "stereo_subroutines.h"
#include "spqr_solve_matrix.h"
#include "regu.h"
#include "frechet.h"
#include "forward_modeling.h"
#include "data_fitting_err.h"
#include "update_model.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" 								",
" 								",
" Required parameters:                                          ",
" invfile=              kinematic invariant file          ",
" datadir=              name of work area(char type)          ",
" modelfile=            initial modelspace                ",
" 								",
" Optional Parameters:                                                  ",
" datatype=             data type(binary/text), 0 for text and 1 for binary          ",
" 								",
" vel_min=              the minimum velocity bound by geology          ",
" vel_max=              the maximum velocity bound by geology          ",
" 								",
" vel_begin=1000        initial value of velocity gradient         ",
" vel_end=4000          final value of velocity gradient          ",

" epos_min=             the minimum epsilon bound by geology          ",
" epos_max=             the maximum epsilon bound by geology          ",
" 								",
" epos_begin=0.1        initial value of epsilon gradient         ",
" epos_end=0.5          final value of epsilon gradient          ",
" 								",
" delta_min=            the minimum delta bound by geology          ",
" delta_max=            the maximum delta bound by geology          ",
" 								",
" delta_begin=0.1       initial value of delta gradient         ",
" delta_end=0.5         final value of delta gradient          ",
"                                                                 ",
" fx=0                  first lateral sample in velocity          ",
" fz=0                  first depth sample in velocity            ",
" dx=                   the grid space for velocity model         ",
" dz=                   the grid space for velocity model         ",
" nx=                   number of lateral samples in velocity     ",
" nz=                   number of depth samples in velocity       ",
" DX=                   the B-spline space for velocity model     ",
" DZ=                   the B-spline space for velocity model     ",
" 								",
" nmax_iter=            number of iteration          ",
" ndip_iter=            number of dip iteration          ",
" 								",
" regu_ray_damp=        damping regularization for ray update      ",
" regu_damp=            damping regularization for model update         ",
" regu_homox=           smoothing regularization for velocity update          ",
" regu_homoz=           smoothing regularization for velocity update   ",
" verbose=              note(value=1 support print)          ",
" 								",
" 								",
" Notes:                                                                ",
"	stereotomography     					",
" 				wyx 1.1_30803			",
NULL };
/**************** end self doc ***********************************/
int main(int argc, char **argv)
{
	/* multithreading parameter */
	int myid, numprocs;
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int  ret;
	ret = MPI_Init(&argc, &argv);
	if (ret) {
		printf("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, ret);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	fprintf(stdout, "Process %d of %d is on %s\n", myid, numprocs, processor_name);
	fflush(stdout);
	initargs(argc, argv);
	requestdoc(0); /* stdin not used */

	/*Initialize*/
	int i;
	int iter;//iteration
	int nPrevious;//num of raypairs in perv iteration
	double tbegin;//starting time of some function
	int bool_prior_vel;
	int bool_prior_epos;
	int bool_prior_delta;
	float **prior_vel;
	float **prior_epos;
	float **prior_delta;//prior background field 
	float *pf,*pn;
	float **vcount;
	float **vcount_ini;
	float **regu_vij;
	float **Hess;
	float *error,*misfit;//residual
	int ndata;//num of raypairs
	int Nmodel;//size of Frechet matrix
	float nMemoryMbytes;//Memory cost
	char vij_filename[500];
	char eij_filename[500];
	char dij_filename[500];
	char velocity_filename[500];
	char epos_filename[500];
	char delta_filename[500];
	char data_space_filename[500];
	char model_space_filename[500];
	char filter_data_filename[500];
	char filter_model_filename[500];
	char bspline_filename[500];
	char cost_fun_filename[500];
	char dip_bar_filename[500];
	char gradient_filename[500];
	char data_misfit_filename[500];
	Geo2d geo2dv;//observing system
	bspline_para bpx, bpz;//bspline node
	basisfun_para basis;//bspline function
	DataSpace dp_true, dp_cal;//dataspace
	ModelSpace mp0, mp1;//modelspace
	BGfield **vel, **epos, **delta;//Background field

	/* input information */
	char *parfile;
	char *modelfile;
	char *datadir;//file path
	char *prior_vel_filename;
	char *prior_epos_filename;
	char *prior_delta_filename;
	int nx, nz;
	int datatype;
	int nmax_iter, ndip_iter;
	int verbose;
	float dt;
	float dstep;
	float fx, fz;
	float dx, dz;
	float DX, DZ;
	float max_dip;
	float regu_ray_damp, regu_damp, regu_homox, regu_homoz;//regu coe
	float vel_begin, vel_end;
	float vel_max, vel_min;
	float epos_begin, epos_end;
	float epos_max, epos_min;//wx
	float delta_begin, delta_end;
	float delta_max, delta_min;

	MUSTGETPARSTRING("invfile", &parfile);
	MUSTGETPARSTRING("modelfile", &modelfile);
	MUSTGETPARSTRING("datadir", &datadir);
	if (!getparint("datatype", &datatype)) datatype = 0; /* 0 for text and 1 for binary*/
	if (getparstring("velocity", &prior_vel_filename))
		bool_prior_vel = 1;
	else
		bool_prior_vel = 0;
	if (getparstring("anisotropic_para_epos", &prior_epos_filename))
		bool_prior_epos = 1;
	else
		bool_prior_epos = 0;
	if (getparstring("anisotropic_para_delta", &prior_delta_filename))
		bool_prior_delta = 1;
	else
		bool_prior_delta = 0;
	if (!getparfloat("vel_begin", &vel_begin))   vel_begin = 1000;
	if (!getparfloat("vel_end", &vel_end))	    vel_end = 4000;
	if (!getparfloat("vel_max", &vel_max))	    vel_max = 6500;
	if (!getparfloat("vel_min", &vel_min))	    vel_min = 500;
	if (!getparfloat("epos_begin", &epos_begin)) epos_begin = 0.;
	if (!getparfloat("epos_end", &epos_end))	    epos_end = 0.5;
	if (!getparfloat("epos_max", &epos_max))	    epos_max = 1.0;
	if (!getparfloat("epos_min", &epos_min))	    epos_min = -1.0;//wx
	if (!getparfloat("delta_begin", &delta_begin)) delta_begin = 0.;
	if (!getparfloat("delta_end", &delta_end))	    delta_end = 0.5;
	if (!getparfloat("delta_max", &delta_max))	    delta_max = 1.0;
	if (!getparfloat("delta_min", &delta_min))	    delta_min = 0.0;//wx
	if (!getparint("nmax_iter", &nmax_iter))	    nmax_iter = 20;  
	if (!getparint("ndip_iter", &ndip_iter))	    ndip_iter = 3;
	if (!getparfloat("regu_ray_damp", &regu_ray_damp))	regu_ray_damp = 0;
	if (!getparfloat("regu_damp", &regu_damp))	 regu_damp = 0;
	if (!getparfloat("regu_homox", &regu_homox)) regu_homox = 0;
	if (!getparfloat("regu_homoz", &regu_homoz)) regu_homoz = 60;
	if (!getparfloat("max_dip", &max_dip))	max_dip = 60;
	if (!getparfloat("dt", &dt))	dt = 0.001;
	if (!getparfloat("dstep", &dstep)) dstep = 1;
	if (!getparfloat("fx", &fx))	fx = 0;
	if (!getparfloat("fz", &fz))	fz = 0;
	if (!getparfloat("dx", &dx))	dx = 1;
	if (!getparfloat("dz", &dz))	dz = 1;
	if (!getparint("nx", &nx))	    nx = 1;
	if (!getparint("nz", &nz))	    nz = 1;
	if (!getparfloat("DX", &DX))	DX = 1;
	if (!getparfloat("DZ", &DZ))	DZ = 1;
	if (!getparint("verbose", &verbose))	verbose = 0;
	checkpars();

	nMemoryMbytes = 0;
	/*define spline function*/
	bparameter_initial(&bpx, DX, dx, nx);
	bparameter_initial(&bpz, DZ, dz, nz);
	bpx.node = ealloc1float(bpx.Nnode);
	bpz.node = ealloc1float(bpz.Nnode);
	nMemoryMbytes += bpx.Nnode * FSIZE;
	nMemoryMbytes += bpz.Nnode * FSIZE;
	node_initial(&bpx);
	node_initial(&bpz);

	basis.bsplinex = ealloc2float(bpx.Nfun, bpx.Nnode);
	basis.dbsplinex = ealloc2float(bpx.Nfun, bpx.Nnode);
	basis.ddbsplinex = ealloc2float(bpx.Nfun, bpx.Nnode);
	basis.bsplinez = ealloc2float(bpz.Nfun, bpz.Nnode);
	basis.dbsplinez = ealloc2float(bpz.Nfun, bpz.Nnode);
	basis.ddbsplinez = ealloc2float(bpz.Nfun, bpz.Nnode);
	nMemoryMbytes += bpx.Nfun * bpx.Nnode * FSIZE;
	nMemoryMbytes += bpx.Nfun * bpx.Nnode * FSIZE;
	nMemoryMbytes += bpx.Nfun * bpx.Nnode * FSIZE;
	nMemoryMbytes += bpz.Nfun * bpz.Nnode * FSIZE;
	nMemoryMbytes += bpz.Nfun * bpz.Nnode * FSIZE;
	nMemoryMbytes += bpz.Nfun * bpz.Nnode * FSIZE;
	basis_spline(bpx, basis.bsplinex, basis.dbsplinex, basis.ddbsplinex); /* calculation */
	basis_spline(bpz, basis.bsplinez, basis.dbsplinez, basis.ddbsplinez);

	float **dip_sum, **dip_count;
	dip_sum = ealloc2float(bpz.Nfun, bpx.Nfun);
	dip_count = ealloc2float(bpz.Nfun, bpx.Nfun);
	nMemoryMbytes += bpz.Nfun * bpx.Nfun *FSIZE;
	nMemoryMbytes += bpz.Nfun * bpx.Nfun *FSIZE;
	
	if (verbose && myid == 0) {
		sprintf(bspline_filename, "%s/basis.bsplinex", datadir);
		write2float(basis.bsplinex, bpx.Nfun, bpx.Nnode, bspline_filename);
		sprintf(bspline_filename, "%s/basis.dbsplinex", datadir);
		write2float(basis.dbsplinex, bpx.Nfun, bpx.Nnode, bspline_filename);
		sprintf(bspline_filename, "%s/basis.bsplinez", datadir);
		write2float(basis.bsplinez, bpz.Nfun, bpz.Nnode, bspline_filename);
		sprintf(bspline_filename, "%s/basis.dbsplinez", datadir);
		write2float(basis.dbsplinez, bpz.Nfun, bpz.Nnode, bspline_filename);
	}
	/*screen output parameter*/
	if (myid == 0) {
		printf("data dir:%s\n", datadir);
		printf("invariants file:%s\n", parfile);
		printf("initail model file:%s\n", modelfile);
		printf("data type:%s\n", datatype ? "binary" : "text");
		printf("velocity allowed min:%.1f max:%.1f\n", vel_min, vel_max);
		if (!bool_prior_vel) printf("velocity begin:%.1f end:%.1f\n", vel_begin, vel_end);
		printf("anisotropic para epos allowed min:%.1f max:%.1f\n", epos_min, epos_max);
		if (!bool_prior_epos) printf("epos begin:%.1f end:%.1f\n", epos_begin, epos_end);
		printf("anisotropic para delta allowed min:%.1f max:%.1f\n", delta_min, delta_max);
		if (!bool_prior_delta) printf("delta begin:%.1f end:%.1f\n", delta_begin, delta_end);
		printf("max dip:%.1f degree\n", max_dip);
		printf("dip iteration:%d max iteration:%d\n", ndip_iter, nmax_iter);
		printf("bspline parameters are:\n");
		printf("bpx.Norigin=%6d, bpz.Norigin=%6d\n", bpx.Norigin, bpz.Norigin);
		printf("bpx.Nfun   =%6d, bpz.Nfun   =%6d\n", bpx.Nfun, bpz.Nfun);
		printf("bpx.Nnode  =%6d, bpz.Nnode  =%6d\n", bpx.Nnode, bpz.Nnode);
		printf("fx=%8.1f\tdx=%5.1f nx=%6d DX=%5.1f xrange=%8.1f xmax=%11.f\n", fx, dx, nx, DX, dx*(nx - 1), fx + dx * (nx - 1));
		printf("fz=%8.1f\tdz=%5.1f nz=%6d DZ=%5.1f zrange=%8.1f zmax=%11.f\n", fz, dz, nz, DZ, dz*(nz - 1), fz + dz * (nz - 1));
	}
	
	/*---------------------------------------------define model space------------------------------------*/
	mp0.vnx = bpx.Norigin;
	mp0.vnz = bpz.Norigin;
	mp0.vij = ealloc2float(bpz.Norigin, bpx.Norigin);
	mp0.eij = ealloc2float(bpz.Norigin, bpx.Norigin);
	mp0.dij = ealloc2float(bpz.Norigin, bpx.Norigin);
	mp1.vnx = bpx.Norigin;
	mp1.vnz = bpz.Norigin;
	mp1.vij = ealloc2float(bpz.Norigin, bpx.Norigin);
	mp1.eij = ealloc2float(bpz.Norigin, bpx.Norigin);
	mp1.dij = ealloc2float(bpz.Norigin, bpx.Norigin);
	regu_vij= ealloc2float(bpz.Norigin, bpx.Norigin);
	nMemoryMbytes += bpz.Norigin * bpx.Norigin * FSIZE;
	nMemoryMbytes += bpz.Norigin * bpx.Norigin * FSIZE;
	nMemoryMbytes += bpz.Norigin * bpx.Norigin * FSIZE;
	nMemoryMbytes += bpz.Norigin * bpx.Norigin * FSIZE;
	nMemoryMbytes += bpz.Norigin * bpx.Norigin * FSIZE;
	nMemoryMbytes += bpz.Norigin * bpx.Norigin * FSIZE;
	nMemoryMbytes += bpz.Norigin * bpx.Norigin * FSIZE;
	/*---------------------------------------------define velocity model------------------------------------*/
	vel = (BGfield**)ealloc2(nz, nx, sizeof(BGfield));
	epos = (BGfield**)ealloc2(nz, nx, sizeof(BGfield));
	delta = (BGfield**)ealloc2(nz, nx, sizeof(BGfield));
	prior_vel = ealloc2float(nz, nx);//wx
	prior_epos = ealloc2float(nz, nx);//wx
	prior_delta = ealloc2float(nz, nx);//wx
	vcount = ealloc2float(nz, nx);
	vcount_ini = ealloc2float(nz, nx);
	nMemoryMbytes += nx * nz * sizeof(BGfield);
	nMemoryMbytes += nx * nz * sizeof(BGfield);
	nMemoryMbytes += nx * nz * sizeof(BGfield);
	nMemoryMbytes += nx * nz *FSIZE;
	nMemoryMbytes += nx * nz *FSIZE;
	nMemoryMbytes += nx * nz *FSIZE;
	nMemoryMbytes += nx * nz *FSIZE;
	/*------------------------------------------define observation system-------------------------------------*/
	geo2dv.dt = dt; 
	geo2dv.nx = nx;	geo2dv.nz = nz;
	geo2dv.dx = dx;	geo2dv.dz = dz;	     
	geo2dv.fx = fx; geo2dv.fz = fz;
	geo2dv.xmin = fx; geo2dv.zmin = fz;
	geo2dv.xmax = fx + (nx - 1)*dx; geo2dv.zmax = fz + (nz - 1)*dz;
	/*-------------------------------------initialize velocity space-------------------------------------------*/
	if (myid == 0) printf("initialize velocity model\n");
	/* initialize velocity model */
	if (bool_prior_vel) {
		if (myid == 0) {
			printf("read prior file:%s\n", prior_vel_filename);
			sprintf(vij_filename, "%s/vij_0.dat", datadir);
			if (fopen(vij_filename, "r") == NULL) {
				read2float(prior_vel, nx, nz, prior_vel_filename);
				bspline_interpolation(prior_vel, bpx, bpz, mp0.vij, basis);
			}
			else
				read2float(mp0.vij, bpx.Norigin, bpz.Norigin, vij_filename);
			cal_vel_field(bpx, bpz, basis, mp0.vij, prior_vel);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&mp0.vij[0][0], bpz.Norigin*bpx.Norigin, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&prior_vel[0][0], bpx.Nfun*bpz.Nfun, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}
	else {
		if (myid == 0) printf("we construct gradient velocity\n");
		initialize_velocity(bpx, bpz, basis, mp0.vij, mp0.vnx, mp0.vnz, prior_vel, vel_begin, vel_end);
	}
	field_assign(prior_vel, vel, bpx.Nfun, bpz.Nfun);
	if (myid == 0) {
		//write velocity 0
		sprintf(velocity_filename, "%s/velocity_0.dat", datadir);
		write2field(vel, bpx.Nfun, bpz.Nfun, velocity_filename);
		//write vij 0
		sprintf(vij_filename, "%s/vij_0.dat", datadir);
		write2float(mp0.vij, bpx.Norigin, bpz.Norigin, vij_filename);
		printf("initialize velocity model successful\n");
	}
	/*-------------------------------------initialize epsilon space-------------------------------------------*/
	if (myid == 0) printf("initialize epos model\n");
	/* initialize epos model */
	if (bool_prior_epos) {
		if (myid == 0) {
			printf("read prior file:%s\n", prior_epos_filename);
			sprintf(eij_filename, "%s/eij_0.dat", datadir);
			if (fopen(eij_filename, "r") == NULL) {
				read2float(prior_epos, nx, nz, prior_epos_filename);
				bspline_interpolation(prior_epos, bpx, bpz, mp0.eij, basis);
			}
			else
				read2float(mp0.eij, bpx.Norigin, bpz.Norigin, eij_filename);
			cal_vel_field(bpx, bpz, basis, mp0.eij, prior_epos);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&mp0.eij[0][0], bpz.Norigin*bpx.Norigin, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&prior_epos[0][0], bpx.Nfun*bpz.Nfun, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}
	else {
		if (myid == 0) printf("we construct gradient epos\n");
		initialize_velocity(bpx, bpz, basis, mp0.eij, mp0.vnx, mp0.vnz, prior_epos, epos_begin, epos_end);
	}//wx
	field_assign(prior_epos, epos, bpx.Nfun, bpz.Nfun);
	if (myid == 0) {
		//write epos 0
		sprintf(epos_filename, "%s/epos_0.dat", datadir);
		write2field(epos, bpx.Nfun, bpz.Nfun, epos_filename);
		//write eij 0
		sprintf(eij_filename, "%s/eij_0.dat", datadir);
		write2float(mp0.eij, bpx.Norigin, bpz.Norigin, eij_filename);
		printf("initialize epos model successful\n");
	}//wx
	/*-------------------------------------initialize delta space-------------------------------------------*/
	if (myid == 0) printf("initialize delta model\n");
	/* initialize delta model */
	if (bool_prior_delta) {
		if (myid == 0) {
			printf("read prior file:%s\n", prior_delta_filename);
			sprintf(dij_filename, "%s/dij_0.dat", datadir);
			if (fopen(dij_filename, "r") == NULL) {
				read2float(prior_delta, nx, nz, prior_delta_filename);
				bspline_interpolation(prior_delta, bpx, bpz, mp0.dij, basis);
			}
			else
				read2float(mp0.dij, bpx.Norigin, bpz.Norigin, dij_filename);
			cal_vel_field(bpx, bpz, basis, mp0.dij, prior_delta);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&mp0.dij[0][0], bpx.Norigin*bpz.Norigin, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&prior_delta[0][0], bpx.Nfun*bpz.Nfun, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}
	else {
		if (myid == 0) printf("we construct gradient delta\n");
		initialize_velocity(bpx, bpz, basis, mp0.dij, mp0.vnx, mp0.vnz, prior_delta, delta_begin, delta_end);
	}//wx
	field_assign(prior_delta, delta, bpx.Nfun, bpz.Nfun);
	if (myid == 0) {
		//write delta 0
		sprintf(delta_filename, "%s/delta_0.dat", datadir);
		write2field(delta, bpx.Nfun, bpz.Nfun, delta_filename);
		//write nij 0
		sprintf(dij_filename, "%s/dij_0.dat", datadir);
		write2float(mp0.dij, bpx.Norigin, bpz.Norigin, dij_filename);
		printf("initialize delta model successful\n");
	}//wx
	/* compute second derivatives of background field */
	dv2(geo2dv, vel);
	dv2(geo2dv, epos);//wx
	dv2(geo2dv, delta);//wx

	sprintf(filter_data_filename, "%s/filter_data.txt", datadir);
	sprintf(filter_model_filename, "%s/filter_model.txt", datadir);
	if (myid == 0){
		check_data(parfile, filter_data_filename, modelfile, filter_model_filename, datatype, &dp_true, geo2dv, vel, epos, delta);
	//	check_data_ini(parfile, filter_data_filename, datatype, &dp_true, geo2dv, vel, epos, delta);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&dp_true.ndata, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ndata = dp_true.ndata;
	dp_cal.ndata = ndata;
	if (ndata < 5 && myid == 0) {
		printf("ndata=%d, program will exit.\n", ndata);
		exit(1);
	}
	/* define dataspace and modelspace */
	dp_true.d = (Data_para*)ealloc1(ndata, sizeof(Data_para));
	dp_cal.d = (Data_para*)ealloc1(ndata, sizeof(Data_para));
	mp0.m = (Model_para*)ealloc1(ndata, sizeof(Model_para));
	mp1.m = (Model_para*)ealloc1(ndata, sizeof(Model_para));
	nMemoryMbytes += ndata * sizeof(Data_para);
	nMemoryMbytes += ndata * sizeof(Data_para);
	nMemoryMbytes += ndata * sizeof(Model_para);
	nMemoryMbytes += ndata * sizeof(Model_para);

	/* input data */
	input_data(filter_data_filename, datatype, &dp_true);
	if (myid == 0) {
		printf("first 5 input data:\n");
		for (i = 0; i < 5; i++) {
			printf("sx=%.1f spx=%0.9f rx=%.1f rpx=%0.9f t=%0.7f\n", dp_true.d[i].sx, dp_true.d[i].spx, dp_true.d[i].rx, dp_true.d[i].rpx, dp_true.d[i].t);
		}
		printf("input data sucessfully\n");
		printf("ndata=%d\n", dp_true.ndata);
	}
	mp0.ndata = ndata;
	mp1.ndata = ndata;
	input_model(filter_model_filename, datatype, &mp0);

	/* define misfit */
	Nmodel = 3 * bpx.Norigin * bpz.Norigin;
	if (myid == 0)
		printf("Nmodel=%d\n", Nmodel);

	misfit = ealloc1float(ndata);
	pf = ealloc1float(Nmodel);
	pn = ealloc1float(Nmodel);
	Hess = ealloc2float(Nmodel, Nmodel);
	error = ealloc1float(nmax_iter);

	zero1float(misfit, ndata);
	zero1float(pf, Nmodel);
	zero1float(pn, Nmodel);
	zero1float(error, nmax_iter);
	zero2float(Hess, Nmodel, Nmodel);

	nMemoryMbytes += ndata * FSIZE;
	nMemoryMbytes += Nmodel * FSIZE;
	nMemoryMbytes += Nmodel * FSIZE;
	nMemoryMbytes += Nmodel * FSIZE;
	nMemoryMbytes += nmax_iter * FSIZE;
	nMemoryMbytes += Nmodel * Nmodel * FSIZE;

	if (myid == 0) printf("Memory cost %.1f G bytes \n", nMemoryMbytes / 1024 / 1024 / 1024);

	/*------------------------initialize model---------------------------*/
//	initialize_model(&dp_true, &mp0, geo2dv, vel, epos, delta, myid, numprocs);
//	if (myid == 0) printf("-------\ninitialize model sucessfully\n");
	
    //write value of cost function with iteration
	forward_modeling(&mp0, &dp_cal, geo2dv, vel, epos, delta, vcount_ini, myid, numprocs);
	cal_cost_fun(&error[0], misfit, &dp_true, &dp_cal);
	sprintf(cost_fun_filename, "%s/cost_func.dat", datadir);
	write1float(&error[0], 1, cost_fun_filename);
	/* write initial dip_bar */
	sprintf(dip_bar_filename, "%s/dip_bar_0.txt", datadir);
	dip_bar(&mp0, dip_bar_filename, myid, numprocs);
	/*-------------------------------start model update to better fit data---------------------------*/
	for (iter = 0; iter < nmax_iter; iter++) {
		if (myid == 0) printf("\niter=%d\n", iter + 1);
		/*calculate surface data and Frechet derivative*/
		tbegin = tic();
		cal_Frechet_matrix(pf, Hess, bpx, bpz, basis, &mp0, &dp_true, iter, ndip_iter, Nmodel, geo2dv, vel, epos, delta, myid, numprocs);
		if (myid == 0) printf("cal_steepest_gradient time: %8.2lf s\n", toc(tbegin));

		if (myid == 0) {
			if(iter >= ndip_iter){
				/*dip statistic*/
				dip_statistic_small(&mp0, geo2dv, bpx, bpz, dip_sum, dip_count);
				printf("none zero dip=%d\n", none_zero_count(dip_sum, bpx.Nfun, bpz.Nfun));

				sprintf(gradient_filename, "%s/dip_sum_%d.dat", datadir, iter + 1);
				if (verbose) write2float(dip_sum, bpx.Nfun, bpz.Nfun, gradient_filename);
				sprintf(gradient_filename, "%s/dip_count_%d.dat", datadir, iter + 1);
				if (verbose) write2float(dip_count, bpx.Nfun, bpz.Nfun, gradient_filename);

				/*calculate homo regularization equations*/
				printf("homo_regu\n");
				cal_homo_regu(bpx, bpz, regu_homox, regu_homoz, Hess, ndata);
				vij_statistic(bpx, bpz, regu_vij, vcount);
				sprintf(gradient_filename, "%s/regu_vij_%d.dat", datadir, iter + 1);
				if (verbose) write2float(regu_vij, bpx.Norigin, bpz.Norigin, gradient_filename);
			}
			printf("damp_regu\n");
			cal_damp_regu(Hess, ndata, bpx, bpz, regu_ray_damp, regu_damp, iter, ndip_iter, regu_vij);
			printf("cal newton gradient\n");
			tbegin = tic();
			cal_Newton_gradient(pf, Hess, pn, Nmodel);
			printf("cal_newton_gradient time: %8.2lf s\n", toc(tbegin));
			if (verbose == 1) {
				//write F Matrix
				sprintf(gradient_filename, "%s/sd_gradient_%d.dat", datadir, iter + 1);
				printf("write sd_gradient:%s\n", gradient_filename);
				write1float(pf, Nmodel, gradient_filename);
				sprintf(gradient_filename, "%s/nt_gradient_%d.dat", datadir, iter + 1);
				printf("write nt_gradient:%s\n", gradient_filename);
				write1float(pn, Nmodel, gradient_filename);
				sprintf(gradient_filename, "%s/Hess_%d.dat", datadir, iter + 1);
				printf("write Hess:%s\n", gradient_filename);
				write1float(&Hess[0][0], Nmodel*Nmodel, gradient_filename);
			}
			//solve tomographic equations
			if (iter < ndip_iter)
				printf("---------\niterations updating ray segment only!\n");
			else
				printf("---------\niterations updating field and ray segments!\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(pn, Nmodel, MPI_FLOAT, 0, MPI_COMM_WORLD);
		tbegin = tic();
		/*update model space while cost function decrease*/
		update_model(&mp0, &mp1, &dp_cal, &dp_true, bpx, bpz, basis, pn, &dstep, misfit, error, iter, ndip_iter, geo2dv, vel, epos, delta, vel_min, vel_max, epos_min, epos_max, delta_min, delta_max, vcount, myid, numprocs);
		if (myid == 0) printf("update_model time: %8.2lf s\n", toc(tbegin));

		if (myid == 0) {
			sprintf(velocity_filename, "%s/velocity_%d.dat", datadir, iter + 1);
			write2field(vel, bpx.Nfun, bpz.Nfun, velocity_filename);
			sprintf(epos_filename, "%s/epos_%d.dat", datadir, iter + 1);
			write2field(epos, bpx.Nfun, bpz.Nfun, epos_filename);
			sprintf(delta_filename, "%s/delta_%d.dat", datadir, iter + 1);
			write2field(delta, bpx.Nfun, bpz.Nfun, delta_filename);
			//write vij eij dij
			sprintf(vij_filename, "%s/vij_%d.dat", datadir, iter+1);
			write2float(mp0.vij, bpx.Norigin, bpz.Norigin, vij_filename);
			sprintf(eij_filename, "%s/eij_%d.dat", datadir, iter+1);
			write2float(mp0.eij, bpx.Norigin, bpz.Norigin, eij_filename);
			sprintf(dij_filename, "%s/dij_%d.dat", datadir, iter+1);
			write2float(mp0.dij, bpx.Norigin, bpz.Norigin, dij_filename);

			//write dip bar
			sprintf(dip_bar_filename, "%s/dip_bar_%d.txt", datadir, iter + 1);
			dip_bar(&mp0, dip_bar_filename, myid, numprocs);
		}
		
		if (iter >= ndip_iter - 2) {
			nPrevious = ndata;
			/*report the error and delete the wrong data*/
		//	error_report(&mp0, &mp1, &dp_true, &dp_cal, &ndata, &Nmodel, geo2dv,bpx,bpz, misfit,max_dip,verbose, myid,numprocs);
			if (myid == 0) {
				if (nPrevious - ndata) printf("new ndata=%d, del=%d\n", ndata, nPrevious - ndata);
				//write filterd file 
				sprintf(filter_data_filename, "%s/filtered_data_%d.txt", datadir, iter + 1);
				write_data(filter_data_filename, datatype, &dp_true);
			}
		}
		
		if (myid == 0) {
			//write value of cost function with iteration
			write1float_append(&error[iter], 1, cost_fun_filename);
			//write data space
			sprintf(data_space_filename, "%s/data_space_%d.dat", datadir, iter + 1);
			write1(&dp_true.d[0], sizeof(Data_para), ndata, data_space_filename);
			//write data space
			sprintf(model_space_filename, "%s/model_space_%d.dat", datadir, iter + 1);
			write1(&mp0.m[0], sizeof(Model_para), ndata, model_space_filename);
		}

		//zero array mF and data_misfit
		zero1float(misfit, ndata);
		zero1float(pf, Nmodel);
		zero1float(pn, Nmodel);
		zero2float(Hess, Nmodel, Nmodel);
		MPI_Barrier(MPI_COMM_WORLD);
	}//end of loop iterations
	if (myid == 0) printf("well done\n");
	MPI_Finalize();

	free1float(bpx.node);
	free1float(bpz.node);
	free1float(pf);
	free1float(pn);
	free2float(Hess);
	free2float(basis.bsplinex);
	free2float(basis.bsplinez);
	free2float(basis.dbsplinex);
	free2float(basis.dbsplinez);
	free2float(basis.ddbsplinex);
	free2float(basis.ddbsplinez);
	free2float(dip_sum);
	free2float(dip_count);
	free2float(mp0.vij);
	free2float(mp1.vij);
	free2float(mp0.eij);
	free2float(mp1.eij);
	free2float(mp0.dij);
	free2float(mp1.dij);
	free2float(regu_vij);
	free2float(vcount);
	free2float(prior_vel);
	free2float(prior_epos);
	free2float(prior_delta);
	free1(dp_true.d);
	free1(dp_cal.d);
	free1(mp0.m);
	free1(mp1.m);
	free1float(misfit);
	free1float(error);
	free2((void**)vel);
	free2((void**)epos);
	free2((void**)delta);
	
	return	EXIT_SUCCESS;
}
