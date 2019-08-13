#include "mpi.h"
#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include "RaySub.h"
#include "bsplines.h"
#include "forward_modeling.h"

void DRayShot(ModelSpace *mp, Data_para *dp_cal, Ray *rayA, Ray *rayB, int i, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	float thetaA, thetaB;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	double VphA, VphB;
	double si2A2, siA2;
	double si2B2, siB2;
	Source source;

	source.x = mp->m[i].xcor;
	source.z = mp->m[i].zcor;
	vel2Interp(geo2dv, vel, source.x, source.z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
	vel2Interp(geo2dv, epos, source.x, source.z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
	vel2Interp(geo2dv, delta, source.x, source.z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);

	thetaA = angle2radian(180 - mp->m[i].theta_s);
	siA2 = sin(thetaA)*sin(thetaA);
	si2A2 = sin(2 * thetaA)*sin(2 * thetaA);
	VphA = vel0 * sqrt(0.5 + epos0 * siA2 + 0.5 * sqrt(pow(1 + 2 * epos0 * siA2, 2) - 2 * (epos0 - delta0) * si2A2));
	source.px = sin(thetaA) / VphA;
	source.pz = cos(thetaA) / VphA;
	dynamic_ray_tracing(rayA, source, geo2dv, vel, epos, delta, 99);

	thetaB = angle2radian(180 - mp->m[i].theta_r);
	siB2 = sin(thetaB)*sin(thetaB);
	si2B2 = sin(2 * thetaB)*sin(2 * thetaB);
	VphB = vel0 * sqrt(0.5 + epos0 * siB2 + 0.5 * sqrt(pow(1 + 2 * epos0 * siB2, 2) - 2 * (epos0 - delta0) * si2B2));
	source.px = sin(thetaB) / VphB;
	source.pz = cos(thetaB) / VphB;
	dynamic_ray_tracing(rayB, source, geo2dv, vel, epos, delta, 99);

	dp_cal->sx = rayA->rs[rayA->nrs - 1].x;
	dp_cal->spx = rayA->rs[rayA->nrs - 1].px;
	dp_cal->rx = rayB->rs[rayB->nrs - 1].x;
	dp_cal->rpx = rayB->rs[rayB->nrs - 1].px;
	dp_cal->t = rayA->rs[rayA->nrs - 1].t + rayB->rs[rayB->nrs - 1].t;
	return;
}

void forward_modeling(ModelSpace *mp1, DataSpace *dp_cal, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float **vcount_float, int myid, int numprocs)
{
	int i, j;
	int ndata;
	int *recvcounts;
	int *displs;

	int **vcount;
	int **vcount_int;
	vcount = ealloc2int(geo2dv.nz, geo2dv.nx);
	vcount_int = ealloc2int(geo2dv.nz, geo2dv.nx);
	zero2int(vcount, geo2dv.nz, geo2dv.nx);
	zero2int(vcount_int, geo2dv.nz, geo2dv.nx);
	
	Ray rayA, rayB;
	ndata = dp_cal->ndata;
	rayA.rs = (RayStep*)ealloc1(MAX_RAYSTEP, sizeof(RayStep));
	rayB.rs = (RayStep*)ealloc1(MAX_RAYSTEP, sizeof(RayStep));
	recvcounts = ealloc1int(numprocs);
	displs = ealloc1int(numprocs);
	mpi_gatherv_construct(0, ndata, numprocs, 1, recvcounts, displs);
	
	for (i = displs[myid], j = 0; j < recvcounts[myid]; i++, j++) {
		KRayShot(mp1, dp_cal, &rayA, &rayB, i, j, geo2dv, vel, epos, delta);
		RayPathAnalysis(&rayA, geo2dv, vcount);
		RayPathAnalysis(&rayB, geo2dv, vcount);
	}
	
	/* begin of type define */
	MPI_Datatype d_stype, d_type[1];
	MPI_Aint d_disp[1];
	int d_blocklen[1];
	d_type[0] = MPI_FLOAT;
	d_blocklen[0] = 5; ///sizeof(Data_para)/sizeof(float)
	d_disp[0] = 0;
	MPI_Type_struct(1, d_blocklen, d_disp, d_type, &d_stype);
	MPI_Type_commit(&d_stype);
    /* end of type define */
	Data_para *dp_tmp;
	if (myid == 0)
		dp_tmp = (Data_para*)ealloc1(ndata, sizeof(Data_para));
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv((void *)&dp_cal->d[0], recvcounts[myid], d_stype, (void *)&dp_tmp[0], recvcounts, displs, d_stype, 0, MPI_COMM_WORLD);

	if (myid == 0) {
		for (i = 0; i < ndata; i++)
			dp_cal->d[i] = dp_tmp[i];
	}
	MPI_Bcast(&dp_cal->d[0], ndata, d_stype, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&vcount[0][0], &vcount_int[0][0], geo2dv.nz*geo2dv.nx, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(myid==0){
		for(i = 0; i < geo2dv.nx; i++)
		  for(j = 0; j < geo2dv.nz; j++)
			vcount_float[i][j] = vcount_int[i][j];
	}

	free2int(vcount_int);
	free2int(vcount);
	free1(rayA.rs);
	free1(rayB.rs);
	free1int(recvcounts);
	free1int(displs);
	if (myid == 0)
		free1(dp_tmp);
	return;
}

void KRayShot(ModelSpace *mp, DataSpace *dp_cal, Ray *rayA, Ray *rayB, int i, int j, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	float thetaA, thetaB;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	double VphA, VphB;
	double si2A2, siA2;
	double si2B2, siB2;
	Source source;

	source.x = mp->m[i].xcor;
	source.z = mp->m[i].zcor;
	vel2Interp(geo2dv, vel, source.x, source.z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
	vel2Interp(geo2dv, epos, source.x, source.z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
	vel2Interp(geo2dv, delta, source.x, source.z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);

	thetaA = angle2radian(180 - mp->m[i].theta_s);
	siA2 = sin(thetaA)*sin(thetaA);
	si2A2 = sin(2 * thetaA)*sin(2 * thetaA);
	VphA = vel0 * sqrt(0.5 + epos0 * siA2 + 0.5 * sqrt(pow(1 + 2 * epos0 * siA2, 2) - 2 * (epos0 - delta0) * si2A2));
	source.px = sin(thetaA) / VphA;
	source.pz = cos(thetaA) / VphA;
	kinematic_ray_tracing(rayA, source, geo2dv, vel, epos, delta, 99);

	thetaB = angle2radian(180 - mp->m[i].theta_r);
	siB2 = sin(thetaB)*sin(thetaB);
	si2B2 = sin(2 * thetaB)*sin(2 * thetaB);
	VphB = vel0 * sqrt(0.5 + epos0 * siB2 + 0.5 * sqrt(pow(1 + 2 * epos0 * siB2, 2) - 2 * (epos0 - delta0) * si2B2));
	source.px = sin(thetaB) / VphB;
	source.pz = cos(thetaB) / VphB;
	kinematic_ray_tracing(rayB, source, geo2dv, vel, epos, delta, 99);
	
	dp_cal->d[j].sx = rayA->rs[rayA->nrs - 1].x;
	dp_cal->d[j].spx = rayA->rs[rayA->nrs - 1].px;
	dp_cal->d[j].rx = rayB->rs[rayB->nrs - 1].x;
	dp_cal->d[j].rpx = rayB->rs[rayB->nrs - 1].px;
	dp_cal->d[j].t = rayA->rs[rayA->nrs - 1].t + rayB->rs[rayB->nrs - 1].t;
	return;
}

void initialize_velocity(bspline_para bpx, bspline_para bpz, basisfun_para basis, float **vij, int vnx, int vnz, float **vel_field, float vel_begin, float vel_end)
{
	int i,j;
	/*-----initialize bspline coe vij----------------*/
	for (i = 0; i < vnx; i++) 
		for (j = 0; j < vnz; j++) {
			vij[i][j] = vel_begin + (vel_end - vel_begin) / (bpz.Norigin) * (j+ 1);
		}
	/*------cal initial vel_field---------*/
	cal_vel_field(bpx, bpz, basis, vij, vel_field);
}

void initialize_model(DataSpace *dp_true, ModelSpace *mp0, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, int myid, int numprocs)
{
	int i, j;
	int ndata;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	float boundary_eps = 1;
	float vph;
	float value_asin;
	Ray rayA, rayB;
	Source source;
		
	ndata = mp0->ndata;
	rayA.rs = (RayStep*)ealloc1(MAX_RAYSTEP, sizeof(RayStep));
	rayB.rs = (RayStep*)ealloc1(MAX_RAYSTEP, sizeof(RayStep));

	int *recvcounts;
	int *displs;
	recvcounts = ealloc1int(numprocs);
	displs = ealloc1int(numprocs);
	mpi_gatherv_construct(0, ndata, numprocs, 1, recvcounts, displs);
	/*--------initialize ray segment with ray tracing--------*/
	for (i = displs[myid], j = 0; j < recvcounts[myid]; i++, j++) {
		//source 
		source.x = dp_true->d[i].sx;
		source.z = 0;
		vel2Interp(geo2dv, vel, source.x, source.z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
		vel2Interp(geo2dv, epos, source.x, source.z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
		vel2Interp(geo2dv, delta, source.x, source.z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
		source.px = -dp_true->d[i].spx;
		source.pz = sqrt((1/vel0/vel0-(1+2*epos0)*source.px*source.px)/(1-2*(epos0-delta0)*vel0*vel0*source.px*source.px));
		kinematic_ray_tracing(&rayA, source, geo2dv, vel, epos, delta, dp_true->d[i].t/2.0);
		vph = 1.0 / sqrt(pow(rayA.rs[rayA.nrs-1].px, 2) + pow(rayA.rs[rayA.nrs-1].pz, 2));
		value_asin = -rayA.rs[rayA.nrs-1].px * vph;
		if (value_asin > 1 || value_asin < -1) {
			printf("id:%d value_asin=%f out of range [-1,1]\n",i,value_asin);
			value_asin = value_asin > 0 ? 1 : -1;
		}
		mp0->m[j].theta_s = radian2angle(asin(value_asin));
		//receiver
		source.x = dp_true->d[i].rx;
		source.z = 0;
		vel2Interp(geo2dv, vel, source.x, source.z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
		vel2Interp(geo2dv, epos, source.x, source.z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
		vel2Interp(geo2dv, delta, source.x, source.z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
		source.px = -dp_true->d[i].rpx;
		source.pz = sqrt((1/vel0/vel0-(1+2*epos0)*source.px*source.px)/(1-2*(epos0-delta0)*vel0*vel0*source.px*source.px));
		kinematic_ray_tracing(&rayB, source, geo2dv, vel, epos, delta, dp_true->d[i].t/2.0);
		vph = 1.0 / sqrt(pow(rayB.rs[rayB.nrs-1].px, 2) + pow(rayB.rs[rayB.nrs-1].pz, 2));
		value_asin = -rayB.rs[rayB.nrs-1].px * vph;
		if (value_asin > 1 || value_asin < -1) {
			printf("id:%d value_asin=%f out of range [-1,1]\n",i,value_asin);
			value_asin=value_asin>0?1:-1;
		}
		mp0->m[j].theta_r = radian2angle(asin(value_asin));
		mp0->m[j].xcor = (rayA.rs[rayA.nrs-1].x + rayB.rs[rayB.nrs-1].x) / 2;
		mp0->m[j].zcor = (rayA.rs[rayA.nrs-1].z + rayB.rs[rayB.nrs-1].z) / 2;
		
		if (mp0->m[j].xcor<geo2dv.xmin + boundary_eps || mp0->m[j].xcor>geo2dv.xmax - boundary_eps || mp0->m[j].zcor > geo2dv.zmax - boundary_eps) {
		printf("ray:%d(model x=%.1f,z=%.1f,theta=%.3f) (data sx=%.1f, ps=%f, rx=%.1f, pr=%f, t=%.3f)pass through the model boundary. just simple warning.\n", i, mp0->m[j].xcor, mp0->m[j].zcor, mp0->m[j].theta_r, dp_true->d[i].sx, dp_true->d[i].spx, dp_true->d[i].rx, dp_true->d[i].rpx, dp_true->d[i].t);
		}
	}
	/* begin of type define */
	MPI_Datatype m_stype, m_type[1];
	MPI_Aint m_disp[1];
	int m_blocklen[1];
	m_type[0] = MPI_FLOAT;
	m_blocklen[0] = 4; ///sizeof(Model_para)/sizeof(float)
	m_disp[0] = 0;
	MPI_Type_struct(1, m_blocklen, m_disp, m_type, &m_stype);
	MPI_Type_commit(&m_stype);
	/* end of type define */
	Model_para *mp_tmp;
	if (myid == 0)
	  mp_tmp = (Model_para*)ealloc1(ndata, sizeof(Model_para));
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv((void *)&mp0->m[0], recvcounts[myid], m_stype, (void *)&mp_tmp[0], recvcounts, displs, m_stype, 0, MPI_COMM_WORLD);
	if (myid == 0) {
		for (i = 0; i < ndata; i++)
		  mp0->m[i] = mp_tmp[i];
	}
	MPI_Bcast(&mp0->m[0], ndata, m_stype, 0, MPI_COMM_WORLD);
	free1int(recvcounts);
	free1int(displs);
	free1(rayA.rs);
	free1(rayB.rs);
	if (myid == 0) free1(mp_tmp);
	return;
}

void RayPathAnalysis(Ray *ray, Geo2d geo2dv, int **vcount)
{
	int i;
	int ix, iz;
	for(i = 0; i < ray->nrs; i++){
		ix = (ray->rs[i].x - geo2dv.fx) / geo2dv.dx + 0.5;
		iz = (ray->rs[i].z - geo2dv.fz) / geo2dv.dz + 0.5;
		if(ix < 0 || ix >= geo2dv.nx || iz < 0 || iz >= geo2dv.nz){
			printf("RayPathAnalysis:\n");
			printf("ix=%d i=%d x=%f fx=%f nrs=%d\n", ix, i, ray->rs[i].x, geo2dv.fx, ray->nrs);
			printf("iz=%d z=%f fz=%f out\n",iz, ray->rs[i].z, geo2dv.fz);
			printf("i=%d x=%f z=%f\n", i - 1, ray->rs[i-1].x, ray->rs[i-1].z);
			printf("i=%d x=%f z=%f\n", i - 2, ray->rs[i-2].x, ray->rs[i-2].z);
			continue;
		}
		vcount[ix][iz]++;
	}
}

