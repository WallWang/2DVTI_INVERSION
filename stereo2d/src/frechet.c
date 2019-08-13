#include "mpi.h"
#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include "RaySub.h"
#include "bsplines.h"
#include "sparse_matrix.h"
#include "frechet.h"
#include "forward_modeling.h"

void cal_Frechet_matrix(float *pf, float **Hess, bspline_para bpx, bspline_para bpz, basisfun_para basis, ModelSpace *mp0, DataSpace *dp_true, int iter, int ndip_iter, int Nmodel, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, int myid, int numprocs)
{
	int i, j, k, m, n;
	int ndata;
	float **Fi;
	float data_misfit[5];
	float *pf_tmp;
	float **Hess_tmp;
	Ray rayA, rayB;
	Data_para dp_cal;
	rayA.rs = (RayStep*)ealloc1(MAX_RAYSTEP, sizeof(RayStep));
	rayB.rs = (RayStep*)ealloc1(MAX_RAYSTEP, sizeof(RayStep));
	Fi = ealloc2float(5, Nmodel);
	pf_tmp = ealloc1float(Nmodel);
	Hess_tmp = ealloc2float(Nmodel, Nmodel);
	ndata = dp_true->ndata;
	zero1float(pf_tmp, Nmodel);
	zero2float(Hess_tmp, Nmodel, Nmodel);

	int *recvcounts;
	int *displs;
	recvcounts = ealloc1int(numprocs);
	displs = ealloc1int(numprocs);
	mpi_gatherv_construct(0, ndata, numprocs, 1, recvcounts, displs);

	for (i = displs[myid], j = 0; j < recvcounts[myid]; i++, j++) {
		zero2float(Fi, 5, Nmodel);
		zero1float(data_misfit, 5);
		DRayShot(mp0, &dp_cal, &rayA, &rayB, i, geo2dv, vel, epos, delta);
		
		/*calculate Frechet derivatives about position and angle*/
//		Frechet_matrix_x_z_sita(dp_cal, Fi, i, &rayA, &rayB, geo2dv, vel, epos, delta);
		
		/*calculate Frechet derivatives about vij eij and nij*/
		if (iter >= ndip_iter) {
			Frechet_matrix_field(Fi, bpx, bpz, basis, mp0, geo2dv, vel, epos, delta, &rayA, &rayB);
		}
		data_fitting_err(&dp_cal, dp_true, data_misfit, i);
		for (m = 0; m < Nmodel; m++){
			for(n = 0; n< Nmodel; n++){
				for(k = 0; k< 5; k++)
				  Hess_tmp[m][n] += Fi[m][k]*Fi[n][k];
			}
			for (k = 0; k < 5; k++)
			  pf_tmp[m] += Fi[m][k] * data_misfit[k];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(pf_tmp, pf, Nmodel, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Hess_tmp[0][0], &Hess[0][0], Nmodel*Nmodel, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	free1int(displs);
	free1int(recvcounts);
	free1float(pf_tmp);
	free2float(Hess_tmp);
	free2float(Fi);
	free1(rayA.rs);
	free1(rayB.rs);
}

void Frechet_matrix_x_z_sita(DataSpace *dp_cal, tsmatrix *F, int index, Ray *rayA, Ray *rayB,Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	int ndata;
	int row;

	float value;
//	float sin_thetas,cos_thetas;
//	float sin_thetar,cos_thetar;
	float pxzA, pxzB;
	float dtA, dtB;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	ndata = dp_cal->ndata;
	
	//extra raytracing raypath a
	vel2Interp(geo2dv, vel, rayA->rs[rayA->nrs - 1].x, rayA->rs[rayA->nrs - 1].z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
	vel2Interp(geo2dv, epos, rayA->rs[rayA->nrs - 1].x, rayA->rs[rayA->nrs - 1].z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
	vel2Interp(geo2dv, delta, rayA->rs[rayA->nrs - 1].x, rayA->rs[rayA->nrs - 1].z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
//	sin_thetas = rayA->rs[0].px*vel0;
//	cos_thetas = -rayA->rs[0].pz*vel0;
	pxzA = (rayA->rs[rayA->nrs-1].px*(1+2*epos0-2*vel0*vel0*(epos0-delta0)*rayA->rs[rayA->nrs-1].pz*rayA->rs[rayA->nrs-1].pz)) / (rayA->rs[rayA->nrs-1].pz*(1-2*vel0*vel0*(epos0-delta0)*rayA->rs[rayA->nrs-1].px*rayA->rs[rayA->nrs-1].px));
	dtA = -1.0/(vel0*vel0*rayA->rs[rayA->nrs-1].pz*(1-2*vel0*vel0*(epos0-delta0)*rayA->rs[rayA->nrs-1].px*rayA->rs[rayA->nrs-1].px));
	//extra raytracing raypath b
	vel2Interp(geo2dv, vel, rayB->rs[rayB->nrs - 1].x, rayB->rs[rayB->nrs - 1].z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
	vel2Interp(geo2dv, epos, rayB->rs[rayB->nrs - 1].x, rayB->rs[rayB->nrs - 1].z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
	vel2Interp(geo2dv, delta, rayB->rs[rayB->nrs - 1].x, rayB->rs[rayB->nrs - 1].z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
//	sin_thetar = rayB->rs[0].px*vel0;
//	cos_thetar = -rayB->rs[0].pz*vel0;
	pxzB = (rayB->rs[rayB->nrs-1].px*(1+2*epos0-2*vel0*vel0*(epos0-delta0)*rayB->rs[rayB->nrs-1].pz*rayB->rs[rayB->nrs-1].pz)) / (rayB->rs[rayB->nrs-1].pz*(1-2*vel0*vel0*(epos0-delta0)*rayB->rs[rayB->nrs-1].px*rayB->rs[rayB->nrs-1].px));
	dtB = -1.0/(vel0*vel0*rayB->rs[rayB->nrs-1].pz*(1-2*vel0*vel0*(epos0-delta0)*rayB->rs[rayB->nrs-1].px*rayB->rs[rayB->nrs-1].px));
	/*----calculate Frechet derivatives about x_z_sita----*/
	row = index * 5;

	value = (rayA->rs[rayA->nrs - 1].Q[1][0] * dtA + rayB->rs[rayB->nrs - 1].Q[1][0] * dtB) / wt;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 0, index + 0 * ndata, value);//t about x
	value = (rayA->rs[rayA->nrs - 1].Q[1][1] * dtA + rayB->rs[rayB->nrs - 1].Q[1][1] * dtB) / wt;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 0, index + 1 * ndata, value);//t about z
	value = -(rayA->rs[0].pz * rayA->rs[rayA->nrs - 1].Q[1][2] - rayA->rs[0].px * rayA->rs[rayA->nrs - 1].Q[1][3]) * dtA / wt;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 0, index + 2 * ndata, value);//t about ts
	value = -(rayB->rs[0].pz * rayB->rs[rayB->nrs - 1].Q[1][2] - rayB->rs[0].px * rayB->rs[rayB->nrs - 1].Q[1][3]) * dtB / wt;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 0, index + 3 * ndata, value);//t about tr

//	value=-(sin_thetas+sin_thetar)/vel0/wt;
//	if(fabsf(value)>1e-7) Matrix_Add(F,row+0,index+0*ndata,value);
//	value=(cos_thetas+cos_thetar)/vel0/wt;
//	if(fabsf(value)>1e-7) Matrix_Add(F,row+0,index+1*ndata,value);


	value = rayA->rs[rayA->nrs - 1].Q[2][0] / wp;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 1, index + 0 * ndata, value);//spx about x
	value = rayA->rs[rayA->nrs - 1].Q[2][1] / wp;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 1, index + 1 * ndata, value);//spx about z
	value = -(rayA->rs[0].pz * rayA->rs[rayA->nrs - 1].Q[2][2] - rayA->rs[0].px * rayA->rs[rayA->nrs - 1].Q[2][3]) / wp;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 1, index + 2 * ndata, value);//spx about theta_s
	
	value = (rayA->rs[rayA->nrs - 1].Q[0][0] - pxzA * rayA->rs[rayA->nrs - 1].Q[1][0]) / wx;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 2, index + 0 * ndata, value);//sx about x
	value = (rayA->rs[rayA->nrs - 1].Q[0][1] - pxzA * rayA->rs[rayA->nrs - 1].Q[1][1]) / wx;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 2, index + 1 * ndata, value);//sx about z
	value = -((rayA->rs[0].pz * rayA->rs[rayA->nrs - 1].Q[0][2] - rayA->rs[0].px * rayA->rs[rayA->nrs - 1].Q[0][3]) - pxzA * (rayA->rs[0].pz * rayA->rs[rayA->nrs - 1].Q[1][2] - rayA->rs[0].px * rayA->rs[rayA->nrs - 1].Q[1][3])) / wx;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 2, index + 2 * ndata, value);//sx about theta_s
	
	value = rayB->rs[rayB->nrs - 1].Q[2][0] / wp;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 3, index + 0 * ndata, value);//rpx about x
	value = rayB->rs[rayB->nrs - 1].Q[2][1] / wp;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 3, index + 1 * ndata, value);//rpx about z
	value = -(rayB->rs[0].pz * rayB->rs[rayB->nrs - 1].Q[2][2] - rayB->rs[0].px * rayB->rs[rayB->nrs - 1].Q[2][3]) / wp;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 3, index + 3 * ndata, value);//rpx about theta_r
	
	value = (rayB->rs[rayB->nrs - 1].Q[0][0] - pxzB * rayB->rs[rayB->nrs - 1].Q[1][0]) / wx;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 4, index + 0 * ndata, value);//rx about x
	value = (rayB->rs[rayB->nrs - 1].Q[0][1] - pxzB * rayB->rs[rayB->nrs - 1].Q[1][1]) / wx;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 4, index + 1 * ndata, value);//rx about z
	value = -((rayB->rs[0].pz * rayB->rs[rayB->nrs - 1].Q[0][2] - rayB->rs[0].px * rayB->rs[rayB->nrs - 1].Q[0][3]) - pxzB * (rayB->rs[0].pz * rayB->rs[rayB->nrs - 1].Q[1][2] - rayB->rs[0].px * rayB->rs[rayB->nrs - 1].Q[1][3])) / wx;
	if (fabsf(value) > 1e-7) Matrix_Add(F, row + 4, index + 3 * ndata, value);//rx about theta_r

	return;
} 

void Frechet_matrix_field(float **F, bspline_para bpx, bspline_para bpz, basisfun_para basis, ModelSpace *mp0, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, Ray *rayA, Ray *rayB)
{
	int i;
	int col;
	int nodex, nodez;
	float stv, sxv, spv, rtv, rxv, rpv;
	int ndata;
	int nbspline;
	DField sf, rf;
	ndata = mp0->ndata;
	nbspline = mp0->vnx * mp0->vnz;

	int **vijA, **vijB;
	vijA = ealloc2int(bpz.Norigin, bpx.Norigin);
	vijB = ealloc2int(bpz.Norigin, bpx.Norigin);
	RayPathZone(geo2dv, bpx, bpz, rayA, vijA);
	RayPathZone(geo2dv, bpx, bpz, rayB, vijB);

	for (nodex = 0; nodex < mp0->vnx; nodex++) 
	  for (nodez = 0; nodez < mp0->vnz; nodez++) {
		  col = mp0->vnz * nodex + nodez;
		  zerofrechet(&sf);
		  zerofrechet(&rf);
		  if (vijA[nodex][nodez]) {
			  RayPathFieldFrechet(rayA, basis, nodex, nodez, geo2dv, vel, epos, delta, &sf);
			  F[col+0*nbspline][1] = sf.pv * vscale;//spx about vij
			  F[col+1*nbspline][1] = sf.pe * escale;//spx about eij
			  F[col+2*nbspline][1] = sf.pd * dscale;//spx about dij
			  F[col+0*nbspline][2] = sf.xv * vscale;//sx about vij
			  F[col+1*nbspline][2] = sf.xe * escale;//sx about eij
			  F[col+2*nbspline][2] = sf.xd * dscale;//sx about dij
		  }
		  if (vijB[nodex][nodez]) {
			  RayPathFieldFrechet(rayB, basis, nodex, nodez, geo2dv, vel, epos, delta, &rf);
			  F[col+0*nbspline][3] = rf.pv * vscale;//rpx about vij
			  F[col+1*nbspline][3] = rf.pe * escale;//rpx about eij
			  F[col+2*nbspline][3] = rf.pd * dscale;//rpx about dij
			  F[col+0*nbspline][4] = rf.xv * vscale;//rx about vij
			  F[col+1*nbspline][4] = rf.xe * escale;//rx about eij
			  F[col+2*nbspline][4] = rf.xd * dscale;//rx about dij
		  }
		  F[col+0*nbspline][0] = (sf.tv + rf.tv) * vscale;//t about vij
		  F[col+1*nbspline][0] = (sf.te + rf.te) * escale;//t about eij
		  F[col+2*nbspline][0] = (sf.td + rf.td) * dscale;//t about dij
	  }
	free2int(vijA);
	free2int(vijB);
	return;
}

void RayPathFieldFrechet(Ray *ray, basisfun_para basis, int nodex, int nodez, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, DField *f)
{
	int i;
	int ix, iz;
	float ds, dt;
	float dx, dz, dl;//distance between each step
	float Px, Pz;
	double A2;
	double sin_fai, cos_fai;
	double delta_w1, delta_w2, delta_w3, delta_w4;
	double delta_f1, delta_f2, delta_f3, delta_f4;
	double delta_o1, delta_o2, delta_o3, delta_o4;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	double s_tv, s_xv, s_zv, s_pxv, s_pzv;
	double s_te, s_xe, s_ze, s_pxe, s_pze;
	double s_td, s_xd, s_zd, s_pxd, s_pzd;
	double e_xv, e_zv;
	double e_xe, e_ze;
	double e_xd, e_zd;
	double DVpertux, DVpertuz, V_pertu;
	double bepos, enda, dhdx, dhdz;
	double v2, v3, v4;
	double Px2, Pz2;
	double bx, bz, dbx, dbz;

	s_tv = 0; s_xv = 0; s_zv = 0; s_pxv = 0; s_pzv = 0;
	s_te = 0; s_xe = 0; s_ze = 0; s_pxe = 0; s_pze = 0;
	s_td = 0; s_xd = 0; s_zd = 0; s_pxd = 0; s_pzd = 0;
	for (i = 0; i < ray->nrs-1; i++) {
		ds = ray->rs[i + 1].s - ray->rs[i].s;
		dt = ray->rs[i + 1].t - ray->rs[i].t;
		dx = ray->rs[i + 1].x - ray->rs[i].x;
		dz = ray->rs[i + 1].z - ray->rs[i].z;
		dl = sqrt(dx * dx + dz * dz);
		sin_fai = dx / dl;
		cos_fai = dz / dl;
		vel2Interp(geo2dv, vel, ray->rs[i].x, ray->rs[i].z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
		vel2Interp(geo2dv, epos, ray->rs[i].x, ray->rs[i].z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
		vel2Interp(geo2dv, delta, ray->rs[i].x, ray->rs[i].z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
		A2 = pow(sin_fai * cos_fai, 2) / (1 + 2 * delta0) + pow(sin_fai, 4) / (1 + 2 * epos0) + pow(cos_fai, 2);

		ix = (int)(roundf((ray->rs[i].x - geo2dv.fx) / geo2dv.dx));
		iz = (int)(roundf((ray->rs[i].z - geo2dv.fz) / geo2dv.dz));
		bx = basis.bsplinex[nodex][ix];
		dbx = basis.dbsplinex[nodex][ix];
		bz = basis.bsplinez[nodez][iz];
		dbz = basis.dbsplinez[nodez][iz];

		DVpertux = dbx * bz * 1.0 + bx * dbz * 0.0; //theta=0, means DV/DX
		DVpertuz = dbx * bz * 0.0 + bx * dbz * 1.0; //theta=90,means DV/DZ
		V_pertu = bx * bz;
		Px = ray->rs[i].px;
		Pz = ray->rs[i].pz;
		Px2 = Px * Px;
		Pz2 = Pz * Pz;
		bepos = 1.0 + 2 * epos0;
		enda = epos0 - delta0;
		v2 = vel0 * vel0;
		v3 = vel0 * v2;
		v4 = v2 * v2;

		delta_w1 = 2*vel0*Px*(bepos-4*enda*v2*Pz2)*V_pertu;
		delta_w2 = 2*vel0*Pz*(1-4*enda*v2*Px2)*V_pertu;
		delta_w3 = ((12*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*dvdx-2*vel0*Px2*(1-2*v2*Pz2)*dedx-4*v3*Px2*Pz2*deldx)*V_pertu + (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*DVpertux;
		delta_w4 = ((12*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*dvdz-2*vel0*Px2*(1-2*v2*Pz2)*dedz-4*v3*Px2*Pz2*deldz)*V_pertu + (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*DVpertuz;
		
		delta_f1 = 2*v2*Px*(1-v2*Pz2)*V_pertu;
		delta_f2 = -2*v4*Px2*Pz*V_pertu;
		delta_f3 = 2*vel0*Px2*(2*v2*Pz2-1)*dvdx*V_pertu+v2*Px2*(v2*Pz2-1)*DVpertux;
		delta_f4 = 2*vel0*Px2*(2*v2*Pz2-1)*dvdz*V_pertu+v2*Px2*(v2*Pz2-1)*DVpertuz;
		
		delta_o1 = 2*v4*Px*Pz2*V_pertu;
		delta_o2 = 2*v4*Px2*Pz*V_pertu;
		delta_o3 = -4*v3*Px2*Pz2*dvdx*V_pertu-v4*Px2*Pz2*DVpertux;
		delta_o4 = -4*v3*Px2*Pz2*dvdz*V_pertu-v4*Px2*Pz2*DVpertuz;
	
		s_tv += -dt * V_pertu / vel0;
		s_xv  += ds * (+ray->rs[i+1].Q[2][2] * delta_w1 + ray->rs[i+1].Q[3][2] * delta_w2 - ray->rs[i+1].Q[0][2] * delta_w3 - ray->rs[i+1].Q[1][2] * delta_w4);
		s_zv  += ds * (+ray->rs[i+1].Q[2][3] * delta_w1 + ray->rs[i+1].Q[3][3] * delta_w2 - ray->rs[i+1].Q[0][3] * delta_w3 - ray->rs[i+1].Q[1][3] * delta_w4);
		s_pxv += ds * (-ray->rs[i+1].Q[2][0] * delta_w1 - ray->rs[i+1].Q[3][0] * delta_w2 + ray->rs[i+1].Q[0][0] * delta_w3 + ray->rs[i+1].Q[1][0] * delta_w4);
		s_pzv += ds * (-ray->rs[i+1].Q[2][1] * delta_w1 - ray->rs[i+1].Q[3][1] * delta_w2 + ray->rs[i+1].Q[0][1] * delta_w3 + ray->rs[i+1].Q[1][1] * delta_w4);
		
		s_te  += -dt * pow(sin_fai, 4) * V_pertu / pow((1 + 2 * epos0), 2) / A2;
		s_xe  += ds * (+ray->rs[i+1].Q[2][2] * delta_f1 + ray->rs[i+1].Q[3][2] * delta_f2 - ray->rs[i+1].Q[0][2] * delta_f3 - ray->rs[i+1].Q[1][2] * delta_f4);
		s_ze  += ds * (+ray->rs[i+1].Q[2][3] * delta_f1 + ray->rs[i+1].Q[3][3] * delta_f2 - ray->rs[i+1].Q[0][3] * delta_f3 - ray->rs[i+1].Q[1][3] * delta_f4);
		s_pxe += ds * (-ray->rs[i+1].Q[2][0] * delta_f1 - ray->rs[i+1].Q[3][0] * delta_f2 + ray->rs[i+1].Q[0][0] * delta_f3 + ray->rs[i+1].Q[1][0] * delta_f4);
		s_pze += ds * (-ray->rs[i+1].Q[2][1] * delta_f1 - ray->rs[i+1].Q[3][1] * delta_f2 + ray->rs[i+1].Q[0][1] * delta_f3 + ray->rs[i+1].Q[1][1] * delta_f4);
		
		s_td += -dt * pow(sin_fai * cos_fai, 2) * V_pertu / pow((1 + 2 * delta0), 2) / A2;
		s_xd  += ds * (+ray->rs[i+1].Q[2][2] * delta_o1 + ray->rs[i+1].Q[3][2] * delta_o2 - ray->rs[i+1].Q[0][2] * delta_o3 - ray->rs[i+1].Q[1][2] * delta_o4);
		s_zd  += ds * (+ray->rs[i+1].Q[2][3] * delta_o1 + ray->rs[i+1].Q[3][3] * delta_o2 - ray->rs[i+1].Q[0][3] * delta_o3 - ray->rs[i+1].Q[1][3] * delta_o4);
		s_pxd += ds * (-ray->rs[i+1].Q[2][0] * delta_o1 - ray->rs[i+1].Q[3][0] * delta_o2 + ray->rs[i+1].Q[0][0] * delta_o3 + ray->rs[i+1].Q[1][0] * delta_o4);
		s_pzd += ds * (-ray->rs[i+1].Q[2][1] * delta_o1 - ray->rs[i+1].Q[3][1] * delta_o2 + ray->rs[i+1].Q[0][1] * delta_o3 + ray->rs[i+1].Q[1][1] * delta_o4);
	}
	e_xv  = ray->rs[ray->nrs-1].Q[0][0] * s_xv + ray->rs[ray->nrs-1].Q[0][1] * s_zv + ray->rs[ray->nrs-1].Q[0][2] * s_pxv + ray->rs[ray->nrs-1].Q[0][3] * s_pzv;
	e_zv  = ray->rs[ray->nrs-1].Q[1][0] * s_xv + ray->rs[ray->nrs-1].Q[1][1] * s_zv + ray->rs[ray->nrs-1].Q[1][2] * s_pxv + ray->rs[ray->nrs-1].Q[1][3] * s_pzv;
	
	e_xe  = ray->rs[ray->nrs-1].Q[0][0] * s_xe + ray->rs[ray->nrs-1].Q[0][1] * s_ze + ray->rs[ray->nrs-1].Q[0][2] * s_pxe + ray->rs[ray->nrs-1].Q[0][3] * s_pze;
	e_ze  = ray->rs[ray->nrs-1].Q[1][0] * s_xe + ray->rs[ray->nrs-1].Q[1][1] * s_ze + ray->rs[ray->nrs-1].Q[1][2] * s_pxe + ray->rs[ray->nrs-1].Q[1][3] * s_pze;
	
	e_xd  = ray->rs[ray->nrs-1].Q[0][0] * s_xd + ray->rs[ray->nrs-1].Q[0][1] * s_zd + ray->rs[ray->nrs-1].Q[0][2] * s_pxd + ray->rs[ray->nrs-1].Q[0][3] * s_pzd;
	e_zd  = ray->rs[ray->nrs-1].Q[1][0] * s_xd + ray->rs[ray->nrs-1].Q[1][1] * s_zd + ray->rs[ray->nrs-1].Q[1][2] * s_pxd + ray->rs[ray->nrs-1].Q[1][3] * s_pzd;
	
	vel2Interp(geo2dv, vel, ray->rs[ray->nrs-1].x, ray->rs[ray->nrs-1].z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
	vel2Interp(geo2dv, epos, ray->rs[ray->nrs-1].x, ray->rs[ray->nrs-1].z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
	vel2Interp(geo2dv, delta, ray->rs[ray->nrs-1].x, ray->rs[ray->nrs-1].z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
	dhdx = ray->rs[ray->nrs-1].px*(1+2*epos0-2*(epos0-delta0)*vel0*vel0*ray->rs[ray->nrs-1].pz*ray->rs[ray->nrs-1].pz);
	dhdz = ray->rs[ray->nrs-1].pz*(1-2*(epos0-delta0)*vel0*vel0*ray->rs[ray->nrs-1].px*ray->rs[ray->nrs-1].px);

	f->tv = s_tv / wt;
	f->te = s_te / wt;
	f->td = s_td / wt;
	f->xv = (e_xv - dhdx * e_zv / dhdz) / wx;
	f->xe = (e_xe - dhdx * e_ze / dhdz) / wx;
	f->xd = (e_xd - dhdx * e_zd / dhdz) / wx;
	f->pv = (ray->rs[ray->nrs-1].Q[2][0] * s_xv + ray->rs[ray->nrs-1].Q[2][1] * s_zv + ray->rs[ray->nrs-1].Q[2][2] * s_pxv + ray->rs[ray->nrs-1].Q[2][3] * s_pzv) / wp;
	f->pe = (ray->rs[ray->nrs-1].Q[2][0] * s_xe + ray->rs[ray->nrs-1].Q[2][1] * s_ze + ray->rs[ray->nrs-1].Q[2][2] * s_pxe + ray->rs[ray->nrs-1].Q[2][3] * s_pze) / wp;
	f->pd = (ray->rs[ray->nrs-1].Q[2][0] * s_xd + ray->rs[ray->nrs-1].Q[2][1] * s_zd + ray->rs[ray->nrs-1].Q[2][2] * s_pxd + ray->rs[ray->nrs-1].Q[2][3] * s_pzd) / wp;
	return;
}

void RayPathZone(Geo2d geo2dv, bspline_para bpx, bspline_para bpz, Ray *ray, int **vij)
{
	int i;
	int k, l;
	int ks, ls;
	int m;

	m = bpx.K;

	zero2int(vij, bpz.Norigin, bpx.Norigin);
	for(i=0; i<ray->nrs; i++){
		ks = (ray->rs[i].x - geo2dv.fx) / bpx.interval;
		ls = (ray->rs[i].z - geo2dv.fz) / bpz.interval;
		if(ks < m) ks = m;
		if(ls < m) ls = m;
		if(ks > bpx.Norigin -1) ks = bpx.Norigin -1;
		if(ls > bpz.Norigin -1) ls = bpz.Norigin -1;
		for(k=ks; k>=ks-m; k--)
		  for(l=ls; l>=ls-m; l--){
			  vij[k][l] = 1;
		  }
	}
	return;
}

//calculate Ax=b
void cal_Newton_gradient(float *b,float **A,float *x,int Nmodel)
{
	int i,j,k;
	double *Ax,*pa,*p0,*r0,*r1;
	double ba,bb,bs;
	double ka,kb,ks;
	Ax=ealloc1double(Nmodel);
	pa=ealloc1double(Nmodel);
	p0=ealloc1double(Nmodel);
	r0=ealloc1double(Nmodel);
	r1=ealloc1double(Nmodel);
	
	//calculate the misfit r0=b-Ax
	for(i=0;i<Nmodel;i++){
		Ax[i]=0.0;
		for(j=0;j<Nmodel;j++)
		  Ax[i]+=A[j][i]*x[j];
		r0[i]=b[i]-Ax[i];
	}
	//p0=r0
	for(i=0;i<Nmodel;i++)
	  p0[i]=r0[i];
	
	for(k=0;k<6667;k++){
		/*a0=r0'*r0/(p0'*A*p0)*/
		ka=0;
		for(i=0;i<Nmodel;i++)
		  ka+=r0[i]*r0[i];//ka=r0'*r0
		//kb=p0'*A*p0
		for(i=0;i<Nmodel;i++){
			pa[i]=0.0;
			for(j=0;j<Nmodel;j++)
			  pa[i]+=A[i][j]*p0[j];
		}
		kb=0;
		for(i=0;i<Nmodel;i++)
		  kb+=pa[i]*p0[i];
		ks=ka/kb;
		/*update x0*/
		for(i=0;i<Nmodel;i++)
		  x[i]+=ks*p0[i];
		
		/*r1=r0-a0*A*p0*/
		for(i=0;i<Nmodel;i++){
			pa[i]=0.0;
			for(j=0;j<Nmodel;j++)
			  pa[i]+=A[j][i]*p0[j];
			r1[i]=r0[i]-ks*pa[i];
		}

		/*b0=r1'*r1/(r0'*r0)*/
		ba=0; bb=0;
		for(i=0;i<Nmodel;i++){
			ba+=r1[i]*r1[i];
			bb+=r0[i]*r0[i];
		}
		bs=ba/bb;
		/*p1=r1+b0*p0*/
		for(i=0;i<Nmodel;i++)
		  p0[i]=r1[i]+bs*p0[i];

		for(i=0;i<Nmodel;i++)
		  r0[i]=r1[i];
	}
	/*calculate the norm2 misfit*/
	float misfit=0;
	for(i=0;i<Nmodel;i++){
		Ax[i]=0.0;
		for(j=0;j<Nmodel;j++)
		  Ax[i]+=A[j][i]*x[j];
		r0[i]=b[i]-Ax[i];
	}
	for(i=0;i<Nmodel;i++)
	  misfit+=r0[i]*r0[i];
	printf("calculate newton grad misfit is %e\n",misfit);
	free1double(Ax);
	free1double(pa);
	free1double(p0);
	free1double(r0);
	free1double(r1);
}

