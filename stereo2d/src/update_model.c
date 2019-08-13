#include "mpi.h"
#include "su.h"
#include <stdlib.h>
#include "stereo_struct.h"
#include "wx.h"
#include "RaySub.h"
#include "bsplines.h"
#include "sparse_matrix.h"
#include "stereo_subroutines.h"
#include "forward_modeling.h"
#include "data_fitting_err.h"
#include "update_model.h"

void fedcpy(BGfield **vel_tmp, BGfield **vel, int nx, int nz)
{
	int i, j;
	for (i = 0; i < nx; i++)
		for (j = 0; j < nz; j++)
			vel_tmp[i][j] = vel[i][j];
}

void update_model(ModelSpace *mp0, ModelSpace *mp1,DataSpace *dp_cal, DataSpace *dp_true, bspline_para bpx, bspline_para bpz, basisfun_para basis, float *grad, float *dstep, float *misfit, float *error, int iter, int ndip_iter, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float vel_min, float vel_max, float epos_min, float epos_max, float delta_min, float delta_max, float **vcount, int myid, int numprocs)
{
	int i, j;
	int ndata;
	int nbspline;
	float weight;
	float acc_error;
	int vel_min_sum, vel_max_sum;
	float **vij, **vel_field_update;
	int epos_min_sum, epos_max_sum;
	float **eij, **epos_field_update;
	int delta_min_sum, delta_max_sum;
	float **dij, **delta_field_update;
	BGfield **vel_tmp, **epos_tmp, **delta_tmp;
	FILE *fv,*fe,*fd;

	ndata = dp_true->ndata;
	nbspline = bpx.Norigin * bpz.Norigin;
	vij = ealloc2float(bpz.Norigin, bpx.Norigin);
	vel_field_update = ealloc2float(bpz.Nfun, bpx.Nfun);
	vel_tmp = (BGfield**)ealloc2(bpz.Nfun, bpx.Nfun, sizeof(BGfield));
	eij = ealloc2float(bpz.Norigin, bpx.Norigin);
	epos_field_update = ealloc2float(bpz.Nfun, bpx.Nfun);
	epos_tmp = (BGfield**)ealloc2(bpz.Nfun, bpx.Nfun, sizeof(BGfield));
	dij = ealloc2float(bpz.Norigin, bpx.Norigin);
	delta_field_update = ealloc2float(bpz.Nfun, bpx.Nfun);
	delta_tmp = (BGfield**)ealloc2(bpz.Nfun, bpx.Nfun, sizeof(BGfield));

	/*---- calculate model parameter weight ----*/
	if (iter >= ndip_iter) {
		for (i = 0; i < mp0->vnx; i++) 
		  for (j = 0; j < mp0->vnz; j++) {
			  // save the update vij eij nij
			  vij[i][j] = grad[i * mp0->vnz + j] * vscale * (*dstep);
			  eij[i][j] = grad[nbspline + i * mp0->vnz + j] * escale * (*dstep);
			  dij[i][j] = grad[2 * nbspline + i * mp0->vnz + j] * dscale * (*dstep);
		  }
		cal_vel_field(bpx, bpz, basis, vij, vel_field_update);
		cal_vel_field(bpx, bpz, basis, eij, epos_field_update);
		cal_vel_field(bpx, bpz, basis, dij, delta_field_update);
		//save the previous full velocity field
		fedcpy(vel_tmp, vel, bpx.Nfun, bpz.Nfun);
		fedcpy(epos_tmp, epos, bpx.Nfun, bpz.Nfun);
		fedcpy(delta_tmp, delta, bpx.Nfun, bpz.Nfun);
	}
	for (weight = 1; weight > 0; weight -= 0.05) {
		if (iter >= ndip_iter) {
			for (i = 0; i < bpx.Nfun; i++) 
			  for (j = 0; j < bpz.Nfun; j++) {
				  vel[i][j].u = vel_tmp[i][j].u + weight * vel_field_update[i][j];
				  epos[i][j].u = epos_tmp[i][j].u + weight * epos_field_update[i][j];
				  delta[i][j].u = delta_tmp[i][j].u + weight * delta_field_update[i][j];
			  }
			/* check if velocity is reasonable */
			Vel_Statistics(vel, bpx.Nfun, bpz.Nfun, vel_min, vel_max, &vel_min_sum, &vel_max_sum);
			Vel_Statistics(epos, bpx.Nfun, bpz.Nfun, epos_min, epos_max, &epos_min_sum, &epos_max_sum);
			Vel_Statistics(delta, bpx.Nfun, bpz.Nfun, delta_min, delta_max, &delta_min_sum, &delta_max_sum);
			if(myid==2){
				fv=fopen("vel.dat","w");
				fe=fopen("epos.dat","w");
				fd=fopen("delta.dat","w");
				for(i=0;i<bpx.Nfun;i++)
				  for(j=0;j<bpz.Nfun;j++){
					  fwrite(&vel[i][j].u,sizeof(float),1,fv);
					  fwrite(&epos[i][j].u,sizeof(float),1,fe);
					  fwrite(&delta[i][j].u,sizeof(float),1,fd);
				  }
				fclose(fv);
				fclose(fe);
				fclose(fd);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		/* update xcor, zcor, theta_s, theta_r */
		if (iter >= ndip_iter) {
			dv2(geo2dv, vel);
			dv2(geo2dv, epos);
			dv2(geo2dv, delta);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		forward_modeling(mp0, dp_cal, geo2dv, vel, epos, delta, vcount, myid, numprocs);
		cal_cost_fun(&acc_error, misfit, dp_true, dp_cal);
		if (myid == 0) printf("myid:%d weight=%f now:acc_error=%e error_iter=%e per=%e n=%d\n", myid, weight, acc_error, error[iter], acc_error / ndata, ndata);
		/* accept updated model */
		if (acc_error < error[iter] || fabs(weight - 0.05) < 1.0E-6 || ((iter == ndip_iter) && iter != 0)) {
			//update bspline coefficients 
			*dstep *= weight;
			error[iter + 1] = acc_error;
			if (iter < ndip_iter)
			  break;
			for (i = 0; i < mp0->vnx; i++)
			  for (j = 0; j < mp0->vnz; j++){
				  mp0->vij[i][j] += vij[i][j] * weight;
				  mp0->eij[i][j] += eij[i][j] * weight;
				  mp0->dij[i][j] += dij[i][j] * weight;
			  }
			break;
		}
	}// end of weight loop
	free2float(vij);
	free2float(vel_field_update);
	free2float(eij);
	free2float(epos_field_update);
	free2float(dij);
	free2float(delta_field_update);
	free(vel_tmp);
	free(epos_tmp);
	free(delta_tmp);
}
