#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include "bsplines.h"
#include "sparse_matrix.h"
#include "regu.h"

void cal_damp_regu(float **Hess, int ndata, bspline_para bpx, bspline_para bpz, float regu_ray_damp, float regu_damp, int iter, int ndip_iter, float **regu_vij)
{
/*
 *add small damp on Hessian 
 */
	int i, j;
	float regu, reguv;
	int nbspline;
	int col, diagonal;
	nbspline = bpx.Norigin*bpz.Norigin;

	if (iter < ndip_iter)
		regu = regu_ray_damp;
	else
		regu = regu_damp;

	//regu for ray segments
//	for (i = 0; i < 4 * ndata; i++)
//	  Hess[i][i]+=regu*regu;
	//regu for Bspline coes
//	col = 4 * ndata;
	col=0;
	if (iter >= ndip_iter) {
		for (i = 0; i < bpx.Norigin; i++) 
		  for (j = 0; j < bpz.Norigin; j++) {
			  diagonal = i * bpx.Norigin + j;
			  reguv = regu_damp * regu_vij[i][j] * vscale;
			  Hess[col+diagonal][col+diagonal]+=pow(reguv,2);  
			  Hess[col+nbspline+diagonal][col+nbspline+diagonal]+=pow(reguv,2);
			  Hess[col+2*nbspline+diagonal][col+2*nbspline+diagonal]+=pow(reguv,2);
		  }
	}
	return;
}

void cal_homo_regu(bspline_para bpx, bspline_para bpz, float regu_homox, float regu_homoz, float **Hess, int ndata)
{
//smooth the update
	int i, j, m, n;
	int col[3];
	float coe;
	float value[3];
	float factorx = 10.0, factorz = 10.0;
	int nbspline;
	nbspline = bpx.Norigin*bpz.Norigin;

	/*--homo regu along x direction--*/
	for (i = 0; i < bpz.Norigin; i++) 
		for (j = 0; j < bpx.Norigin - 2; j++) {
			coe = regu_homox * vscale;
			value[0]= coe;
			value[1]=-2*coe;
			value[2]= coe;
			col[0] = j * bpz.Norigin + i;
			col[1] = (j + 1) * bpz.Norigin + i;
			col[2] = (j + 2) * bpz.Norigin + i;
			for (m = 0; m < 3; m++)
			  for (n = 0; n < 3; n++){
				  Hess[col[m]][col[n]]+=value[m]*value[n];
				  Hess[col[m]+nbspline][col[n]+nbspline]+=value[m]*value[n];
				  Hess[col[m]+2*nbspline][col[n]+2*nbspline]+=value[m]*value[n];
			  }
		}

	/*--homo regu along z direction--*/
	for (i = 0; i < bpx.Norigin; i++) 
		for (j = 0; j < bpz.Norigin - 2; j++) {
			coe = regu_homoz * vscale;
			value[0]= coe;
			value[1]=-2*coe;
			value[2]= coe;
			col[0] = i * bpz.Norigin + j;
			col[1] = i * bpz.Norigin + j + 1;
			col[2] = i * bpz.Norigin + j + 2;
			for (m = 0; m < 3; m++)
			  for (n = 0; n < 3; n++){
				  Hess[col[m]][col[n]]+=value[m]*value[n];
				  Hess[col[m]+nbspline][col[n]+nbspline]+=value[m]*value[n];
				  Hess[col[m]+2*nbspline][col[n]+2*nbspline]+=value[m]*value[n];
			  }
		}
}

void vij_statistic(bspline_para bpx, bspline_para bpz, float **regu_vij, float **vcount)
{
	int i, j, k, l;
	int ks, ls;
	float **vijcount;

	vijcount = ealloc2float(bpz.Norigin, bpx.Norigin);
	zero2float(vijcount, bpz.Norigin, bpx.Norigin);

	for (i = 0; i < bpx.Nfun; i++) 
		for (j = 0; j < bpz.Nfun; j++) {
			ks = i * bpx.dsample / bpx.interval;
			ls = j * bpz.dsample / bpz.interval;

			if (ks < 3) ks = 3;
			if (ls < 3) ls = 3;
			if (ks > bpx.Norigin - 1) ks = bpx.Norigin - 1;
			if (ls > bpz.Norigin - 1) ls = bpz.Norigin - 1;

			for (k = ks; k >= ks - 3; k--) 
				for (l = ls; l >= ls - 3; l--) {
					vijcount[k][l] += vcount[i][j];
				}
			
		}
	

	for (k = 0; k < bpx.Norigin; k++) 
		for (l = 0; l < bpz.Norigin; l++) {
			if (vijcount[k][l] > 0)
				regu_vij[k][l] = 1;
			else
				regu_vij[k][l] = 10;
		}
	free2float(vijcount);
}

void dip_statistic(ModelSpace *mp0, Geo2d geo2dv, bspline_para bpx, bspline_para bpz, float **dip_sum, float **dip_count)
{
	int i, j, k, l;
	int ks, ls;
	zero2float(dip_sum, bpz.Norigin, bpx.Norigin);
	zero2float(dip_count, bpz.Norigin, bpx.Norigin);

	for (i = 0; i < mp0->ndata; i++) {
		ks = (mp0->m[i].xcor - geo2dv.fx) / bpx.interval;
		ls = (mp0->m[i].zcor - geo2dv.fz) / bpz.interval;

		if (ks < 3) ks = 3;
		if (ls < 3) ls = 3;
		if (ks > bpx.Norigin - 1) ks = bpx.Norigin - 1;
		if (ls > bpz.Norigin - 1) ls = bpz.Norigin - 1;

		for (k = ks; k >= ks - 3; k--) 
			for (l = ls; l >= ls - 3; l--) {
				dip_sum[k][l] += (mp0->m[i].theta_s + mp0->m[i].theta_r) / 2.0;
				dip_count[k][l]++;
			}
	}
	for (i = 0; i < bpx.Norigin; i++)
		for (j = 0; j < bpz.Norigin; j++) {
			if (dip_count[i][j] >= 1) 
				dip_sum[i][j] /= dip_count[i][j];
		}
}

void dip_statistic_small(ModelSpace *mp0, Geo2d geo2dv, bspline_para bpx, bspline_para bpz, float **dip_sum, float **dip_count)
{
	int i, j;
	int ix, iz;
	zero2float(dip_sum, bpz.Nfun, bpx.Nfun);
	zero2float(dip_count, bpz.Nfun, bpx.Nfun);

	for (i = 0; i < mp0->ndata; i++) {
		ix = (mp0->m[i].xcor - geo2dv.fx) / geo2dv.dx + 0.5;
		iz = (mp0->m[i].zcor - geo2dv.fz) / geo2dv.dx + 0.5;
		dip_sum[ix][iz] += (mp0->m[i].theta_s + mp0->m[i].theta_r) / 2.0;
		dip_count[ix][iz]++;
	}
	for (i = 0; i < bpx.Nfun; i++)
		for (j = 0; j < bpz.Nfun; j++) {
			if (dip_count[i][j] >= 1) 
				dip_sum[i][j] /= dip_count[i][j];
		}
}

int none_zero_count(float **data, int nx, int nz)
{
	int i, j;
	int count;
	count = 0;
	for (i = 0; i < nx; i++) 
		for (j = 0; j < nz; j++) {
			if (fabsf(data[i][j]) > 1e-7)
				count++;
		}
	return count;
}
