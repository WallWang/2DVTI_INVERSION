#include <stdlib.h>
#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include "RaySub.h"

void draw_raypath(FILE *fp_ray_path, Ray ray)
{
	int i;
	for (i = 0; i < ray.nrs; i++) {
		if (i != 0)
			fprintf(fp_ray_path, "PA %15f %15f\n", ray.rs[i].x, ray.rs[i].z);
		else
			fprintf(fp_ray_path, "MA %15f %15f\n", ray.rs[i].x, ray.rs[i].z);
	}
}

void dv2(Geo2d geo2dv, BGfield **vel)
/*calculate second derivatives from a 2D velocity grid by finite defference*/
{
	int i, j;
	int nx, nz;
	float dx, dz;
	nx = geo2dv.nx;
	nz = geo2dv.nz;
	dx = geo2dv.dx;
	dz = geo2dv.dz;
	/*--initialize--*/
	for (i = 0; i < nx; i++)
		for (j = 0; j < nz; j++) {
			vel[i][j].uxx = 0.0;
			vel[i][j].uxz = 0.0;
			vel[i][j].uzz = 0.0;
		}
	/*--finite difference--*/
	for (i = 1; i < nx - 1; i++)
		for (j = 1; j < nz - 1; j++) {
			vel[i][j].uxx = (vel[i + 1][j].u - 2.0*vel[i][j].u + vel[i - 1][j].u) / dx / dx;
			vel[i][j].uxz = (vel[i + 1][j + 1].u - vel[i + 1][j - 1].u - vel[i - 1][j + 1].u + vel[i - 1][j - 1].u) / 4.0 / dx / dz;
			vel[i][j].uzz = (vel[i][j + 1].u - 2.0*vel[i][j].u + vel[i][j - 1].u) / dz / dz;
		}
	//for boundary
	for (i = 0; i < nx; i++) {
		vel[i][0].uxx = vel[i][1].uxx;
		vel[i][0].uxz = vel[i][1].uxz;
		vel[i][0].uzz = vel[i][1].uzz;
		vel[i][nz - 1].uxx = vel[i][nz - 2].uxx;
		vel[i][nz - 1].uxz = vel[i][nz - 2].uxz;
		vel[i][nz - 1].uzz = vel[i][nz - 2].uzz;
	}
	for (j = 1; j < nz - 1; j++) {
		vel[0][j].uxx = vel[1][j].uxx;
		vel[0][j].uxz = vel[1][j].uxz;
		vel[0][j].uzz = vel[1][j].uzz;
		vel[nx - 1][j].uxx = vel[nx - 2][j].uxx;
		vel[nx - 1][j].uxz = vel[nx - 2][j].uxz;
		vel[nx - 1][j].uzz = vel[nx - 2][j].uzz;
	}
}

void vel2Interp(Geo2d geo2dv, BGfield **vel, float x, float z, double *u, double *ux, double *uz, double *uxx, double *uxz, double *uzz)
/*vel2Interp - Function to support interpolation of velocity and its derivatives*/
{
	int jx, jz;
	int nx = geo2dv.nx, nz = geo2dv.nz;
	float fx = geo2dv.fx, fz = geo2dv.fz;
	float dx = geo2dv.dx, dz = geo2dv.dz;
	double ax, az, sx, sz, sxx, szz, a0, a1, a2, a3;

	/* determine interpolate coefficients */
	ax = (x - fx) / dx;
	jx = ax;//jx is integer
	if (jx < 0) jx = 0;
	if (jx > nx - 2) jx = nx - 2;
	sx = ax - jx;
	az = (z - fz) / dz;
	jz = az;
	if (jz < 0) jz = 0;
	if (jz > nz - 2) jz = nz - 2;
	sz = az - jz;

	sxx = 0.5*sx*(1.0 - sx)*dx*dx;
	szz = 0.5*sz*(1.0 - sz)*dz*dz;

	a0 = (1.0 - sx)*(1.0 - sz);
	a1 = (1.0 - sx)*sz;
	a2 = sx * (1.0 - sz);
	a3 = sx * sz;
	
	/*--Bilinear interpolation for second velocity derivative--*/
	*uxx = a0 * vel[jx][jz].uxx + a1 * vel[jx][jz + 1].uxx + a2 * vel[jx + 1][jz].uxx + a3 * vel[jx + 1][jz + 1].uxx;
	*uxz = a0 * vel[jx][jz].uxz + a1 * vel[jx][jz + 1].uxz + a2 * vel[jx + 1][jz].uxz + a3 * vel[jx + 1][jz + 1].uxz;
	*uzz = a0 * vel[jx][jz].uzz + a1 * vel[jx][jz + 1].uzz + a2 * vel[jx + 1][jz].uzz + a3 * vel[jx + 1][jz + 1].uzz;
	
	*u = a0 * vel[jx][jz].u + a1 * vel[jx][jz + 1].u + a2 * vel[jx + 1][jz].u + a3 * vel[jx + 1][jz + 1].u - (sxx * (*uxx) + szz * (*uzz));
	*ux = ((1.0 - sz)*(vel[jx + 1][jz].u - vel[jx][jz].u - sxx * (vel[jx + 1][jz].uxx - vel[jx][jz].uxx) - szz * (vel[jx + 1][jz].uzz
		- vel[jx][jz].uzz)) + sz * (vel[jx + 1][jz + 1].u - vel[jx][jz + 1].u - sxx * (vel[jx + 1][jz + 1].uxx - vel[jx][jz + 1].uxx)
		- szz * (vel[jx + 1][jz + 1].uzz - vel[jx][jz + 1].uzz))) / dx + (sx - 0.5)*dx*(*uxx);
	*uz = ((1.0 - sx)*(vel[jx][jz + 1].u - vel[jx][jz].u - sxx * (vel[jx][jz + 1].uxx - vel[jx][jz].uxx) - szz * (vel[jx][jz + 1].uzz
		- vel[jx][jz].uzz)) + sx * (vel[jx + 1][jz + 1].u - vel[jx + 1][jz].u - sxx * (vel[jx + 1][jz + 1].uxx - vel[jx + 1][jz].uxx)
		- szz * (vel[jx + 1][jz + 1].uzz - vel[jx + 1][jz].uzz))) / dz + (sz - 0.5)*dz*(*uzz);
}

void save_ray_path(Ray *ray, double *y)
{
	ray->rs[ray->nrs].x = y[0];
	ray->rs[ray->nrs].z = y[1];
	ray->rs[ray->nrs].px = y[2];
	ray->rs[ray->nrs].pz = y[3];
	ray->rs[ray->nrs].t = y[4]; //travel time
	ray->rs[ray->nrs].s = y[5]; //total integral step
	ray->rs[ray->nrs].Q[0][0] = y[6];
	ray->rs[ray->nrs].Q[0][1] = y[7];
	ray->rs[ray->nrs].Q[0][2] = y[8];
	ray->rs[ray->nrs].Q[0][3] = y[9];
	ray->rs[ray->nrs].Q[1][0] = y[10];
	ray->rs[ray->nrs].Q[1][1] = y[11];
	ray->rs[ray->nrs].Q[1][2] = y[12];
	ray->rs[ray->nrs].Q[1][3] = y[13];
	ray->rs[ray->nrs].Q[2][0] = y[14];
	ray->rs[ray->nrs].Q[2][1] = y[15];
	ray->rs[ray->nrs].Q[2][2] = y[16];
	ray->rs[ray->nrs].Q[2][3] = y[17];
	ray->rs[ray->nrs].Q[3][0] = y[18];
	ray->rs[ray->nrs].Q[3][1] = y[19];
	ray->rs[ray->nrs].Q[3][2] = y[20];
	ray->rs[ray->nrs].Q[3][3] = y[21];
	ray->nrs = ray->nrs + 1;
}

void out_of_model(int nEquations, double *y, double *y0, float xmin, float xmax, float zmin, float zmax, double err)
/*Interpolation location and slowness on model boundary*/
{
	int i;
	double x0 = -9999, z0 = -9999;
	double s0, s1;
	
	if (fabs(y[0] - y0[0]) < err) {
		if (y[3] < 0) {
			x0 = y[0];
			z0 = zmin;
		}
		else {
			x0 = y[0];
			z0 = zmax;
		}
	}
	else if (fabs(y[1] - y0[1]) < err) {
		if (y[2] < 0) {
			x0 = xmin;
			z0 = y[1];
		}
		else {
			x0 = xmax;
			z0 = y[1];
		}
	}
	else {
		if (y[0] < xmin) {
			x0 = xmin;
			z0 = (x0 * y[1] - x0 * y0[1] + y[0] * y0[1] - y[1] * y0[0]) / (y[0] - y0[0]);
		}
		else if (y[0] > xmax) {
			x0 = xmax;
			z0 = (x0 * y[1] - x0 * y0[1] + y[0] * y0[1] - y[1] * y0[0]) / (y[0] - y0[0]);
		}
		else if (y[1] < zmin) {
			z0 = zmin;
			x0 = (z0 * y[0] - z0 * y0[0] + y[1] * y0[0] - y[0] * y0[1]) / (y[1] - y0[1]);
		}
		else {
			z0 = zmax;
			x0 = (z0 * y[0] - z0 * y0[0] + y[1] * y0[0] - y[0] * y0[1]) / (y[1] - y0[1]);
		}
	}
	s0 = sqrt(pow((y0[0] - x0), 2) + pow((y0[1] - z0), 2));
	s1 = sqrt(pow((y[0] - y0[0]), 2) + pow((y[1] - y0[1]), 2));
	y[0] = x0;
	y[1] = z0;
	for (i = 2; i < nEquations; i++)
		y[i] = y0[i] + (y[i] - y0[i])*s0/s1;
}

void dynamic_ray_tracing(Ray *ray, Source sc, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float timeout)
/*
 * Input: init condition, observation system, velcoity, anisotropic parameters and ending condition
 * Calculation: Raypath and it's propogator matrix
 */
{
	int nEquations = 22;
	int i, j;
	double err = 1.0e-6;
	double z[22][2] = {{0.0}}, y[22], y0[22];
	float xmin = geo2dv.xmin, xmax = geo2dv.xmax;
	float zmin = geo2dv.zmin, zmax = geo2dv.zmax;
	ray->nrs = 0;

	//initial conditions
	y[0] = sc.x;
	y[1] = sc.z;
	y[2] = sc.px;
	y[3] = sc.pz;
	y[4] = 0;//travel time
	y[5] = 0;//integral step
	y[6] = 1.0;   y[7] = 0.0;   y[8] = 0.0;   y[9] = 0.0;
	y[10] = 0.0;  y[11] = 1.0;  y[12] = 0.0;  y[13] = 0.0;
	y[14] = 0.0;  y[15] = 0.0;  y[16] = 1.0;  y[17] = 0.0;
	y[18] = 0.0;  y[19] = 0.0;  y[20] = 0.0;  y[21] = 1.0;//propogator matrix init as unit matirx

	for (i = 0; i < MAX_RAYSTEP; i++) {
		//save raypath information
		save_ray_path(ray, &y[0]);
		for (j = 0; j < nEquations; j++) y0[j] = y[j];
		rkt1_for_dynamic(y, nEquations, &z[0][0], geo2dv, vel, epos, delta);
		//ray end condition 1: timeout
		if (y[4] >= timeout) {
			//save raypath information
			save_ray_path(ray, &y[0]);
			break;
		}
		//ray end condition 2: be out of model
		if (y[0] > xmax + err || y[1] > zmax + err || y[0] < xmin - err || y[1] < zmin - err) {
			out_of_model(nEquations, y, y0, xmin, xmax, zmin, zmax, err);
			//save raypath information
			save_ray_path(ray, &y[0]);
			break;
		}
	}
}

void rkt1_for_dynamic(double *y, int n, double *z, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	int i, j;
	double *b, *d;
	float a[4];
	b = ealloc1double(n);
	d = ealloc1double(n);
	a[0] = geo2dv.dt / 2.0; a[2] = geo2dv.dt;
	a[1] = geo2dv.dt / 2.0; a[3] = geo2dv.dt;
	for (i = 0; i < n; i++) z[i * 2] = y[i];
	func_for_dynamic(y, d, geo2dv, vel, epos, delta);
	for (i = 0; i < n; i++) b[i] = y[i];
	for (j = 0; j <= 2; j++) {
		for (i = 0; i < n; i++) {
			y[i] = z[i * 2] + a[j] * d[i];
			b[i] = b[i] + a[j + 1] * d[i] / 3.0;
		}
		func_for_dynamic(y, d, geo2dv, vel, epos, delta);
	}
	for (i = 0; i < n; i++)
		y[i] = b[i] + geo2dv.dt*d[i] / 6.0;
	for (i = 0; i < n; i++)
		z[i * 2 + 1] = y[i];
	free1double(b);
	free1double(d);
	return;
}

void func_for_dynamic(double *y, double *d, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	int i, j, k;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	double v2, v3, v4, bepos, enda;
	double Px2, Pz2;
	float x, z, Px, Pz;
	float Q[4][4],S[4][4];
	x = y[0]; Px = y[2];
	z = y[1]; Pz = y[3];
	for (i = 0; i < 4; i++) 
		for (j = 0; j < 4; j++) 
		  Q[i][j] = y[i * 4 + j + 6];
		
	vel2Interp(geo2dv, vel, x, z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
	vel2Interp(geo2dv, epos, x, z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
	vel2Interp(geo2dv, delta, x, z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
	v2 = vel0 * vel0;
	v3 = vel0 * v2;
	v4 = v2 * v2;
	bepos = 1.0 + 2 * epos0;
	enda = epos0 - delta0;
	Px2 = Px * Px;
	Pz2 = Pz * Pz;

	/* integration along dt ! */
	d[0] = v2*Px*(bepos-2*enda*v2*Pz2);        /// for X
	d[1] = v2*Pz*(1-2*enda*v2*Px2);         /// for Z
	d[2] = (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*dvdx-v2*Px2*(1-v2*Pz2)*dedx-v4*Px2*Pz2*deldx;
	d[3] = (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*dvdz-v2*Px2*(1-v2*Pz2)*dedz-v4*Px2*Pz2*deldz;
	d[4] = 1-2*enda*v4*Px2*Pz2;               /// for travel time
	d[5] = 1;
	
	S[0][0] = 2*vel0*Px*(bepos-4*enda*v2*Pz2)*dvdx + 2*v2*Px*(1-v2*Pz2)*dedx + 2*v4*Px*Pz2*deldx;
	S[0][1] = 2*vel0*Px*(bepos-4*enda*v2*Pz2)*dvdz + 2*v2*Px*(1-v2*Pz2)*dedz + 2*v4*Px*Pz2*deldz;
	S[0][2] = v2*(bepos-2*enda*v2*Pz2);
	S[0][3] = -4*enda*v4*Px*Pz;

	S[1][0] = 2*vel0*Pz*(1-4*enda*v2*Px2)*dvdx - 2*v4*Px2*Pz*dedx + 2*v4*Px2*Pz*deldx;
	S[1][1] = 2*vel0*Pz*(1-4*enda*v2*Px2)*dvdz - 2*v4*Px2*Pz*dedz + 2*v4*Px2*Pz*deldz;
	S[1][2] = -4*enda*v4*Px*Pz;
	S[1][3] = v2*(1-2*enda*v2*Px2);

	S[2][0] = (12*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*dvdx*dvdx - v2*Px2*(1-v2*Pz2)*epxx - 8*v3*Px2*Pz2*dvdx*deldx + (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*uxx - v4*Px2*Pz2*delxx - 4*vel0*Px2*(1-2*v2*Pz2)*dedx*dvdx;
	S[2][1] = (12*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*dvdx*dvdz - v2*Px2*(1-v2*Pz2)*epxz - 4*v3*Px2*Pz2*(dvdx*deldz+deldx*dvdz) + (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*uxz - v4*Px2*Pz2*delxz - 2*vel0*Px2*(1-2*v2*Pz2)*(dvdx*dedz+dedx*dvdz);
	S[2][2] = 2*vel0*Px*(4*enda*v2*Pz2-bepos)*dvdx - 2*v2*Px*(1-v2*Pz2)*dedx - 2*v4*Px*Pz2*deldx;
	S[2][3] = 2*vel0*Pz*(4*enda*v2*Px2-1)*dvdx + 2*v4*Px2*Pz*dedx - 2*v4*Px2*Pz*deldx;

	S[3][0] = (12*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*dvdz*dvdx - v2*Px2*(1-v2*Pz2)*epxz - 4*v3*Px2*Pz2*(dvdz*deldx+dvdx*deldz) + (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*uxz - v4*Px2*Pz2*delxz - 2*vel0*Px2*(1-2*v2*Pz2)*(dedz*dvdx+dedx*dvdz);
	S[3][1] = (12*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*dvdz*dvdz - v2*Px2*(1-v2*Pz2)*epzz - 8*v3*Px2*Pz2*deldz*dvdz + (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*uzz - v4*Px2*Pz2*delzz - 4*vel0*Px2*(1-2*v2*Pz2)*dvdz*dedz;
	S[3][2] = 2*vel0*Px*(4*enda*v2*Pz2-bepos)*dvdz - 2*v2*Px*(1-v2*Pz2)*dedz - 2*v4*Px*Pz2*deldz;
	S[3][3] = 2*vel0*Pz*(4*enda*v2*Px2-1)*dvdz + 2*v4*Px2*Pz*dedz - 2*v4*Px2*Pz*deldz;

	/* progagator matrix */
	for (i = 0; i < 4; i++) 
		for (j = 0; j < 4; j++) {
			d[i * 4 + j + 6] = 0;
			for (k = 0; k < 4; k++) {
				d[i * 4 + j + 6] += S[i][k] * Q[k][j];
			}
		}
}

void kinematic_ray_tracing(Ray *ray, Source sc, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float timeout)
{
/*
 * Input: init condition, observation system, velcoity, anisotropic parameters and ending condition
 * Calculation: Raypath
 */
	int nEquations = 6;
	double err = 1.0e-6;
	int i, j;
	double z[6][2] = {{0.0}}, y[6], y0[6];
	float xmin = geo2dv.xmin, xmax = geo2dv.xmax;
	float zmin = geo2dv.zmin, zmax = geo2dv.zmax;
	ray->nrs = 0;
	//initial conditions
	y[0] = sc.x;
	y[1] = sc.z;
	y[2] = sc.px;
	y[3] = sc.pz;
	y[4] = 0;
	y[5] = 0;
	for (i = 0; i < MAX_RAYSTEP; i++)
	{
		ray->rs[ray->nrs].x = y[0];
		ray->rs[ray->nrs].z = y[1];
		ray->rs[ray->nrs].px = y[2];
		ray->rs[ray->nrs].pz = y[3];
		ray->rs[ray->nrs].t = y[4];
		ray->rs[ray->nrs].s = y[5];
		ray->nrs = ray->nrs + 1;
		for (j = 0; j < nEquations; j++) y0[j] = y[j];
		rkt1_for_kinematic(y, nEquations, &z[0][0], geo2dv, vel, epos, delta);
		if (y[4] >= timeout) {
			ray->rs[ray->nrs].x = y[0];
			ray->rs[ray->nrs].z = y[1];
			ray->rs[ray->nrs].px = y[2];
			ray->rs[ray->nrs].pz = y[3];
			ray->rs[ray->nrs].t = y[4];
		    ray->rs[ray->nrs].s = y[5];
			ray->nrs = ray->nrs + 1;
			break;
		}
		/* out of the model, then passback this ray and sace the end point. */
		if (y[0] > xmax + err || y[1] > zmax + err || y[0] < xmin - err || y[1] < zmin - err) {
			out_of_model(nEquations, y, y0, xmin, xmax, zmin, zmax, err);
			ray->rs[ray->nrs].x = y[0];
			ray->rs[ray->nrs].z = y[1];
			ray->rs[ray->nrs].px = y[2];
			ray->rs[ray->nrs].pz = y[3];
			ray->rs[ray->nrs].t = y[4];
		    ray->rs[ray->nrs].s = y[5];
			ray->nrs = ray->nrs + 1;
			break;
		}
	}
}

void rkt1_for_kinematic(double *y, int n, double *z, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	int i, j;
	double *b, *d;
	float a[4];
	b = ealloc1double(n);
	d = ealloc1double(n);
	a[0] = geo2dv.dt / 2.0; a[2] = geo2dv.dt;
	a[1] = geo2dv.dt / 2.0; a[3] = geo2dv.dt;
	for (i = 0; i < n; i++) z[i * 2] = y[i];
	func_for_kinematic(y, d, geo2dv, vel, epos, delta);
	for (i = 0; i < n; i++) b[i] = y[i];
	for (j = 0; j <= 2; j++) {
		for (i = 0; i < n; i++) {
			y[i] = z[i * 2] + a[j] * d[i];
			b[i] = b[i] + a[j + 1] * d[i] / 3.0;
		}
		func_for_kinematic(y, d, geo2dv, vel, epos, delta);
	}
	for (i = 0; i < n; i++)
		y[i] = b[i] + geo2dv.dt*d[i] / 6.0;
	for (i = 0; i < n; i++)
		z[i * 2 + 1] = y[i];
	free1double(b);
	free1double(d);
	return;
}

void func_for_kinematic(double *y, double *d, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	double v2, v4, bepos, enda;
	float x, z, Px, Pz;
	double Px2, Pz2;

	x = y[0];   Px = y[2];
	z = y[1];   Pz = y[3];

	vel2Interp(geo2dv, vel, x, z, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
	vel2Interp(geo2dv, epos, x, z, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);
	vel2Interp(geo2dv, delta, x, z, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);
	
	v2 = vel0 * vel0;
	v4 = v2 * v2;
	bepos = 1.0 + 2 * epos0;
	enda = epos0 - delta0;
	Px2 = Px * Px;
	Pz2 = Pz * Pz;
	/* integration along dt ! */
	
	d[0] = v2*Px*(bepos-2*enda*v2*Pz2);        /// for X
	d[1] = v2*Pz*(1-2*enda*v2*Px2);         /// for Z
	d[2] = (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*dvdx-v2*Px2*(1-v2*Pz2)*dedx-v4*Px2*Pz2*deldx;
	d[3] = (4*enda*v2*Px2*Pz2-bepos*Px2-Pz2)*vel0*dvdz-v2*Px2*(1-v2*Pz2)*dedz-v4*Px2*Pz2*deldz;
	d[4] = 1-2*enda*v4*Px2*Pz2;               /// for travel time
	d[5] = 1;
	
}
