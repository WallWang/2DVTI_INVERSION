#pragma once
#ifndef STEREO_STRUCT_H_
#define STEREO_STRUCT_H_
#define MAX_RAYSTEP 99999
static float wt = 1.0e-3;
static float wp = 1.0e-6;
static float wx = 1.0;
static float xcor_scale = 1.0e3;
static float zcor_scale = 1.0e3;
static float sita_scale = 0.5;
static float vscale = 1.0E3;
static float escale = 1.0;
static float dscale = 1.0;

typedef struct {
	float xcor;
	float zcor;
	float theta_s;//theta to source, unit in angle degree
	float theta_r;//theta to receiver, unit in angle degree
}Model_para;

typedef struct {
	int vnx, vnz;
	int ndata;
	float **vij;
	float **eij;
	float **dij;
	Model_para *m;
}ModelSpace;

typedef struct {
	float sx, rx;
	float spx, rpx;
	float t;
}Data_para;

typedef struct {
	int ndata;
	Data_para *d;
}DataSpace;

typedef struct {
	float t;
	float s;
	float x, z;
	float px, pz;
	float Q[4][4];
}RayStep;

typedef struct {
	int nrs;
	RayStep *rs;
}Ray;

typedef struct Geo2DStruct{
	int nx;
	float fx;
	float dx;
	int nz;
	float fz;
	float dz;
	float dt;
	float xmin, xmax, zmin, zmax;
}Geo2d;

typedef struct {
	float u;
	double uxx;
	double uxz;
	double uzz;
}BGfield;

typedef struct {
	float x, z;
	float px, pz;
	float t;
}Source;

typedef struct{
	float tv, pv, xv;
	float te, pe, xe;
	float td, pd, xd;
}DField;
#endif // !STEREO_STRUCT_H_
