#pragma once
#ifndef BSPLINES_H_
#define BSPLINES_H_
typedef struct 
{
	int Norigin;
	int K;
	int Nnode;
	int interval;
	float dsample;
	int Nfun;
	float *node;
}bspline_para;

typedef struct 
{
	float **bsplinex, **dbsplinex, **ddbsplinex;
	float **bsplinez, **dbsplinez, **ddbsplinez;
}basisfun_para;

void bparameter_initial(bspline_para *bp, float interval, float dsamle, int nsample);
void node_initial(bspline_para *bp);
void basis_spline(bspline_para bp, float **bspline, float **dbspline, float **ddbspline);
void basis_analytical_value(bspline_para bp, float *bspline, float *dbspline, float *ddbspline, float x);
void cal_vel_field(bspline_para bpx, bspline_para bpz, basisfun_para basis, float **vij, float **vel_out);
void bspline_interpolation(float **vel, bspline_para bpx, bspline_para bpz, float **vij, basisfun_para basis);

#endif // !BSPLINES_H_
