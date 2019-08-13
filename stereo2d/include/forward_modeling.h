#ifndef _FORWARD_MODELING_H
#define _FORWARD_MODELING_H


void DRayShot(ModelSpace *mp, Data_para *dp_cal, Ray *rayA, Ray *rayB, int i, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);

void forward_modeling(ModelSpace *mp1, DataSpace *dp_cal, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float **vcount_float, int myid, int numprocs);

void KRayShot(ModelSpace *mp, DataSpace *dp_cal, Ray *rayA, Ray *rayB, int i, int j, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);

void initialize_velocity(bspline_para bpx, bspline_para bpz, basisfun_para basis, float **vij, int vnx, int vnz, float **vel_field, float vel_begin, float vel_end);

void initialize_model(DataSpace *dp_true, ModelSpace *mp0, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, int myid, int numprocs);

void RayPathAnalysis(Ray *ray, Geo2d geo2dv, int **vcount);

#endif
