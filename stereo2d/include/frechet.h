#pragma once
#ifndef FRECHET_H_
#define FRECHET_H_
void cal_Frechet_matrix(float *pf, float **Hess, bspline_para bpx, bspline_para bpz, basisfun_para basis, ModelSpace *mp0, DataSpace *dp_true, int iter, int ndip_iter, int Nmodel, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, int myid, int numprocs);
void Frechet_matrix_x_z_sita(DataSpace *dp_cal, tsmatrix *F, int index, Ray *rayA, Ray *rayB, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);
void Frechet_matrix_field(float **F, bspline_para bpx, bspline_para bpz, basisfun_para basis, ModelSpace *mp0, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, Ray *rayA, Ray *rayB);
void RayPathFieldFrechet(Ray *ray, basisfun_para basis, int nodex, int nodez, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, DField *f);
void RayPathZone(Geo2d geo2dv, bspline_para bpx, bspline_para bpz, Ray *ray, int **vij);
void cal_Newton_gradient(float *b,float **A,float *x,int Nmodel);
#endif // !FRECHET_H_
