#pragma once
#ifndef RAYSUB_H_
#define RAYSUB_H_
void draw_raypath(FILE *fp_ray_path, Ray ray); //record raypath coordinate

void dv2(Geo2d geo2dv, BGfield **vel); // calculate second derivatives

void vel2Interp(Geo2d geo2dv, BGfield **vel, float x, float z, double *u, double *ux, double *uz, double *uxx, double *uxz, double *uzz); //calculate derivatives in grid

void dynamic_ray_tracing(Ray *ray, Source sc, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float timeout);

void rkt1_for_dynamic(double *y, int n, double *z, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);//classical Runge-Kutta to calculate dynamic ray tracing

void func_for_dynamic(double *y, double *d, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);

void kinematic_ray_tracing(Ray *ray, Source sc, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float timeout);

void rkt1_for_kinematic(double *y, int n, double *z, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);//classical Runge-Kutta to calculate kinematic raytracing

void func_for_kinematic(double *y, double *d, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);
#endif // !RAYSUB_H_
