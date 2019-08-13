#pragma once
#ifndef UPDATE_MODEL_H_
#define UPDATE_MODEL_H_
void fedcpy(BGfield **vel_tmp, BGfield **vel, int nx, int nz);
/* update model while cost function decrease */
void update_model(ModelSpace *mp0, ModelSpace *mp1, DataSpace *dp_cal, DataSpace *dp_true, bspline_para bpx, bspline_para bpz, basisfun_para basis, float *grad, float *dstep, float *misfit,float *error, int iter, int ndip_iter, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta, float vel_min, float vel_max, float epos_min, float epos_max, float delta_min, float delta_max, float **vcount, int myid, int numprocs);
#endif // !UPDATE_MODEL_H_
