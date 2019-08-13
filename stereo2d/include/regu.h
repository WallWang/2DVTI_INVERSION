#pragma once
#ifndef REGU_H_
#define REGU_H_
void cal_damp_regu(float **Hess, int ndata, bspline_para bpx, bspline_para bpz, float regu_ray_damp, float regu_damp, int iter, int ndip_iter, float **regu_vij);
void cal_homo_regu(bspline_para bpx, bspline_para bpz, float regu_homox, float regu_homoz, float **Hess, int ndata);
void vij_statistic(bspline_para bpx, bspline_para bpz, float **regu_vij, float **vcount);
void dip_statistic(ModelSpace *mp0, Geo2d geo2dv, bspline_para bpx, bspline_para bpz, float **dip_sum, float **dip_count);
void dip_statistic_small(ModelSpace *mp0, Geo2d geo2dv, bspline_para bpx, bspline_para bpz, float **dip_sum, float **dip_count);
int none_zero_count(float **data, int nx, int nz);
#endif // !REGU_H_
