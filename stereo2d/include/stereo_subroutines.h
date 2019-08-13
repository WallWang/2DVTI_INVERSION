#ifndef _STEREO_SUBROUTINES_H
#define _STEREO_SUBROUTINES_H

void input_data(char *parfile,int datatype,DataSpace *dp_true);

void input_model(char *modelfile,int datatype,ModelSpace *mp0);

void check_data(char *parfile,char *filter_data_filename,char *modelfile,char *filter_model_filename,int datatype,DataSpace *dp_true,Geo2d geo2dv,BGfield **vel,BGfield **epos,BGfield **delta);

void write_data(char *parfile,int datatype,DataSpace *dp_true);

void write_model(char *modelfile, int datatype, ModelSpace *mp);

void dip_bar(ModelSpace *mp,char *filename,int myid,int numprocs);

void error_report(ModelSpace *mp0, ModelSpace *mp1, DataSpace *dp_true, DataSpace *dp_cal, int *ndata_globe, int *Nmodel, Geo2d geo2dv, bspline_para bpx, bspline_para bpz, float *data_misfit, float *misfit, float max_dip, int verbose, int myid, int numprocs);

void delete_wrong_data(ModelSpace *mp0, ModelSpace *mp1, DataSpace *dp_true, DataSpace *dp_cal, int *ndata_globe, int index, float *data_misfit, float *misfit, int *Nmodel, bspline_para bpx, bspline_para bpz);

void seabed_parse(char *seabed_filename,float fx,float fz,bspline_para bpx,bspline_para bpz,float **seabed,int *seabed_zdepth,int verbose,int myid,int numprocs);

void Vel_Statistics(BGfield **vel_field,int nx,int nz,float vel_min,float vel_max,int *vel_min_sum,int *vel_max_sum);

void check_data_ini(char *parfile, char *fileter_data_filename, int datatype, DataSpace *dp_true, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta);
#endif
