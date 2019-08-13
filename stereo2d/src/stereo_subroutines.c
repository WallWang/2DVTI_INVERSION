#include "mpi.h"
#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include "RaySub.h"
#include "bsplines.h"
#include "sparse_matrix.h"
#include "stereo_subroutines.h"
#include "forward_modeling.h"

void check_data(char *parfile, char *filter_data_filename, char *modelfile, char *filter_model_filename, int datatype, DataSpace *dp_true, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta) 
{
	FILE *parfp;
	FILE *mfp;
	int i, nLine;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;//wx
	double delta0, deldx, deldz, delxx, delxz, delzz;//wx
	float sx, spx, rx, rpx, t;
	float xcor, zcor, theta_s, theta_r;
	int iCount;
	DataSpace dp_tmp;
	ModelSpace mp_tmp;
	
	parfp = efopen(parfile, "r");
	mfp = efopen(modelfile, "r");
	if (datatype)
		nLine = nStructCount(parfp, FSIZE * 5);
	else
		nLine = nLineCount(parfp);

	dp_tmp.d = (Data_para*)ealloc1(nLine, sizeof(Data_para));
	mp_tmp.m = (Model_para*)ealloc1(nLine, sizeof(Model_para));

	iCount = 0;
	for (i = 0; i < nLine; i++) {	
		if (datatype) {
			efread(&sx, FSIZE, 1, parfp);
			efread(&spx, FSIZE, 1, parfp);
			efread(&rx, FSIZE, 1, parfp);
			efread(&rpx, FSIZE, 1, parfp);
			efread(&t, FSIZE, 1, parfp);
			efread(&xcor, FSIZE, 1, mfp);
			efread(&zcor, FSIZE, 1, mfp);
			efread(&theta_s, FSIZE, 1, mfp);
			efread(&theta_r, FSIZE, 1, mfp);
		}
		else {
			fscanf(parfp, "%f%f%f%f%f", &sx, &spx, &rx, &rpx, &t);
			fscanf(mfp, "%f%f%f%f", &xcor, &zcor, &theta_s, &theta_r);
		}
		vel2Interp(geo2dv, vel, sx, 0, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
		vel2Interp(geo2dv, epos, sx, 0, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);//wx
		vel2Interp(geo2dv, delta, sx, 0, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);//wx
		if((1.0/vel0/vel0-(1+2*epos0)*spx*spx)/(1-2*(epos0-delta0)*vel0*vel0*spx*spx)<0){
			printf("%d line : sx=%f spx=%f rx=%f rpx=%f t=%f is wrong, v=%f e=%f n=%f. program delete it.\n",i,sx,spx,rx,rpx,t,vel0,epos0,delta0);
			continue;
		}
		vel2Interp(geo2dv, vel, rx, 0, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
		vel2Interp(geo2dv, epos, rx, 0, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);//wx
		vel2Interp(geo2dv, delta, rx, 0, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);//wx
		if((1.0/vel0/vel0-(1+2*epos0)*rpx*rpx)/(1-2*(epos0-delta0)*vel0*vel0*rpx*rpx)<0){
			printf("%d line : sx=%f spx=%f rx=%f rpx=%f t=%f is wrong, v=%f e=%f n=%f. program delete it.\n",i,sx,spx,rx,rpx,t,vel0,epos0,delta0);
			continue;
		}
		dp_tmp.d[iCount].sx = sx;
		dp_tmp.d[iCount].spx = spx;
		dp_tmp.d[iCount].rx = rx;
		dp_tmp.d[iCount].rpx = rpx;
		dp_tmp.d[iCount].t = t;
		mp_tmp.m[iCount].xcor = xcor;
		mp_tmp.m[iCount].zcor = zcor;
		mp_tmp.m[iCount].theta_s = theta_s;
		mp_tmp.m[iCount].theta_r = theta_r;
		iCount++;
	}
	dp_true->ndata = iCount;
	efclose(parfp);
	efclose(mfp);
	if (datatype) {// for binary data
		parfp = efopen(filter_data_filename, "wb");
		mfp = efopen(filter_model_filename, "wb");
		for (i = 0; i < iCount; i++) {
			fwrite(&dp_tmp.d[i].sx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].spx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].rx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].rpx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].t, FSIZE, 1, parfp);
			fwrite(&mp_tmp.m[i].xcor, FSIZE, 1, mfp);
			fwrite(&mp_tmp.m[i].zcor, FSIZE, 1, mfp);
			fwrite(&mp_tmp.m[i].theta_s, FSIZE, 1, mfp);
			fwrite(&mp_tmp.m[i].theta_r, FSIZE, 1, mfp);
		}
	}else{
		parfp = efopen(filter_data_filename, "wt");
		mfp = efopen(filter_model_filename, "wt");
		for (i = 0; i < iCount; i++) {
			fprintf(parfp, "%0.5f\t%0.10f\t%0.5f\t%0.10f\t%0.10f\n", dp_tmp.d[i].sx, dp_tmp.d[i].spx, dp_tmp.d[i].rx, dp_tmp.d[i].rpx, dp_tmp.d[i].t);
			fprintf(mfp, "%f %f %f %f\n", mp_tmp.m[i].xcor, mp_tmp.m[i].zcor, mp_tmp.m[i].theta_s, mp_tmp.m[i].theta_r);
		}
	}
	efclose(parfp);
	efclose(mfp);
	free1(dp_tmp.d);
	free1(mp_tmp.m);
}

void input_data(char *parfile, int datatype, DataSpace *dp_true)
{
	FILE *parfp;
	int i, nLine;
	float sx, spx, rx, rpx, t;
	
	parfp = efopen(parfile, "r");
	if (datatype)
		nLine = nStructCount(parfp, FSIZE * 5);
	else
		nLine = nLineCount(parfp);

	for(i=0; i<nLine; i++)
	{	
		if(datatype) {
			efread(&sx, FSIZE, 1, parfp);
			efread(&spx, FSIZE, 1, parfp);
			efread(&rx, FSIZE, 1, parfp);
			efread(&rpx, FSIZE, 1, parfp);
			efread(&t, FSIZE, 1, parfp);
		} else 
			fscanf(parfp, "%f%f%f%f%f", &sx, &spx, &rx, &rpx, &t);
		dp_true->d[i].sx = sx;
		dp_true->d[i].spx = spx;
		dp_true->d[i].rx = rx;
		dp_true->d[i].rpx = rpx;
		dp_true->d[i].t = t;
	}
	dp_true->ndata = nLine;
	efclose(parfp);
}

void input_model(char *modelfile, int datatype, ModelSpace *mp0)
{
	FILE *mfp;
	int i, nLine;
	float xcor, zcor, theta_s, theta_r;
	
	mfp = efopen(modelfile, "r");
	if (datatype)
		nLine = nStructCount(mfp, FSIZE * 4);
	else
		nLine = nLineCount(mfp);

	for (i = 0; i < nLine; i++)
	{
		if (datatype) {
			efread(&xcor, FSIZE, 1, mfp);
			efread(&zcor, FSIZE, 1, mfp);
			efread(&theta_s, FSIZE, 1, mfp);
			efread(&theta_r, FSIZE, 1, mfp);
		}
		else
			fscanf(mfp, "%f%f%f%f", &xcor, &zcor, &theta_s, &theta_r);
		mp0->m[i].xcor = xcor;
		mp0->m[i].zcor = zcor;
		mp0->m[i].theta_s = theta_s;
		mp0->m[i].theta_r = theta_r;
	}
	mp0->ndata = nLine;
	efclose(mfp);
}

void write_data(char *parfile, int datatype, DataSpace *dp_true)
{
	FILE *parfp;
	int i;
	if(datatype) {
		parfp = efopen(parfile, "wb");
		for(i=0; i<dp_true->ndata; i++) {
			fwrite(&dp_true->d[i].sx, FSIZE, 1, parfp);
			fwrite(&dp_true->d[i].spx, FSIZE, 1, parfp);
			fwrite(&dp_true->d[i].rx, FSIZE, 1, parfp);
			fwrite(&dp_true->d[i].rpx, FSIZE, 1, parfp);
			fwrite(&dp_true->d[i].t, FSIZE, 1, parfp);
		}
	}else{
		parfp = efopen(parfile, "wt");
		for (i = 0; i < dp_true->ndata; i++)
			fprintf(parfp, "%f %f %f %f %f\n", dp_true->d[i].sx, dp_true->d[i].spx, dp_true->d[i].rx, dp_true->d[i].rpx, dp_true->d[i].t);
	}
	efclose(parfp);
}

void write_model(char *modelfile, int datatype, ModelSpace *mp)
{
	FILE *mfp;
	int i;

	if (datatype) {
		mfp = efopen(modelfile, "wb");
		for (i = 0; i < mp->ndata; i++) {
			fwrite(&mp->m[i].xcor, FSIZE, 1, mfp);
			fwrite(&mp->m[i].zcor, FSIZE, 1, mfp);
			fwrite(&mp->m[i].theta_s, FSIZE, 1, mfp);
			fwrite(&mp->m[i].theta_r, FSIZE, 1, mfp);
		}
	}
	else {
		mfp = efopen(modelfile, "wt");
		for (i = 0; i < mp->ndata; i++)
			fprintf(mfp, "%f %f %f %f\n", mp->m[i].xcor, mp->m[i].zcor, mp->m[i].theta_s, mp->m[i].theta_r);
	}
	efclose(mfp);
}

void dip_bar(ModelSpace *mp, char *filename, int myid, int numprocs)
{
	int i;
	float xcor, zcor;
	float dip, range;
	FILE *fp;
	float err;
	fp = efopen(filename, "wt");
	
	range = 50;
	for (i = 0; i < mp->ndata; i++){
		dip = angle2radian((mp->m[i].theta_s + mp->m[i].theta_r) / 2.0);
		xcor = mp->m[i].xcor + range;
		zcor = mp->m[i].zcor + range * tan(dip);
		fprintf(fp, "%5s   %18.8f   %18.8f\n", "MA", xcor, zcor);
		xcor = mp->m[i].xcor - range;
		zcor = mp->m[i].zcor - range * tan(dip);
		fprintf(fp, "%5s   %18.8f   %18.8f\n", "PA", xcor, zcor);
	}
	efclose(fp);
	return;
}

void error_report(ModelSpace *mp0, ModelSpace *mp1, DataSpace *dp_true, DataSpace *dp_cal, int *ndata_globe, int *Nmodel, Geo2d geo2dv, bspline_para bpx, bspline_para bpz, float *data_misfit, float *misfit, float max_dip, int verbose, int myid, int numprocs)
{
	int row, i;
	int index;
	float dip;
	float boundary_eps;
	float xmin, xmax, zmin, zmax;
	int cLarge_Misfit = 0;
	int cLarge_Dip = 0;
	int cOut_Model = 0;
	int cAbove_Seabed = 0;
	int cAngle = 0;

	boundary_eps = 10;
	xmin = geo2dv.xmin + boundary_eps;
	xmax = geo2dv.xmax - boundary_eps;
	zmin = geo2dv.zmin + boundary_eps;
	zmax = geo2dv.zmax - boundary_eps;
	
	for (i = 0; i < mp0->ndata; i++){
		dip = fabsf((mp0->m[i].theta_s + mp0->m[i].theta_r) / 2.0);
		row = i * 5;
		if (misfit[i] > 1e+7){
			if (myid == 0 && verbose > 1) {
				printf("error_report for large misfit, index:%d, dip:%.1f\n", i, dip);
				printf("true:");
				printf("sx=%.1f ", dp_true->d[i].sx);
				printf("ps=%f ", dp_true->d[i].spx);
				printf("rx=%.1f ", dp_true->d[i].rx);
				printf("pr=%f ", dp_true->d[i].rpx);
				printf("t=%.4f ", dp_true->d[i].t);
				printf("\n");

				printf("cal: ");
				printf("sx=%.1f ", dp_cal->d[i].sx);
				printf("ps=%f ", dp_cal->d[i].spx);
				printf("rx=%.1f ", dp_cal->d[i].rx);
				printf("pr=%f ", dp_cal->d[i].rpx);
				printf("t=%.4f ", dp_cal->d[i].t);
				printf("\n");
			
				printf("err: ");
				printf("sx=%.1f ", data_misfit[row + 2]);
				printf("ps=%f ", data_misfit[row + 1]);
				printf("rx=%.1f ", data_misfit[row + 4]);
				printf("pr=%f ", data_misfit[row + 3]);
				printf("t=%.1f ", data_misfit[row + 0]);
				printf("err=%e ", misfit[i]);
				printf("\n");
			}
			delete_wrong_data(mp0, mp1, dp_true, dp_cal, ndata_globe, i, data_misfit, misfit, Nmodel, bpx, bpz);
			cLarge_Misfit++;
			i--;
			continue;
		}// for large misfit
		
		if (dip > max_dip && misfit[i] > 1e+4){
			if (myid == 0 && verbose > 1){
				printf("error_report for large dip, index:%d, dip:%.1f\n", i, dip);
				printf("sx=%.1f ", data_misfit[row + 2]);
				printf("ps=%f ", data_misfit[row + 1]);
				printf("rx=%.1f ", data_misfit[row + 4]);
				printf("pr=%f ", data_misfit[row + 3]);
				printf("t=%.1f ", data_misfit[row + 0]);
				printf("err=%e ", misfit[i]);
				printf("\n");
			}
			delete_wrong_data(mp0, mp1, dp_true, dp_cal, ndata_globe, i, data_misfit, misfit, Nmodel, bpx, bpz);
			cLarge_Dip++;
			i--;
			continue;
		}// for large dip
		
		if (mp0->m[i].xcor<xmin || mp0->m[i].xcor>xmax || mp0->m[i].zcor<zmin || mp0->m[i].zcor>zmax) {
			if (myid == 0 && verbose > 1) {
				printf("error_report for out of model, index:%d, dip:%.1f\n", i, dip);
				printf("mp.x=%.1f ", mp0->m[i].xcor);
				printf("mp.z=%.1f ", mp0->m[i].zcor);
				printf("xmin=%.1f ", xmin);
				printf("xmin=%.1f ", xmax);
				printf("\n");
			}
			delete_wrong_data(mp0, mp1, dp_true, dp_cal, ndata_globe, i, data_misfit, misfit, Nmodel, bpx, bpz);
			cOut_Model++;
			i--;
			continue;
		}// for out of model

		if (ABS(mp0->m[i].theta_s - mp0->m[i].theta_r) > 200) {
			if (myid == 0 && verbose > 1) {
				printf("error_report for angle, index:%d, dip:%.1f\n", i, dip);
				printf("mp.x=%.1f ", mp0->m[i].xcor);
				printf("mp.z=%.1f ", mp0->m[i].zcor);
				printf("mp.theta_s=%.1f ", mp0->m[i].theta_s);
				printf("mp.theta_r=%.1f ", mp0->m[i].theta_r);
				printf("\n");
			}
			delete_wrong_data(mp0, mp1, dp_true, dp_cal, ndata_globe, i, data_misfit, misfit, Nmodel, bpx, bpz);
			cAngle++;
			i--;
			continue;
		}// for angle
	}
	
	if(myid==0 && verbose) {
		if(cLarge_Misfit)	printf("Large misfit del=%d\n", cLarge_Misfit);
		if(cLarge_Dip)		printf("Large Dip    del=%d\n", cLarge_Dip);
		if(cOut_Model)		printf("Out of Model del=%d\n", cOut_Model);
		if(cAngle)		    printf("Angle 	     del=%d\n", cAngle);
	}
}

void delete_wrong_data(ModelSpace *mp0, ModelSpace *mp1, DataSpace *dp_true, DataSpace *dp_cal, int *ndata_globe, int index, float *data_misfit, float *misfit, int *Nmodel, bspline_para bpx, bspline_para bpz)
{
	int i, ndata;
	int row1, row2;
	ndata = *ndata_globe;
	for (i = index + 1; i < ndata; i++) {
		dp_true->d[i - 1] = dp_true->d[i];
		dp_cal->d[i - 1] = dp_cal->d[i];
		mp0->m[i - 1] = mp0->m[i];
		mp1->m[i - 1] = mp1->m[i];
		misfit[i - 1] = misfit[i];
		row1 = (i - 1) * 5;
		row2 = i * 5;
		data_misfit[row1 + 0] = data_misfit[row2 + 0];
		data_misfit[row1 + 1] = data_misfit[row2 + 1];
		data_misfit[row1 + 2] = data_misfit[row2 + 2];
		data_misfit[row1 + 3] = data_misfit[row2 + 3];
		data_misfit[row1 + 4] = data_misfit[row2 + 4];
	}
	ndata = ndata - 1;
	mp0->ndata = ndata;
	mp1->ndata = ndata;
	dp_true->ndata = ndata;
	dp_cal->ndata = ndata;
	*ndata_globe = ndata;
	*Nmodel = 3 * bpx.Norigin * bpz.Norigin;
	return;
}

void seabed_parse(char *seabed_filename, float fx, float fz, bspline_para bpx, bspline_para bpz, float **seabed, int *seabed_zdepth, int verbose, int myid, int numprocs)
{
	FILE *fp;
	int nLine;
	float center;
	int index, i;
	
	fp = efopen(seabed_filename, "r");
	nLine = nLineCount(fp);
	for (i = 0; i < nLine; i++) {
		//1 for z, 0 for z
		fscanf(fp, "%f %f", &seabed[0][i], &seabed[1][i]);
	}
	efclose(fp);

	for (i = 0; i < bpx.Norigin; i++) {
		center = fx + (i + 2) * bpx.interval;
		index = binary_find(&seabed[1][0], 0, nLine - 1, center);
		seabed_zdepth[i] = (seabed[0][index] - fz) / bpz.interval - bpz.K;

		if (seabed_zdepth[i] < 0) seabed_zdepth[i] = 0;
		if (seabed_zdepth[i] >= bpz.Norigin) seabed_zdepth[i] = bpz.Norigin - 1;
		if (verbose > 1 && myid == 0) printf("i=%d index_of_center=%d c=%f x=%f z=%f iz=%d\n", i, index, center, seabed[1][index], seabed[0][index], seabed_zdepth[i]);
	}
}

void Vel_Statistics(BGfield **vel_field, int nx, int nz, float vel_min, float vel_max, int *vel_min_sum, int *vel_max_sum)
{
	int i, j;
	*vel_min_sum = 0;
	*vel_max_sum = 0;
	for (i = 0; i < nx; i++) 
		for (j = 0; j < nz; j++) {
			if (vel_field[i][j].u < vel_min)
				(*vel_min_sum)++;
			else if (vel_field[i][j].u > vel_max)
				(*vel_max_sum)++;
		}
}

void check_data_ini(char *parfile, char *filter_data_filename, int datatype, DataSpace *dp_true, Geo2d geo2dv, BGfield **vel, BGfield **epos, BGfield **delta)
{
	FILE *parfp;
	int i, nLine;
	double vel0, dvdx, dvdz, uxx, uxz, uzz;
	double epos0, dedx, dedz, epxx, epxz, epzz;
	double delta0, deldx, deldz, delxx, delxz, delzz;
	float sx, spx, rx, rpx, t;
	int iCount;
	DataSpace dp_tmp;

	parfp = efopen(parfile, "r");
	if (datatype)
		nLine = nStructCount(parfp, FSIZE * 5);
	else
		nLine = nLineCount(parfp);
	dp_tmp.d = (Data_para*)ealloc1(nLine, sizeof(Data_para));

	iCount = 0;
	for (i = 0; i < nLine; i++)
	{
		if (datatype) {
			efread(&sx, FSIZE, 1, parfp);
			efread(&spx, FSIZE, 1, parfp);
			efread(&rx, FSIZE, 1, parfp);
			efread(&rpx, FSIZE, 1, parfp);
			efread(&t, FSIZE, 1, parfp);
		}
		else 
			fscanf(parfp, "%f%f%f%f%f", &sx, &spx, &rx, &rpx, &t);
		vel2Interp(geo2dv, vel, sx, 0, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
		vel2Interp(geo2dv, epos, sx, 0, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);//wx
		vel2Interp(geo2dv, delta, sx, 0, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);//wx
		if ((1.0/vel0/vel0-(1+2*epos0)*spx*spx)/(1-2*(epos0-delta0)*vel0*vel0*spx*spx)<0) {
			printf("%d line : sx=%f spx=%f rx=%f rpx=%f t=%f is wrong, v=%f e=%f n=%f. program delete it.\n", i, sx, spx, rx, rpx, t, vel0, epos0, delta0);
			continue;
		}
		vel2Interp(geo2dv, vel, rx, 0, &vel0, &dvdx, &dvdz, &uxx, &uxz, &uzz);
		vel2Interp(geo2dv, epos, rx, 0, &epos0, &dedx, &dedz, &epxx, &epxz, &epzz);//wx
		vel2Interp(geo2dv, delta, rx, 0, &delta0, &deldx, &deldz, &delxx, &delxz, &delzz);//wx
		if ((1.0/vel0/vel0-(1+2*epos0)*rpx*rpx)/(1-2*(epos0-delta0)*vel0*vel0*rpx*rpx)<0) {
			printf("%d line : sx=%f spx=%f rx=%f rpx=%f t=%f is wrong, v=%f e=%f n=%f. program delete it.\n", i, sx, spx, rx, rpx, t, vel0, epos0, delta0);
			continue;
		}
		dp_tmp.d[iCount].sx = sx;
		dp_tmp.d[iCount].spx = spx;
		dp_tmp.d[iCount].rx = rx;
		dp_tmp.d[iCount].rpx = rpx;
		dp_tmp.d[iCount].t = t;
		iCount++;
	}
	dp_true->ndata = iCount;
	efclose(parfp);
	
	if (datatype) {//for binary data
		parfp = efopen(filter_data_filename, "wb");
		for (i = 0; i < iCount; i++) {
			fwrite(&dp_tmp.d[i].sx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].spx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].rx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].rpx, FSIZE, 1, parfp);
			fwrite(&dp_tmp.d[i].t, FSIZE, 1, parfp);
		}
	}
	else {//for text data
		parfp = efopen(filter_data_filename, "wt");
		for (i = 0; i < iCount; i++)
		  fprintf(parfp, "%f %f %f %f %f\n", dp_tmp.d[i].sx, dp_tmp.d[i].spx, dp_tmp.d[i].rx, dp_tmp.d[i].rpx, dp_tmp.d[i].t);
	}
	efclose(parfp);
	free1(dp_tmp.d);
}
