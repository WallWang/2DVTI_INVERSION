#include "su.h"
#include "stereo_struct.h"
#include "data_fitting_err.h"

void cal_cost_fun(float *fit_error, float *misfit, DataSpace *dp_true, DataSpace *dp_cal)
{//compute the cost function
	int ndata, i;
	double error;

	error= 0;
	ndata=dp_true->ndata;
	
	for(i=0;i<ndata;i++){
		error += pow((dp_cal->d[i].t   - dp_true->d[i].t)   / wt, 2);
//		error += pow((dp_cal->d[i].spx - dp_true->d[i].spx) / wp, 2);
		error += pow((dp_cal->d[i].sx  - dp_true->d[i].sx)  / wx, 2);
//		error += pow((dp_cal->d[i].rpx - dp_true->d[i].rpx) / wp, 2);
    	error += pow((dp_cal->d[i].rx  - dp_true->d[i].rx)  / wx, 2);
		misfit[i] = pow((dp_cal->d[i].t - dp_true->d[i].t) / wt, 2) + pow((dp_cal->d[i].spx - dp_true->d[i].spx) / wp, 2) + pow((dp_cal->d[i].sx - dp_true->d[i].sx) / wx, 2) + pow((dp_cal->d[i].rpx - dp_true->d[i].rpx) / wp, 2) + pow((dp_cal->d[i].rx - dp_true->d[i].rx) / wx, 2);
	}
	*fit_error=error;
}

void data_fitting_err(Data_para *dp_cal, DataSpace *dp_true, float *data_misfit, int index)
{//compute the residual 

	
	data_misfit[0] = (dp_true->d[index].t - dp_cal->t) / wt;
//	data_misfit[1] = (dp_true->d[index].spx - dp_cal->spx) / wp;
	data_misfit[2] = (dp_true->d[index].sx - dp_cal->sx) / wx;
//	data_misfit[3] = (dp_true->d[index].rpx - dp_cal->rpx) / wp;
	data_misfit[4] = (dp_true->d[index].rx - dp_cal->rx) / wx;
}
