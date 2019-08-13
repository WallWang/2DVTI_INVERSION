#pragma once
#ifndef DATA_FITTING_ERR_H_
#define DATA_FITTING_ERR_H_
//calculate cost function
void cal_cost_fun(float *fit_error, float *misfit, DataSpace *dp_true, DataSpace *dp_cal);
//calculate data_fitting error
void data_fitting_err(Data_para *dp_cal, DataSpace *dp_true, float *data_misfit, int index);
#endif // !DATA_FITTING_ERR_H_
