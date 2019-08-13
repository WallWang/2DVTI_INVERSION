#pragma once
#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

/*struct triple--save sparse matrix in compact way
i,j---row and cloumn position in sparse matrix 
e---value of unzero element in sparse matrix*/
typedef struct {
	int i, j;
	float e;
}triple;

/*rpos---the first 10 element's position chart in every row 
mu,nu,tu---row,column and !0 element's number
*num---the number of !0 elements in every row*/
typedef struct {
	triple *data;
	int mu, nu, tu;
}tsmatrix;

void cal_F1_F2_F3(tsmatrix *F, float **F1, int **F2, int *F3);
void Matrix_Add(tsmatrix *F, int row, int col, float value);
void Matrix_Print(tsmatrix *F, char *file);

#endif // !SPARSE_MATRIX_H_


