#include "su.h"
#include "sparse_matrix.h"

void Matrix_Add(tsmatrix *F, int row, int col, float value)
{
	F->data[F->tu].e= value;
	F->data[F->tu].i= row;
	F->data[F->tu].j= col;
	F->tu+=1;	
}

void Matrix_Print(tsmatrix *F, char *file)
{
	FILE *fp;
	int i;
	fp = efopen(file, "wt");
	for(i=0; i<F->tu; i++) {
		fprintf(fp, "%d %d %e\n", F->data[i].i, F->data[i].j, F->data[i].e);		
	}
	efclose(fp);	
}

void cal_F1_F2_F3(tsmatrix *F, float **F1, int **F2, int *F3 )
{
	int i, j, k=0, rowf, columf, new_count;
	float scale_colum, scale_row, vscale0,tem1,tem2;
	int *temi=alloc1int(F->mu);

	for (i=0; i<F->mu; i++)
		temi[i]=0;
		
	//order F include 2 steps:\
	//step 1:
	for(k=0;k<F->tu;k++)
	{
		rowf=F->data[k].i;	
		columf=F->data[k].j;

		F2[rowf][ temi[rowf] ]=columf;
		F1[rowf][ temi[rowf] ]=F->data[k].e;
		temi[rowf]=temi[rowf]+1;
	}

	//calculate F1,F2
	new_count=-1;
	for (i=0; i<F->mu; i++)
	{
		for(j=0; j<F3[i]; j++)
		{	new_count++;
			for(k=j+1;k<F3[i];k++)
			{
		if(F2[i][j]>F2[i][k]){	//printf("pass step2\n");
			tem2=F2[i][j];          tem1=F1[i][j];
			F2[i][j]=F2[i][k];      F1[i][j]=F1[i][k];
			F2[i][k]=tem2;          F1[i][k]=tem1;  }
			}
			//step 2 ..over

		F->data[new_count].i=i;
		F->data[new_count].j=F2[i][j];
		F->data[new_count].e=F1[i][j];

		}
	//if(j==0)	printf("row=%-3d,  j==0,  e=%-e",i,F1[i][j]);
	}
	
	free1int(temi);

	return;
}
