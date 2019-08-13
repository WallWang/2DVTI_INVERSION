#ifndef _WX_H
#define _WX_H

typedef struct PointTag
{
	int lX;
	int lY;
}Point;

double tic (void);
double toc (double t);//relate to record running time

int nLineCount(FILE *infile);
int nLineCount_filename(char *filename);
int nFileSize(FILE *infile);
int nFileSize_filename(char *filename);
int nStructCount(FILE *infile, int StructSize);//relate to file input

float radian2angle(float radian);
float angle2radian(float angle);//degrees and radians transform

void mpi_gatherv_construct(int start, int end, int numprocs, int ntask, int *recvcounts, int *displs);

void normalization(float *trace, int lt);

void float1max(float *d, int n, float *max, int *pos);
void float1min(float *d, int n, float *min, int *pos);

void zero1float(float *p, int n1);
void zero2float(float **p, int n1, int n2);
void zero2int(int **p, int n1, int n2);//init matrix

float average(float *data, int n);
float average_abs(float *data, int n);

int binary_find(float * arr, int low, int high, float key);
void fieldcpy(BGfield **vel_tmp,BGfield **vel,int n1,int n2);
void field_assign(float **vel_field,BGfield **vel,int nx,int nz);

#endif

