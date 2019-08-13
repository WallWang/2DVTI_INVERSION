#ifndef _IOFILE_H
#define _IOFILE_H

void write1(void *data, int size, int N, char *name);
void read1(void *data, int size, int Nx, char *name);

void write1float(float *data, int N, char *name);
void write1float_append(float *data, int N, char *name);
void write2float(float **data, int Nx, int Nz, char *name);
void write3float(float ***data, int Nx, int Ny, int Nz, char *name);
void write2int(int **data, int Nx, int Nz, char *name);
void write1txt(float *data, int N, char *name);
void write1double(double *data, int N, char *name);
void write2double(double **data, int Nx, int Nz, char *name);
void read1float(float *data, int Nx, char *name);
void read2float(float **data, int Nx, int Nz, char *name);
void read3float(float ***data, int Nx, int Ny, int Nz, char *name);
void write2field(BGfield **vel,int Nx,int Ny,char *name);

#endif

