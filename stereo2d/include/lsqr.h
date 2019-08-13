#ifndef _LSQR_H
#define _LSQR_H

void avpu(int m,int n,float *u,float *v,float **a,int **ax,int *axl);
void atupv(int m,int n,float *u,float *v,float **a,int **ax,int *axl);
void normlz(int n,float *x,float *s);
void lsqr(int m,int n,float *x,float  *u,int itmax,float **a,int **ax,int *axl,int rank);

#endif
