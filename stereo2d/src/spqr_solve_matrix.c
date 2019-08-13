#include "su.h"
#include "SuiteSparseQR_C.h"
#define Long SuiteSparse_long
#include "stereo_struct.h"
#include "wx.h"
#include "sparse_matrix.h"

void spqr_solve_matrix(float *model_update, float *data_misfit, tsmatrix *F)
{
	int i;
        cholmod_sparse *A;
        cholmod_triplet T;
        cholmod_dense b;
        cholmod_dense *x;
        cholmod_common Common, *cc ;

        double rnorm, one [2] = {1,0}, minusone [2] = {-1,0} ;
        cholmod_dense *Residual;

        cc = &Common ; /* initialize cholmod */
        cholmod_l_start (cc) ;

	double t;
	t = tic () ;
        T.nzmax=F->tu;
        T.nrow=F->mu;
        T.ncol=F->nu;
        T.nnz=T.nzmax;
        T.i=(SuiteSparse_long *)malloc(sizeof(SuiteSparse_long)*T.nzmax);
        if(T.i==NULL) { printf("T.i memory allocation failed\n"); exit(1); }
        T.j=(SuiteSparse_long *)malloc(sizeof(SuiteSparse_long)*T.nzmax);
        if(T.j==NULL) { printf("T.j memory allocation failed\n"); exit(1); }
        T.x=(double *)malloc(sizeof(double)*T.nzmax);
        if(T.x==NULL) { printf("T.x memory allocation failed\n"); exit(1); }
        T.z=NULL;
        T.stype=0; /* matrix is "unsymmetric" */
        T.itype=CHOLMOD_LONG;
        T.xtype=CHOLMOD_REAL;
        T.dtype=CHOLMOD_DOUBLE;

        b.nrow=T.nrow;
        b.ncol=1;
        b.nzmax=T.nrow;
        b.d=T.nrow; /* leading dimension (d >= nrow must hold) */
        b.x=(double *)malloc(sizeof(double)*T.nrow);
        if(b.x==NULL) { printf("b.x memory allocation failed\n"); exit(1); }
        b.z=NULL;
        b.xtype=CHOLMOD_REAL; /* pattern, real, complex, or zomplex */
        b.dtype=CHOLMOD_DOUBLE; /* x and z double or float */

        printf("T: nrow=%6ld ncol=%6ld\nnzmax=%6ld\n", T.nrow, T.ncol, T.nzmax);
        printf("Valid data percentage %.2lf\n", ((double)T.nzmax)/T.nrow/T.ncol*100.0);
//	printf("itype=%d xtype=%d dtype=%d\n", T.itype, T.xtype, T.dtype);
	//printf("b: nrow=%6ld ncol=%6ld\nnzmax=%6ld\nxtype=%d dtype=%d\n", b.nrow, b.ncol, b.nzmax, b.xtype, b.dtype);

        for (i=0; i<T.nzmax; i++){
		if(F->data[i].i<0 || F->data[i].i >=F->mu || F->data[i].j<0 ||F->data[i].j>=F->nu) 
			printf("F's index out of range: %d %d   mu=%d,nu=%d\n", F->data[i].i, F->data[i].j,F->mu,F->nu);
		
		if(isnan(F->data[i].e)) printf("F[%d][%d]=nan\n", F->data[i].i, F->data[i].j);

                *((SuiteSparse_long *)T.i+i)=F->data[i].i;
                *((SuiteSparse_long *)T.j+i)=F->data[i].j;
                *((double *)T.x+i)=F->data[i].e;
        }

//	for (i=0; i<50; i++){
//		printf("F[%d][%d]=%f\n", F->data[i].i, F->data[i].j, F->data[i].e);
//	}
//	for (i=0; i<50; i++){
//		printf("data_misfit[%d]=%f\n", i, data_misfit[i]);
//	}

        for(i=0; i<T.nrow; i++) {
		*((double *)b.x+i) = data_misfit[i];
		if(isnan(data_misfit[i])) printf("data_misfit[%d]=nan\n", i);
	}

        A = cholmod_l_triplet_to_sparse(&T, 0, cc) ;
        x = SuiteSparseQR_C_backslash_default (A, &b, cc) ;

        for(i=0; i<T.ncol; i++) model_update[i]= *((double *)x->x+i);

	//printf("x: nrow=%6ld ncol=%6ld\nnzmax=%6ld\nxtype=%d dtype=%d\n", x->nrow, x->ncol, x->nzmax, x->xtype, x->dtype);

        /* rnorm = norm (B-A*X) */
        Residual = cholmod_l_copy_dense (&b, cc) ;
        cholmod_l_sdmult (A, 0, minusone, one, x, Residual, cc) ;
        rnorm = cholmod_l_norm_dense (Residual, 2, cc) ;
        printf ("2-norm of residual: %8.1e\n", rnorm) ;
        printf ("number of TBB tasks %ld\n", cc->SPQR_istat [3]) ;
        //printf ("rank %ld\n", cc->SPQR_istat [4]) ;

        free(T.i);
        free(T.j);
        free(T.x);
        free(b.x);

        cholmod_l_free_sparse (&A, cc) ;
        cholmod_l_free_dense (&x, cc) ;
        cholmod_l_free_dense (&Residual, cc);

	printf ("time: %8.2f s\n", toc (t)) ;

        cholmod_l_finish (cc) ;
}
