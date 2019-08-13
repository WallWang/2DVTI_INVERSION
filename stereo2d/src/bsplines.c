#include "su.h"
#include "bsplines.h"
#include "sparse_matrix.h"
#include "spqr_solve_matrix.h"
#include "lsqr.h"

void bparameter_initial(bspline_para  *bp, float DX, float dx, int nx)
{
	int nrange;
	bp->K=3;
	bp->interval=DX;
	bp->dsample=dx;
	
	nrange=(nx-1)*dx;
	bp->Nnode=nrange/bp->interval +1;
	if(nrange%bp->interval!=0)
		bp->Nnode+=1;
	bp->Norigin=bp->Nnode -1 - bp->K;
	bp->Nfun=nx;
}

void node_initial( bspline_para *bp )
{
	int i;
	for(i=1; i< bp->Nnode-1; i++)
	{
		bp->node[i]=i* bp->interval;
	}
	bp->node[0]=-1.0e+8;
	//bp->node[0]=0;
	bp->node[bp->Nnode-1]=1.0e+8;	/* deal with right boudary */
}

void basis_spline(bspline_para  bp, float **bspline, float **dbspline, float **ddbspline)
{
	float *b, *db, *ddb, x=0.0;
	b=  (float*)calloc(bp.Nnode, sizeof(float));
	db= (float*)calloc(bp.Nnode, sizeof(float));
	ddb=(float*)calloc(bp.Nnode, sizeof(float));
	int i, j;	
	for(j=0; j<bp.Nfun; j++)
	{
		x=j*bp.dsample;//add
		basis_analytical_value (bp, b, db, ddb, x);
		for(i=0; i<bp.Nnode; i++)
		{
			bspline[i][j]=b[i];
			dbspline[i][j]=db[i];
			ddbspline[i][j]= ddb[i];
		}
	}
	free1float(b);
	free1float(db);
	free1float(ddb);
}

void basis_analytical_value (bspline_para  bp, float *bspline, float *dbspline, float *ddbspline, float x)
{
	int i, m;
	float coe1, coe2, dcoe1, dcoe2;
	float x0, bx0, dbx0;
	if ( x<0 )
		x= -x;
	/*------------deal with left boundary----------------*/
	x0= (bp.Nfun-1)*bp.dsample- x;	
	/*--initial bpline for order zero----*/
	for(i=0; i<bp.Nnode-1; i++)
	{
		if(x0>=bp.node[i] && x0<bp.node[i+1])
			bspline[i]=1.0;
		else
			bspline[i]=0.0;
	}
	/*---calculate bspline and its derivatives of order m----*/
	for(m=1; m<=bp.K; m++)
	{
		for(i=0; i<bp.Nnode-m-1; i++)
		{
			coe1=(x0-bp.node[i])/ (bp.node[i+m]- bp.node[i]);
			coe2=(bp.node[i+m+1]-x0)/ (bp.node[i+m+1]- bp.node[i+1]);
			dcoe1=bp.node[i+m]- bp.node[i];
			dcoe2=bp.node[i+m+1]- bp.node[i+1];
			if (m==bp.K-1 )
				dbspline[i]=m*(bspline[i]/dcoe1- bspline[i+1]/dcoe2);
			else if(m==bp.K)
			{
				ddbspline[i]=m*(dbspline[i]/dcoe1- dbspline[i+1]/dcoe2);
				dbspline[i]=m*(bspline[i]/dcoe1- bspline[i+1]/dcoe2);
			}
			bspline[i]=coe1*bspline[i]+ coe2*bspline[i+1];
		}	
	}
	bx0= bspline[bp.Norigin- 1];
	dbx0= -dbspline[bp.Norigin- 1];

	/*--initial bpline for order zero----*/
	for(i=0; i<bp.Nnode-1; i++)
	{
		if(x>=bp.node[i] && x<bp.node[i+1])
			bspline[i]=1.0;
		else
			bspline[i]=0.0;
	}
	/*---calculate bspline and its derivatives of order m----*/
	for(m=1; m<=bp.K; m++)
	{
		for(i=0; i<bp.Nnode-m-1; i++)
		{
			coe1=(x-bp.node[i])/ (bp.node[i+m]- bp.node[i]);
			coe2=(bp.node[i+m+1]-x)/ (bp.node[i+m+1]- bp.node[i+1]);
			dcoe1=bp.node[i+m]- bp.node[i];
			dcoe2=bp.node[i+m+1]- bp.node[i+1];
			if (m==bp.K-1 )
				dbspline[i]=m*(bspline[i]/dcoe1- bspline[i+1]/dcoe2);
			else if(m==bp.K)
			{
				ddbspline[i]=m*(dbspline[i]/dcoe1- dbspline[i+1]/dcoe2);
				dbspline[i] =m*(bspline[i]/dcoe1- bspline[i+1]/dcoe2);
			}
			bspline[i]=coe1*bspline[i]+ coe2*bspline[i+1];
		}	
	}
	
	//bspline[0]=  bx0;
	//dbspline[0]= dbx0;
}

void cal_vel_field(bspline_para bpx, bspline_para bpz, basisfun_para basis, float **vij, float **vel_out)
/*
	vij: the bspline velocity model(input) 
	vel_out:the whole velocity model(output)
 */
{
	int i, j, k, l;
	int ks,ls;

	zero2float(vel_out, bpz.Nfun, bpx.Nfun);

	for(i=0; i<bpx.Nfun; i++) {
		for(j=0; j<bpz.Nfun; j++) {

			ks=i*bpx.dsample/bpx.interval;
			ls=j*bpz.dsample/bpz.interval;

			if(ks<3) ks=3;
			if(ls<3) ls=3;
			if(ls>bpz.Norigin-1) ls=bpz.Norigin-1;
			if(ks>bpx.Norigin-1) ks=bpx.Norigin-1;
			
			for(k=ks; k>=ks-3; k--) {
				for(l=ls; l>=ls-3; l--) { 
			//for(k=0; k<bpx.Norigin; k++) {
			//	for(l=0; l<bpz.Norigin; l++) { 
					vel_out[i][j] += vij[k][l] * basis.bsplinex[k][i] * basis.bsplinez[l][j];
				}
			}			
			
		}
	}
}

void bspline_interpolation(float **vel, bspline_para bpx,  bspline_para bpz, float **vij, basisfun_para basis)
/*
	vij: the bspline velocity model(output) 
*/
{
	float *m, *d;
	int i, j, k, l, scr_r, scr_c;
	float bx, bz;
	int ks,ls;
	int nvx,nvz,nv,nvijx,nvijz,nvij;

	nvx=bpx.Nfun;
	nvz=bpz.Nfun; 
	nv =nvx*nvz;
	nvijx= bpx.Norigin;
	nvijz= bpz.Norigin;
	nvij = nvijx*nvijz; 

	m= ealloc1float(nvij);
	d= ealloc1float(nv);
	zero1float(m, nvij);
	zero1float(d, nv);

	tsmatrix F;
	F.data = (triple *)( ealloc1( nv*16, sizeof(triple) ) );
	F.mu=nv;
	F.nu=nvij;
	F.tu=0;

	//calculte Frechet derivatives of interpolated velocity
	for(i=0; i<bpx.Nfun; i++) {
		for(j=0; j<bpz.Nfun; j++) {
			ks=i*bpx.dsample/bpx.interval;
			ls=j*bpz.dsample/bpz.interval;

			if(ks<3) ks=3;
			if(ls<3) ls=3;
			if(ls>bpz.Norigin-1) ls=bpz.Norigin-1;
			if(ks>bpx.Norigin-1) ks=bpx.Norigin-1;
			for(k=ks; k>=ks-3; k--) {
				for(l=ls; l>=ls-3; l--) {
					bx = basis.bsplinex[k][i];
					bz = basis.bsplinez[l][j];
					scr_r = i*nvz+ j;
					scr_c = k*nvijz+ l;
					Matrix_Add(&F, scr_r, scr_c, bx*bz);
				}
			}
		}	
	}

	//calculate right hand iterm of interpolation equation
	for( i=0; i<bpx.Nfun; i++ ) {
		for( j=0; j<bpz.Nfun; j++ ) {
			scr_r= i*nvz+ j;
			d[scr_r]= 0.0- vel[i][j]; 
		}
	}


        float **F1; 
        int **F2; 
        int *F3; 

	F1=(float**)malloc(F.mu* sizeof(float*) );
	F2=(int**)malloc(F.mu* sizeof(int*) );
	F3=ealloc1int(F.mu);
	//calculate F3 and allocate memory for F1,F2
	for (i=0; i<F.mu; i++) {
		F3[i]= 16;
		F1[i]= ealloc1float( F3[i] );
		F2[i]= ealloc1int( F3[i] );
	}

	//calculate interpolation coes vij
	if(0) {
		spqr_solve_matrix(m, d, &F);
	} else {
        	cal_F1_F2_F3(&F, F1, F2, F3);
		lsqr( F.mu, F.nu, m, d, 401, F1, F2, F3, 1);//calculate the solution
	}

	//free memory
        for (i=0; i<F.mu; i++) {
              free1float(F1[i]);
              free1int(F2[i]);
        }
	free(F1); free(F2); 
	free1int(F3);

	//save interpolation coes
	for(i=0; i<bpx.Norigin; i++) {
		for(j=0; j<bpz.Norigin; j++) {
			vij[i][j]= 0.0- m[i* nvijz+ j];
		}
	}

	//free memory
	free1(F.data);
	free1float(m);
	free1float(d);
}

