#include "su.h"
#include "lsqr.h"


void lsqr(int m,int n,float *x,float *u,int itmax,float **a,int **ax,int *axl,int rank)
{
	printf("-----lsqr begins-----\n");
	int i,iter;
	float beta,alfa;
	float b1,aa,b,c,r,s,t1,t2;
	float phibar,rho,rhobar,phi,teta;
	float *v,*w;
	float derr=10E-6;
	v = ealloc1float(n);
	w = ealloc1float(n);

	for(i=0;i<n;i++) 
	{ 
	    x[i]=0;
	    v[i]=0;
	}
	normlz(m,u,&beta); 
	b1=beta;
	atupv(m,n,u,v,a,ax,axl); 
	normlz(n,v,&alfa);

	rhobar=alfa;
	phibar=beta;
	for(i=0;i<n;i++)
	    w[i]=v[i];
	r=phibar/b1;
	for(iter=0;iter<itmax;iter++)
	{	if(iter%200==0)
		printf("iter == %5d\n",iter);
		aa=-alfa; 
		for(i=0;i<m;i++) 
		    u[i]=aa*u[i]; 
		avpu(m,n,u,v,a,ax,axl); 
		normlz(m,u,&beta);
		b=-beta; 
		for(i=0;i<n;i++)
	       	    v[i]=b*v[i]; 
		atupv(m,n,u,v,a,ax,axl);
	       	normlz(n,v,&alfa);
		rho=sqrt(rhobar*rhobar+beta*beta);
		c=rhobar/rho;
		s=beta/rho;
		teta=s*alfa;
		rhobar=-c*alfa;
		phi=c*phibar;
		phibar=s*phibar;
		t1=phi/rho;
		t2=-teta/rho;
		for(i=0;i<n;i++)
	       	{ 
		    x[i]=t1*w[i]+x[i];
		    w[i]=t2*w[i]+v[i];
		}
		r=phibar/b1;
		if(phibar<derr || r<derr)
		{
		   printf("normal exit lsqr\n");
		   goto mmm;
		}
	}
mmm:	free1float(v);
	free1float(w);
	printf("-----lsqr completed-----\n");
}

void avpu(int m,int n,float *u,float *v,float **a,int **ax,int *axl)
{
	//u=u+a*v
	int i,j;
	float temp;
	for(i=0;i<m;i++)
	{
		temp=0;
		for(j=0;j<axl[i];j++)
			temp=temp+a[i][j]*v[ax[i][j]];
		u[i]=u[i]+temp;
	}
}

void atupv(int m,int n,float *u,float *v,float **a,int **ax,int *axl)
{
	int i,j,k;
	float *temp;
	temp= alloc1float(n);
	for(i=0;i<n;i++)
	    temp[i]=0;

	for(i=0;i<m;i++)
	    for(j=0;j<axl[i];j++)
		temp[ax[i][j]]+=a[i][j]*u[i];

	for(i=0;i<n;i++)
	    v[i]=v[i]+temp[i];

	free1float(temp);
}

void normlz(int n,float *x,float *s)
{
	int i;
	*s=0;
	for(i=0;i<n;i++)
	{
		*s=*s+pow(x[i],(int)(2));
	}
	*s=sqrt(*s);
	for(i=0;i<n;i++)
	{
		x[i]=x[i]/(*s);
	}
}
