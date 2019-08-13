#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include <time.h>


double tic(void)
{
/*
 * record current time
 */
	return (clock () / (double) CLOCKS_PER_SEC) ;
}

double toc(double t)
{
	double s = tic(); 
	return (s-t>0?s-t:0) ;
}


int nLineCount(FILE *infile)
{
/*
	count the lines of a text file
	without empty lines.
	effective lines has number in it.
	return number of non-empty lines.
*/
	char ch;
	int nLine;
	int bSymbol;
	bSymbol = 0;
	nLine = 0;
	while ((ch = getc(infile)) != EOF) {
		if (ch >= '0' && ch <= '9')
			bSymbol = 1;
		else if (ch == '\n' && bSymbol){
				nLine++;
				bSymbol = 0;
			}
	}
	fseek(infile, 0, SEEK_SET);
	return nLine;
}

int nLineCount_filename(char *filename)
{
	FILE *infile;
	int nLine;
	
	infile=efopen(filename, "r");
	nLine = nLineCount(infile);
	efclose(infile);	
	
	return nLine;
}

int nFileSize(FILE *infile)
{        
	int nSize;
	fseek(infile, 0, SEEK_END);
	nSize=ftell(infile);
        fseek(infile, 0, SEEK_SET);
        return nSize;
}

int nFileSize_filename(char *filename)
{        
	FILE *infile;
	int nSize;

	infile=efopen(filename, "r");
	nSize = nFileSize(infile);
	efclose(infile);

        return nSize;
}

int nStructCount(FILE *infile, int StructSize)
{
	int nCount;
	int nSize;
	nSize=nFileSize(infile);
	
	if(nSize%StructSize)
		printf("warning: file to count struct may wrong!\n");
	
	nCount=nSize/StructSize;
	
	return nCount;
}

float radian2angle(float radian)
{
	return radian/PI*180;
}

float angle2radian(float angle)
{
	return angle/180*PI;
}

void mpi_gatherv_construct(int start, int end, int numprocs, int ntask, int *recvcounts, int *displs)
{
/* for MPI gatherv function,  */	
	int i;
	int quotient,remainder;
	quotient=(end-start)/numprocs;
	remainder=(end-start)%numprocs;
	for(i=0;i<numprocs;i++) 
		recvcounts[i]=quotient;
	for(i=0;i<remainder;i++) 
		recvcounts[i]++;
	for(i=0;i<numprocs;i++) 
		recvcounts[i]=recvcounts[i]*ntask;
	displs[0]=0;	
	for(i=1;i<numprocs;i++) 
		displs[i]=displs[i-1]+recvcounts[i-1];
}

void normalization(float *trace, int lt)
{
	float amp=0.0;
	int i;
	for(i=0;i<lt;i++)
	{
		amp+=trace[i]*trace[i];
	}
	amp=(float)sqrt(amp);
	if(amp<1.0E-06)	
	{
		amp=1;
	}
	for(i=0;i<lt;i++)
	{
		trace[i]=trace[i]/amp;
	} 
}

void float1max(float *d, int n, float *max, int *pos)
{
	int i;
	float maxv;
	*max=d[0];
	*pos=0;
	for (i=1;i<n;i++){
		if(*max<d[i]){
			*max=d[i];
			*pos=i;
		}
	}
}

void float1min(float *d, int n, float *min, int *pos)
{
	int i;
	float minv;
	*min=d[0];
	*pos=0;
	for (i=1;i<n;i++){
		if(*min>d[i]){
			*min=d[i];
			*pos=i;
		}
	}
}

void zero1float(float *p, int n1)
{
	int i;
	for(i=0;i<n1;i++)
		p[i]=0.0;
}
void zero2float(float **p, int n1, int n2)
{
	int i, j;
	for(i=0;i<n2;i++)
		for(j=0;j<n1;j++)
			p[i][j]=0.0;
}

void zero2int(int **p, int n1, int n2)
{
	int i, j;
	for(i=0;i<n2;i++)
		for(j=0;j<n1;j++)
			p[i][j]=0;
}

float average(float *data, int n)
{
	int i;
	float sum;
	sum=0;
	for(i=0;i<n;i++)
		sum+=data[i];
	return sum/n;
}


float average_abs(float *data, int n)
{
	int i;
	float sum;
	sum=0;
	for(i=0;i<n;i++)
		sum+=fabs(data[i]);
	return sum/n;
}

int binary_find(float * arr, int low, int high, float key)
/*
	binary search for key
	http://www.cnblogs.com/coser/archive/2011/04/11/2013013.html 
*/
{
    if(arr[low] >= key) return low;
    if (low > high) return -1;
    int mid = (low + high) / 2;
    if (arr[mid] < key) return binary_find(arr,mid+1,high,key);
    else if(arr[mid] >= key){
        if(mid >= low && arr[mid-1]>=key)
            return binary_find(arr,low, mid-1, key);
        return mid;
    }
}

void fieldcpy(BGfield **vel_tmp,BGfield **vel,int n1,int n2)
{
	int i,j;
	for(i=0;i<n1;i++)
	  for(j=0;j<n2;j++)
		vel_tmp[i][j]=vel[i][j];
}

void field_assign(float **vel_field,BGfield **vel,int nx,int nz)
{
	int i,j;
	for(i=0;i<nx;i++)
	  for(j=0;j<nz;j++)
		vel[i][j].u=vel_field[i][j];
}

void zerofrechet(DField *sf)
{
	sf->pv = 0; sf->pe = 0; sf->pd = 0;
	sf->xv = 0; sf->xe = 0; sf->xd = 0;
	sf->tv = 0; sf->te = 0; sf->td = 0;
}
