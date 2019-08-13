#include "su.h"
#include "stereo_struct.h"

void write1(void *data, int size, int N, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "wb+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(data, size, N, fp);
	fclose(fp);
}

void write1float(float *data, int N, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "wb+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(&data[0], sizeof(float), N, fp);
	fclose(fp);
}

void write1float_append(float *data, int N, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "ab+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(&data[0], sizeof(float), N, fp);
	fclose(fp);
}

void write2float(float **data, int Nx, int Nz, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "wb+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(&data[0][0], sizeof(float), Nx*Nz, fp);
	fclose(fp);
}

void write2int(int **data, int Nx, int Nz, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "wb+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(&data[0][0], sizeof(int), Nx*Nz, fp);
	fclose(fp);
}

void write3float(float ***data, int Nx, int Ny, int Nz, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "wb+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(&data[0][0][0], sizeof(float), Nx*Ny*Nz, fp);
	fclose(fp);
}

void write1txt(float *data, int N, char *name)
{
	FILE *fp;
	int i;
	if((fp=fopen(name, "w+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	for(i=0; i<N; i++)
	fprintf(fp, "%f\n", data[i]);
	fclose(fp);
}

void write1double(double *data, int N, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "wb+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(&data[0], sizeof(double), N, fp);
	fclose(fp);
}

void write2double(double **data, int Nx, int Nz, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "wb+"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fwrite(&data[0][0], sizeof(double), Nx*Nz, fp);
	fclose(fp);
}

void read1(void *data, int size, int Nx, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "r"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fread(data, size, Nx, fp);
	fclose(fp);
}

void read1float(float *data, int Nx, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "r"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fread(&data[0], sizeof(float), Nx, fp);
	fclose(fp);
}

void read2float(float **data, int Nx, int Nz, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "r"))==NULL)
	{
		printf("open file    %s   with r mode failed\n", name);
		exit(0);
	}
	fread(&data[0][0], sizeof(float), Nx*Nz, fp);
	fclose(fp);
}

void read3float(float ***data, int Nx, int Ny, int Nz, char *name)
{
	FILE *fp;
	if((fp=fopen(name, "r"))==NULL)
	{
		printf("open file    %s   failed\n", name);
		exit(0);
	}
	fread(&data[0][0][0], sizeof(float), Nx*Ny*Nz, fp);
	fclose(fp);
}

void write2field(BGfield **vel,int Nx,int Ny,char *name)
{
	int i,j;
	FILE *fp;
	fp=fopen(name,"wb+");
	for(i=0;i<Nx;i++)
	  for(j=0;j<Ny;j++)
		fwrite(&vel[i][j].u,sizeof(float),1,fp);
	fclose(fp);
}




