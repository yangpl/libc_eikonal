#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cstd.h"

extern float h1, h2, h3;
extern int n1, n2, n3;
extern float x_source, y_source, z_source;
extern int niter, nfpi;
extern float epsilon;

void eikonal_solver(float*** Vnmo, float*** V0, float*** eta, float*** T);

float shooting_method(float x0, float y0, float z0, float x1, float y1, float z1, float Vnmoc, float V0c, float eta);

void main(int argc, char *argv[]){
	/*
	h1 = 10; h2 = 10; h3 = 5;
	n1 = 201; n2 =201; n3 = 101;
	x_source = 1202; y_source =886.2; z_source = 333.3;
	niter = 8; nfpi = 6;
	epsilon = 1.0e-2;
	*/
	initargs(argc,argv);
	if(!getparint("n1",&n1)) n1=201;
	if(!getparint("n2",&n2)) n2=201;
	if(!getparint("n3",&n3)) n3=101;
	if(!getparfloat("h1",&h1)) h1=10.0;
	if(!getparfloat("h2",&h2)) h2=10.0;
	if(!getparfloat("h3",&h3)) h3=5.0;
	if(!getparfloat("x_source",&x_source)) x_source=1202.0;
	if(!getparfloat("y_source",&y_source)) y_source=886.2;
	if(!getparfloat("z_source",&z_source)) z_source=333.3;
	if(!getparint("niter",&niter)) niter=3;
	if(!getparint("nfpi",&nfpi)) nfpi=3;
	if(!getparfloat("epsilon",&epsilon)) epsilon=1.0e-10;
	float*** Vnmo = alloc3float(n3,n2,n1);
	float*** V0 = alloc3float(n3,n2,n1);
	float*** eta = alloc3float(n3,n2,n1);
	float*** t = alloc3float(n3,n2,n1);
	float*** tt = alloc3float(n3,n2,n1);
	for(int i=0; i<n1; i++){
		for(int j=0; j<n2; j++){
			for(int k=0; k<n3; k++){
				Vnmo[i][j][k] = 2500.;
				V0[i][j][k] = 2000.;
				eta[i][j][k] = 0.2;
			}
		}
	}
	
	clock_t start,end;
	float CPU_time_used;
	start=clock();
	
	eikonal_solver(Vnmo, V0, eta, t);
			
	end=clock();
	CPU_time_used=(float)(end-start)/CLOCKS_PER_SEC;
	printf("CPU_time_used=%f\n",CPU_time_used);
	
	FILE *fp = fopen("traveltime.bin","wb+");
	fwrite(&t[0][0][0], n3*n2*n1*sizeof(float), 1, fp);
        fclose(fp);

	float x0 = x_source, y0= y_source, z0 = z_source;
	for(int i=0; i<n1; i++){
		for(int j=0; j<n2; j++){
			for(int k=0; k<n3; k++){
				float x1 = i*h1, y1 = j*h2, z1 = k*h3;
				float dx = x1-x0, dy = y1-y0, dz = z1-z0;
				float v = shooting_method(x0,y0,z0,x1,y1,z1,Vnmo[i][j][k],V0[i][j][k],eta[i][j][k]);
				float r = sqrt(dx*dx+dy*dy+dz*dz);
				tt[i][j][k] = r/v;		
			}
		}
	}
	
	float err = 0;
	for(int i=0; i<n1; i++)
		for(int j=0; j<n2; j++)
			for(int k=0; k<n3; k++)
				if(tt[i][j][k]>1.0e-10)
					err = fmax(err, fabs(t[i][j][k]-tt[i][j][k])/tt[i][j][k]);
	printf("err = %lf \n",err);

	free3float(Vnmo);
	free3float(V0);
	free3float(eta);
	free3float(t);
	free3float(tt);
}
