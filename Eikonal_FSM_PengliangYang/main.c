/* Demo implmentation for fast sweeping method in 2D
 * Reference: Zhao, Hongkai. "A fast sweeping method for eikonal equations." 
 *            Mathematics of computation 74.250 (2005): 603-627.
 *-----------------------------------------------------------------------
 *
 * Copyright (c) 2024 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------
 */
#include "cstd.h"

int main(int argc, char *argv[])
{
  int nx, nz;
  float dx, dz;
  int i, j, iter, niter;
  
  initargs(argc, argv);//
  if(!getparint("niter", &niter)) niter = 1;/* niter */
  if(!getparint("nx", &nx)) nx = 1000;
  if(!getparint("nz", &nz)) nz = 1000;
  if(!getparfloat("dx", &dx)) dx = 10.;
  if(!getparfloat("dz", &dz)) dz = 10.;

  int x0 = nx/2;
  int z0 = nz/2;
  float **s = alloc2float(nz, nx);
  float **T = alloc2float(nz, nx);
  for(i=0; i<nx; i++){
    for(j=0; j<nz; j++){
      s[i][j] = 1./1000.0;//slowness
      T[i][j] = 1e5;//set large numbers for all unknowns
    }
  }

  float a, b, a_b, tmp;
  int im1, ip1, jm1, jp1;
  for(iter=0; iter<niter; iter++){
    T[x0][z0] = 0.0;//set 0 at the source

    for(i=0; i<nx; i++){//forward sweep in x
      im1 = MAX(i-1, 0);
      ip1 = MIN(i+1, nx-1);
      for(j=0; j<nz; j++){//forward sweep in z
	jm1 = MAX(j-1, 0);
	jp1 = MIN(j+1, nz-1);
	  
	a = MIN(T[im1][j], T[ip1][j]);
	b = MIN(T[i][jm1], T[i][jp1]);
	tmp = s[i][j]*dx;
	a_b = a-b;
	if(fabs(a_b)>= tmp) T[i][j] = MIN(a, b) + tmp;
	else T[i][j] = 0.5*(a+b + sqrt(2.*tmp*tmp - (a_b)*(a_b)));	
      }
    }

    for(i=nx-1; i>=0; i--){//backward sweep in x
      im1 = MAX(i-1, 0);
      ip1 = MIN(i+1, nx-1);
      for(j=0; j<nz; j++){//forward sweep in z
	jm1 = MAX(j-1, 0);
	jp1 = MIN(j+1, nz-1);
	  
	a = MIN(T[im1][j], T[ip1][j]);
	b = MIN(T[i][jm1], T[i][jp1]);
	tmp = s[i][j]*dx;
	a_b = a-b;
	if(fabs(a_b)>= tmp) T[i][j] = MIN(a, b) + tmp;
	else T[i][j] = 0.5*(a+b + sqrt(2.*tmp*tmp - (a_b)*(a_b)));	
      }
    }

    for(i=0; i<nx; i++){//forward sweep in x
      im1 = MAX(i-1, 0);
      ip1 = MIN(i+1, nx-1);
      for(j=nz-1; j>=0; j--){//backward sweep in z
	jm1 = MAX(j-1, 0);
	jp1 = MIN(j+1, nz-1);
	  
	a = MIN(T[im1][j], T[ip1][j]);
	b = MIN(T[i][jm1], T[i][jp1]);
	tmp = s[i][j]*dx;
	a_b = a-b;
	if(fabs(a_b)>= tmp) T[i][j] = MIN(a, b) + tmp;
	else T[i][j] = 0.5*(a+b + sqrt(2.*tmp*tmp - (a_b)*(a_b)));	
      }
    }

    for(i=nx-1; i>=0; i--){//backward sweep in x
      im1 = MAX(i-1, 0);
      ip1 = MIN(i+1, nx-1);
      for(j=nz-1; j>=0; j--){//backward sweep in z
	jm1 = MAX(j-1, 0);
	jp1 = MIN(j+1, nz-1);
	  
	a = MIN(T[im1][j], T[ip1][j]);
	b = MIN(T[i][jm1], T[i][jp1]);
	tmp = s[i][j]*dx;
	a_b = a-b;
	if(fabs(a_b)>= tmp) T[i][j] = MIN(a, b) + tmp;
	else T[i][j] = 0.5*(a+b + sqrt(2.*tmp*tmp - (a_b)*(a_b)));	
      }
    }

  }//end for it

  FILE *fp = fopen("traveltime.bin","wb+");
  fwrite(&T[0][0], nz*nx*sizeof(float), 1, fp);
  fclose(fp);
  
  free2float(s);
  free2float(T);
}
