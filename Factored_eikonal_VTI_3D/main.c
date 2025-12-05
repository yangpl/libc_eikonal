#include "cstd.h"
#include "eikonal.h"

void eikonal_solver(eikonal_t *eik);

float shooting_method(float x0, float y0, float z0, float x1, float y1, float z1, float Vnmoc, float V0c, float eta);

void main(int argc, char *argv[])
{
  int i, j, k;
  eikonal_t *eik;
  
  initargs(argc,argv);

  eik = malloc(sizeof(eikonal_t));
  
  if(!getparint("n1",&eik->n1)) eik->n1=201;
  if(!getparint("n2",&eik->n2)) eik->n2=201;
  if(!getparint("n3",&eik->n3)) eik->n3=101;
  if(!getparfloat("h1",&eik->h1)) eik->h1=10.0;
  if(!getparfloat("h2",&eik->h2)) eik->h2=10.0;
  if(!getparfloat("h3",&eik->h3)) eik->h3=5.0;
  if(!getparfloat("x_source",&eik->x_source)) eik->x_source=1202.0;
  if(!getparfloat("y_source",&eik->y_source)) eik->y_source=886.2;
  if(!getparfloat("z_source",&eik->z_source)) eik->z_source=333.3;
  if(!getparint("niter",&eik->niter)) eik->niter=3;
  if(!getparint("nfpi",&eik->nfpi)) eik->nfpi=3;
  if(!getparfloat("epsilon",&eik->epsilon)) eik->epsilon=1.0e-10;

  eik->Vnmo = alloc3float(eik->n3,eik->n2,eik->n1);
  eik->V0 = alloc3float(eik->n3,eik->n2,eik->n1);
  eik->eta = alloc3float(eik->n3,eik->n2,eik->n1);
  eik->T = alloc3float(eik->n3,eik->n2,eik->n1);
  eik->TT = alloc3float(eik->n3,eik->n2,eik->n1);
  for(i=0; i<eik->n1; i++){
    for(j=0; j<eik->n2; j++){
      for(k=0; k<eik->n3; k++){
	eik->Vnmo[i][j][k] = 2500.;
	eik->V0[i][j][k] = 2000.;
	eik->eta[i][j][k] = 0.2;
      }
    }
  }
	
  clock_t start, end;
  float CPU_time_used;
  start=clock();	
  eikonal_solver(eik);			
  end=clock();
  CPU_time_used=(float)(end-start)/CLOCKS_PER_SEC;
  printf("CPU_time_used=%f\n",CPU_time_used);
	
  FILE *fp = fopen("traveltime.bin","wb+");
  fwrite(&eik->T[0][0][0], eik->n3*eik->n2*eik->n1*sizeof(float), 1, fp);
  fclose(fp);

  for(i=0; i<eik->n1; i++){
    for(j=0; j<eik->n2; j++){
      for(k=0; k<eik->n3; k++){
	float x1 = i*eik->h1;
	float y1 = j*eik->h2;
	float z1 = k*eik->h3;
	float dx = x1-eik->x_source;
	float dy = y1-eik->y_source;
	float dz = z1-eik->z_source;
	float v = shooting_method(eik->x_source, eik->y_source, eik->z_source, x1,y1,z1, eik->Vnmo[i][j][k], eik->V0[i][j][k], eik->eta[i][j][k]);
	float r = sqrt(dx*dx+dy*dy+dz*dz);
	eik->TT[i][j][k] = r/v;		
      }
    }
  }
	
  float error = 0;
  for(i=0; i<eik->n1; i++)
    for(j=0; j<eik->n2; j++)
      for(k=0; k<eik->n3; k++)
	if(eik->TT[i][j][k]>1.0e-10)
	  error = fmax(error, fabs(eik->T[i][j][k]-eik->TT[i][j][k])/eik->TT[i][j][k]);
  printf("err = %lf \n", error);

  free3float(eik->Vnmo);
  free3float(eik->V0);
  free3float(eik->eta);
  free3float(eik->T);
  free3float(eik->TT);
}
