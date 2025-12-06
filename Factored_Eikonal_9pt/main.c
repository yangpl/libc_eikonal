#include "cstd.h"
#include "eikonal.h"

float shooting_method(float x0, float y0, float z0, float x1, float y1, float z1, float Vnmoc, float V0c, float eta);

void eikonal_solver(eikonal_t* eik);

void main(int argc, char *argv[])
{
  int i, j, k;
  eikonal_t* eik;

  initargs(argc,argv);

  eik = (eikonal_t*)malloc(sizeof(eikonal_t));
  if(!getparint("n1",&eik->n1)) eik->n1=101;
  if(!getparint("n2",&eik->n2)) eik->n2=101;
  if(!getparint("n3",&eik->n3)) eik->n3=51;
  if(!getparfloat("h1",&eik->h1)) eik->h1=100.0;
  if(!getparfloat("h2",&eik->h2)) eik->h2=100.0;
  if(!getparfloat("h3",&eik->h3)) eik->h3=100.0;
  if(!getparfloat("x_source",&eik->x_source)) eik->x_source=6022.31;
  if(!getparfloat("y_source",&eik->y_source)) eik->y_source=5223.77;
  if(!getparfloat("z_source",&eik->z_source)) eik->z_source=1764.3;
  if(!getparint("niter",&eik->niter)) eik->niter=1;
  if(!getparfloat("epsilon",&eik->epsilon)) eik->epsilon=1.0e-4;

  int n1 = eik->n1, n2 = eik->n2, n3 = eik->n3;
  float h1 = eik->h1, h2 = eik->h2, h3 = eik->h3;
  float x_source = eik->x_source, y_source = eik->y_source, z_source = eik->z_source;

  eik->Vnmo = alloc3float(n3,n2,n1);
  eik->V0 = alloc3float(n3,n2,n1);
  eik->eta = alloc3float(n3,n2,n1);
  eik->T = alloc3float(n3,n2,n1);
  eik->TT = alloc3float(n3,n2,n1);

  char *Vnmo = NULL, *V0 = NULL, *eta = NULL;
  if(!getparstring("Vnmo",&Vnmo)) Vnmo="Vnmo.bin";
  if(!getparstring("V0",&V0)) V0="V0.bin";
  if(!getparstring("eta",&eta)) eta="eta.bin";
  FILE* fp;

  fp = fopen(Vnmo,"rb");
  if(fp==NULL){
    printf("can not open %s\n",Vnmo);
    printf("use default constant Vnmo=2500 m/s\n");
    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  eik->Vnmo[i][j][k] = 2500.0;
  }else{
    fread(&eik->Vnmo[0][0][0], n3*n2*n1*sizeof(float), 1, fp);
    fclose(fp);
  }

  fp = fopen(V0,"rb");
  if(fp==NULL){
    printf("can not open %s\n",V0);
    printf("use default constant V0=2000 m/s\n");
    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  eik->V0[i][j][k] = 2000.0;
  }else{
    fread(&eik->V0[0][0][0], n3*n2*n1*sizeof(float), 1, fp);
    fclose(fp);
  }

  fp = fopen(eta,"rb");
  if(fp==NULL){
    printf("can not open %s\n",eta);
    printf("use default constant eta=0.2\n");
    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  eik->eta[i][j][k] = 0.2;
  }else{
    fread(&eik->eta[0][0][0], n3*n2*n1*sizeof(float), 1, fp);
    fclose(fp);
  }
	
  clock_t start,end;
  float CPU_time_used;
  start=clock();
	
  eikonal_solver(eik);
			
  end=clock();
  CPU_time_used=(float)(end-start)/CLOCKS_PER_SEC;
  printf("CPU_time_used=%f\n",CPU_time_used);
	
  fp = fopen("traveltime.bin","wb+");
  fwrite(&eik->T[0][0][0], n3*n2*n1*sizeof(float), 1, fp);
  fclose(fp);


  // calculate reference solution by shooting method
  for(i=0; i<n1; i++)
    for(j=0; j<n2; j++)
      for(k=0; k<n3; k++){
	float x1 = i*h1, y1 = j*h2, z1 = k*h3;
	float dx = x1 - x_source, dy = y1 - y_source, dz = z1 - z_source;
	float v = shooting_method(x_source, y_source, z_source, x1, y1, z1, eik->Vnmo[i][j][k], eik->V0[i][j][k], eik->eta[i][j][k]);
	float r = sqrt(dx*dx + dy*dy + dz*dz);
	eik->TT[i][j][k] = r / v;
      }
	
  // calculate MAX_ERR and AVERAGE_ERR
  float err_l_infinity = 0, err_l1 = 0;
  for(i=0; i<n1; i++){
    for(j=0; j<n2; j++){
      for(k=0; k<n3; k++){
	float err = fabs(eik->T[i][j][k] - eik->TT[i][j][k]);
	err_l_infinity = fmax(err_l_infinity, err);
	err_l1 += err;
      }
    }
  }
  err_l1 = err_l1 / (n1 * n2 * n3);
  printf("MAX_ERR = %f, AVERAGE_ERR = %f\n", err_l_infinity, err_l1);
	
	
  free3float(eik->Vnmo); free3float(eik->V0); free3float(eik->eta);
  free3float(eik->T); free3float(eik->TT);
  free(eik);
}
