#include "cstd.h"
#include "eikonal.h"

#define HUGE 1e10

void compute_ray_parameters(float*** t, float*** dtdx, float*** dtdy, float*** dtdz, float h1, float h2, float h3, int n1, int n2, int n3)
{
  int i, j, k;
  
  for(i=0; i<n1; i++){
    for(j=0; j<n2; j++){
      for(k=1; k<n3-1; k++){
	dtdz[i][j][k] = (t[i][j][k+1]-t[i][j][k-1])/(2*h3);//1st derivative approximated by centered FD
      }
      dtdz[i][j][0] = (-3*t[i][j][0]+4*t[i][j][1]-t[i][j][2])/(2*h3);//1st derivative approximated by one-sided FD
      dtdz[i][j][n3-1] = (3*t[i][j][n3-3]-4*t[i][j][n3-2]+t[i][j][n3-1])/(2*h3);//1st derivative approximated by one-sided FD
    }
  }

  for(i=0; i<n1; i++){
    for(k=0; k<n3; k++){
      for(j=1; j<n2-1; j++){
	dtdy[i][j][k] = (t[i][j+1][k]-t[i][j-1][k])/(2*h2);//1st derivative approximated by centered FD
      }
      dtdy[i][0][k] = (-3*t[i][0][k]+4*t[i][1][k]-t[i][2][k])/(2*h2);//1st derivative approximated by one-sided FD
      dtdy[i][n2-1][k] = (3*t[i][n2-3][k]-4*t[i][n2-2][k]+t[i][n2-1][k])/(2*h2);//1st derivative approximated by one-sided FD
    }
  }

  for(j=0; j<n2; j++){
    for(k=0; k<n3; k++){
      for(i=1; i<n1-1; i++){
	dtdx[i][j][k] = (t[i+1][j][k]-t[i-1][j][k])/(2*h1);//1st derivative approximated by centered FD
      }
      dtdx[0][j][k] = (-3*t[0][j][k]+4*t[1][j][k]-t[2][j][k])/(2*h1);//1st derivative approximated by one-sided FD
      dtdx[n1-1][j][k] = (3*t[n1-3][j][k]-4*t[n1-2][j][k]+t[n1-1][j][k])/(2*h1);//1st derivative approximated by one-sided FD
    }
  }
}

//solving 2nd order polynomial equation: x*x+q1*x+q0=0
//root:only record the real root
//nd:the number of real root
void find_root_polynomial2(float q1,float q0,float* root,int* nd)
{
  float delta = q1*q1-4*q0;
  root[0] = 0;
  root[1] = 0;
  if(delta>0){
    root[0] = (-q1+sqrt(delta))/2.;
    root[1] = (-q1-sqrt(delta))/2.;
    *nd = 2;
  }else if(fabs(delta)<1.0e-12){
    root[0] = -q1/2.;
    *nd = 1;
  }else{
    *nd = 0;
  }
}



void find_upwind_stencil(eikonal_t *eik, int xc, int yc, int zc, float* T0x, float* T0y, float* T0z,
			 float* taux, float* tauy, float* tauz, int* sx, int* sy, int* sz)
{
  if(xc == 0){
    *T0x = eik->T0[xc+1][yc][zc];
    *taux = eik->tau[xc+1][yc][zc];
    *sx = -1;
  }else if(xc == eik->n1-1){
    *T0x = eik->T0[xc-1][yc][zc];
    *taux = eik->tau[xc-1][yc][zc];
    *sx = 1;
  }else{
    if((eik->tau[xc-1][yc][zc]*eik->T0[xc-1][yc][zc]) < (eik->tau[xc+1][yc][zc]*eik->T0[xc+1][yc][zc])){
      *T0x = eik->T0[xc-1][yc][zc];
      *taux = eik->tau[xc-1][yc][zc];
      *sx = 1;
    }else{
      *T0x = eik->T0[xc+1][yc][zc];
      *taux = eik->tau[xc+1][yc][zc];
      *sx = -1;
    }
  }

  if(yc == 0){
    *T0y = eik->T0[xc][yc+1][zc];
    *tauy = eik->tau[xc][yc+1][zc];
    *sy = -1;
  }else if(yc == eik->n2-1){
    *T0y = eik->T0[xc][yc-1][zc];
    *tauy = eik->tau[xc][yc-1][zc];
    *sy = 1;
  }else{
    if((eik->tau[xc][yc-1][zc]*eik->T0[xc][yc-1][zc]) < (eik->tau[xc][yc+1][zc]*eik->T0[xc][yc+1][zc])){
      *T0y = eik->T0[xc][yc-1][zc];
      *tauy = eik->tau[xc][yc-1][zc];
      *sy = 1;
    }else{
      *T0y = eik->T0[xc][yc+1][zc];
      *tauy = eik->tau[xc][yc+1][zc];
      *sy = -1;
    }
  }

  if(zc == 0){
    *T0z = eik->T0[xc][yc][zc+1];
    *tauz = eik->tau[xc][yc][zc+1];
    *sz = -1;
  }else if(zc == eik->n3-1){
    *T0z = eik->T0[xc][yc][zc-1];
    *tauz = eik->tau[xc][yc][zc-1];
    *sz = 1;
  }else{
    if((eik->tau[xc][yc][zc-1]*eik->T0[xc][yc][zc-1]) < (eik->tau[xc][yc][zc+1]*eik->T0[xc][yc][zc+1])){
      *T0z = eik->T0[xc][yc][zc-1];
      *tauz = eik->tau[xc][yc][zc-1];
      *sz = 1;
    }else{
      *T0z = eik->T0[xc][yc][zc+1];
      *tauz = eik->tau[xc][yc][zc+1];
      *sz = -1;
    }
  }

}

bool is_Causal_root_2D(float root, float taux, float tauy, float h1, float h2, float px0c, float py0c, int sx, int sy, float T0c){
  float a = T0c*(root-taux)/h1 + px0c*root*sx;
  float b = T0c*(root-tauy)/h2 + py0c*root*sy;
  if(a>0 && b>0) return true;
  else return false;
}

bool is_Causal_root_3D(float root, float taux, float tauy, float tauz, float h1, float h2, float h3, float px0c, float py0c, float pz0c, int sx, int sy, int sz, float T0c){
  float a = T0c*(root-taux)/h1 + px0c*root*sx;
  float b = T0c*(root-tauy)/h2 + py0c*root*sy;
  float c = T0c*(root-tauz)/h3 + pz0c*root*sz;
  if(a>0 && b>0 && c>0) return true;
  else return false;
}


float trisolver_xy(eikonal_t *eik, float Vnmoc, float V0c, float etac, float T0c, float px0c, float py0c, float pz0c, float rhsc,
		   float taux, float tauy, float tauz, float T0x, float T0y, float T0z, int sx, int sy, int sz){
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c, T0c2 = T0c*T0c;
  float px0c2 = px0c*px0c, py0c2 = py0c*py0c, pz0c2 = pz0c*pz0c;
  float taucx = (T0c*taux + eik->h1*sqrt(rhsc/(Vnmoc2*(1+2*etac))))/(T0c + eik->h1*fabs(px0c));
  float taucy = (T0c*tauy + eik->h2*sqrt(rhsc/(Vnmoc2*(1+2*etac))))/(T0c + eik->h2*fabs(py0c));
  if(taux == HUGE) return taucy;
  if(tauy == HUGE) return taucx;

  float root[2];
  float result = HUGE;
  int nd = 0, flag = 1;
  float a = Vnmoc2*(1+2*etac)*( T0c2/(eik->h1*eik->h1) + px0c2 +2*sx*T0c*px0c/eik->h1 )
    + Vnmoc2*(1+2*etac)*( T0c2/(eik->h2*eik->h2) + py0c2 +2*sy*T0c*py0c/eik->h2 );
  float b = -2*Vnmoc2*(1+2*etac)*taux*( T0c2/(eik->h1*eik->h1) +sx*T0c*px0c/eik->h1 )
    -2*Vnmoc2*(1+2*etac)*tauy*( T0c2/(eik->h2*eik->h2) +sy*T0c*py0c/eik->h2 );
  float c = Vnmoc2*(1+2*etac)*taux*taux*T0c2/(eik->h1*eik->h1) + Vnmoc2*(1+2*etac)*tauy*tauy*T0c2/(eik->h2*eik->h2) - rhsc;
  find_root_polynomial2(b/a, c/a, root, &nd);
  for(int i=0; i<nd; i++){
    if(is_Causal_root_2D(root[i], taux, tauy, eik->h1, eik->h2, px0c, py0c, sx, sy, T0c)){
      flag = 0;
      result = fmin(result, root[i]);
    }
  }
  if(flag) result = fmin(taucx,taucy);

  return result;
}

float trisolver_xz(eikonal_t *eik, float Vnmoc, float V0c, float etac, float T0c, float px0c, float py0c, float pz0c, float rhsc,
		   float taux, float tauy, float tauz, float T0x, float T0y, float T0z, int sx, int sy, int sz){
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c, T0c2 = T0c*T0c;
  float px0c2 = px0c*px0c, py0c2 = py0c*py0c, pz0c2 = pz0c*pz0c;
  float taucx = (T0c*taux + eik->h1*sqrt(rhsc/(Vnmoc2*(1+2*etac))))/(T0c + eik->h1*fabs(px0c));
  float taucz = (T0c*tauz + eik->h3*sqrt(rhsc/V0c2))/(T0c + eik->h3*fabs(pz0c));
  if(taux == HUGE) return taucz;
  if(tauz == HUGE) return taucx;

  float root[2];
  float result = HUGE;
  int nd = 0, flag = 1;
  float a = Vnmoc2*(1+2*etac)*( T0c2/(eik->h1*eik->h1) + px0c2 +2*sx*T0c*px0c/eik->h1 )
    + V0c2*( T0c2/(eik->h3*eik->h3) + pz0c2 + 2*sz*T0c*pz0c/eik->h3 );
  float b = -2*Vnmoc2*(1+2*etac)*taux*( T0c2/(eik->h1*eik->h1) +sx*T0c*px0c/eik->h1 )
    -2*V0c2*tauz*( T0c2/(eik->h3*eik->h3) + sz*T0c*pz0c/eik->h3 );
  float c = Vnmoc2*(1+2*etac)*taux*taux*T0c2/(eik->h1*eik->h1) + V0c2*tauz*tauz*T0c2/(eik->h3*eik->h3) - rhsc;
  find_root_polynomial2(b/a, c/a, root, &nd);
  for(int i=0; i<nd; i++){
    if(is_Causal_root_2D(root[i], taux, tauz, eik->h1, eik->h3, px0c, pz0c, sx, sz, T0c)){
      flag = 0;
      result = fmin(result, root[i]);
    }
  }
  if(flag) result = fmin(taucx,taucz);

  return result;
}

float trisolver_yz(eikonal_t *eik, float Vnmoc, float V0c, float etac, float T0c, float px0c, float py0c, float pz0c, float rhsc,
		   float taux, float tauy, float tauz, float T0x, float T0y, float T0z, int sx, int sy, int sz){
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c, T0c2 = T0c*T0c;
  float px0c2 = px0c*px0c, py0c2 = py0c*py0c, pz0c2 = pz0c*pz0c;
  float taucy = (T0c*tauy + eik->h2*sqrt(rhsc/(Vnmoc2*(1+2*etac))))/(T0c + eik->h2*fabs(py0c));
  float taucz = (T0c*tauz + eik->h3*sqrt(rhsc/V0c2))/(T0c + eik->h3*fabs(pz0c));
  if(tauy == HUGE) return taucz;
  if(tauz == HUGE) return taucy;

  float root[2];
  float result = HUGE;
  int nd = 0, flag = 1;
  float a = Vnmoc2*(1+2*etac)*( T0c2/(eik->h2*eik->h2) + py0c2 +2*sy*T0c*py0c/eik->h2 )
    + V0c2*( T0c2/(eik->h3*eik->h3) + pz0c2 + 2*sz*T0c*pz0c/eik->h3 );
  float b = -2*Vnmoc2*(1+2*etac)*tauy*( T0c2/(eik->h2*eik->h2) +sy*T0c*py0c/eik->h2 )
    -2*V0c2*tauz*( T0c2/(eik->h3*eik->h3) + sz*T0c*pz0c/eik->h3 );
  float c = Vnmoc2*(1+2*etac)*tauy*tauy*T0c2/(eik->h2*eik->h2) + V0c2*tauz*tauz*T0c2/(eik->h3*eik->h3) - rhsc;
  find_root_polynomial2(b/a, c/a, root, &nd);
  for(int i=0; i<nd; i++){
    if(is_Causal_root_2D(root[i], tauy, tauz, eik->h2, eik->h3, py0c, pz0c, sy, sz, T0c)){
      flag = 0;
      result = fmin(result, root[i]);
    }
  }
  if(flag) result = fmin(taucy,taucz);

  return result;
}


void stencil_solver_3D(eikonal_t *eik, int xc, int yc, int zc)
{
  float Vnmoc = eik->Vnmo[xc][yc][zc], V0c = eik->V0[xc][yc][zc], etac = eik->eta[xc][yc][zc], rhsc = eik->rhs[xc][yc][zc];
  float T0c = eik->T0[xc][yc][zc], px0c = eik->px0[xc][yc][zc], py0c = eik->py0[xc][yc][zc], pz0c = eik->pz0[xc][yc][zc];
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c, T0c2 = T0c*T0c;
  float px0c2 = px0c*px0c, py0c2 = py0c*py0c,pz0c2 = pz0c*pz0c;
  int sx, sy, sz;
  float T0x, T0y, T0z, taux, tauy, tauz, tauc = HUGE;

  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
	if((xc == eik->shotx[i]) && yc == eik->shoty[j] && zc == eik->shotz[k]) return;

  find_upwind_stencil(eik, xc, yc, zc, &T0x, &T0y, &T0z,  &taux, &tauy, &tauz, &sx, &sy, &sz);

  if((taux == HUGE) && (tauy == HUGE) && (tauz == HUGE)) return;

  if(taux == HUGE)
    tauc = trisolver_yz(eik, Vnmoc, V0c, etac, T0c, px0c, py0c, pz0c, rhsc, taux, tauy, tauz, T0x, T0y, T0z, sx, sy, sz);
  else if(tauy == HUGE)
    tauc = trisolver_xz(eik, Vnmoc, V0c, etac, T0c, px0c, py0c, pz0c, rhsc, taux, tauy, tauz, T0x, T0y, T0z, sx, sy, sz);
  else if(tauz == HUGE)
    tauc = trisolver_xy(eik, Vnmoc, V0c, etac, T0c, px0c, py0c, pz0c, rhsc, taux, tauy, tauz, T0x, T0y, T0z, sx, sy, sz);
  else{
    float root[2];
    int nd = 0, flag = 1;
    float a = Vnmoc2*(1+2*etac)*( T0c2/(eik->h1*eik->h1) + px0c2 +2*sx*T0c*px0c/eik->h1 )
      + Vnmoc2*(1+2*etac)*( T0c2/(eik->h2*eik->h2) + py0c2 +2*sy*T0c*py0c/eik->h2 )
      + V0c2*( T0c2/(eik->h3*eik->h3) + pz0c2 + 2*sz*T0c*pz0c/eik->h3 );
    float b = -2*Vnmoc2*(1+2*etac)*taux*( T0c2/(eik->h1*eik->h1) +sx*T0c*px0c/eik->h1 )
      -2*Vnmoc2*(1+2*etac)*tauy*( T0c2/(eik->h2*eik->h2) +sy*T0c*py0c/eik->h2 )
      -2*V0c2*tauz*( T0c2/(eik->h3*eik->h3) + sz*T0c*pz0c/eik->h3 );
    float c = Vnmoc2*(1+2*etac)*taux*taux*T0c2/(eik->h1*eik->h1) + Vnmoc2*(1+2*etac)*tauy*tauy*T0c2/(eik->h2*eik->h2)
      + V0c2*tauz*tauz*T0c2/(eik->h3*eik->h3) - rhsc;
    find_root_polynomial2(b/a, c/a, root, &nd);
    for(int i=0; i<nd; i++){
      if(is_Causal_root_3D(root[i], taux, tauy, tauz, eik->h1, eik->h2, eik->h3, px0c, py0c, pz0c, sx, sy, sz, T0c)){
	tauc = fmin(tauc,root[i]);
	flag = 0;
      }
    }

    if(flag){
      float r1,r2,r3;
      r1 = trisolver_yz(eik, Vnmoc, V0c, etac, T0c, px0c, py0c, pz0c, rhsc, taux, tauy, tauz, T0x, T0y, T0z, sx, sy, sz);
      r2 = trisolver_xz(eik, Vnmoc, V0c, etac, T0c, px0c, py0c, pz0c, rhsc, taux, tauy, tauz, T0x, T0y, T0z, sx, sy, sz);
      r3 = trisolver_xy(eik, Vnmoc, V0c, etac, T0c, px0c, py0c, pz0c, rhsc, taux, tauy, tauz, T0x, T0y, T0z, sx, sy, sz);
      tauc = fmin(r1, fmin(r2,r3));
    }
  }
  eik->tau[xc][yc][zc] = fmin(eik->tau[xc][yc][zc], tauc);
}

void fast_sweeping(eikonal_t *eik)
{
  float maxval;
  int i, j, k;

  int n1 = eik->n1;
  int n2 = eik->n2;
  int n3 = eik->n3;
  
  float*** arr = alloc3float(n3,n2,n1);
  for(int iter=0; iter<eik->niter; iter++){
    printf("# %d-th inner loop\n",iter+1);
    memcpy(&arr[0][0][0], &eik->tau[0][0][0], n1*n2*n3*sizeof(float));
    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  stencil_solver_3D(eik, i, j, k);

    for(i=0; i<n1; i++)
      for(j=n2-1; j>=0; j--)
	for(k=0; k<n3; k++)
	  stencil_solver_3D(eik, i, j, k);

    for(i=n1-1; i>=0; i--)
      for(j=n2-1; j>=0; j--)
	for(k=0; k<n3; k++)
	  stencil_solver_3D(eik, i, j, k);

    for(i=n1-1; i>=0; i--)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  stencil_solver_3D(eik, i, j, k);

    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=n3-1; k>=0; k--)
	  stencil_solver_3D(eik, i, j, k);

    for(i=0; i<n1; i++)
      for(j=n2-1; j>=0; j--)
	for(k=n3-1; k>=0; k--)
	  stencil_solver_3D(eik, i, j, k);

    for(i=n1-1; i>=0; i--)
      for(j=n2-1; j>=0; j--)
	for(k=n3-1; k>=0; k--)
	  stencil_solver_3D(eik, i, j, k);

    for(i=n1-1; i>=0; i--)
      for(j=0; j<n2; j++)
	for(k=n3-1; k>=0; k--)
	  stencil_solver_3D(eik, i, j, k);

    maxval = 0;
    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  maxval = fmax(maxval, fabs(arr[i][j][k]-eik->tau[i][j][k]));
    if(maxval<eik->epsilon) break;
  }
  free3float(arr);
}


void eikonal_solver(eikonal_t *eik)
{
  float maxval;
  float ***arr;

  int n1 = eik->n1;
  int n2 = eik->n2;
  int n3 = eik->n3;
  int i = (int)(eik->x_source/eik->h1);
  int j = (int)(eik->y_source/eik->h2);
  int k = (int)(eik->z_source/eik->h3);

  eik->shotx[0] = i; eik->shotx[1] = i+1;
  eik->shoty[0] = j; eik->shoty[1] = j+1;
  eik->shotz[0] = k; eik->shotz[1] = k+1;
  eik->T0 = alloc3float(n3,n2,n1);
  eik->px0 = alloc3float(n3,n2,n1);
  eik->py0 = alloc3float(n3,n2,n1);
  eik->pz0 = alloc3float(n3,n2,n1);
  eik->tau = alloc3float(n3,n2,n1);
  eik->px = alloc3float(n3,n2,n1);
  eik->py = alloc3float(n3,n2,n1);
  eik->pz = alloc3float(n3,n2,n1);
  eik->rhs = alloc3float(n3,n2,n1);
  arr = alloc3float(n3,n2,n1);

  for(i=0; i<n1; i++){
    for(j=0; j<n2; j++){
      for(k=0; k<n3; k++){
	eik->rhs[i][j][k] = 1.0;//initialize rhs=1.
	eik->T[i][j][k] = 0.0;
      }
    }
  }

  for(int loop=0; loop<eik->nfpi; loop++){
    printf("------# %d-th outer loop----\n", loop+1);
    memcpy(&arr[0][0][0], &eik->T[0][0][0], n1*n2*n3*sizeof(float));
    // calculate T0 && grad((T0)
    float Vnmos = eik->Vnmo[eik->shotx[0]][eik->shoty[0]][eik->shotz[0]];
    float V0s = eik->V0[eik->shotx[0]][eik->shoty[0]][eik->shotz[0]];
    float etas = eik->eta[eik->shotx[0]][eik->shoty[0]][eik->shotz[0]];
    float Vnmos2 = Vnmos*Vnmos;
    float V0s2 = V0s*V0s;
    float a0 = Vnmos2*(1+2*etas), b0 = a0, c0 = V0s2;
    for(i=0; i<n1; i++){
      for(j=0; j<n2; j++){
	for(k=0; k<n3; k++){
	  float dx = i*eik->h1-eik->x_source;
	  float dy = j*eik->h2-eik->y_source;
	  float dz = k*eik->h3-eik->z_source;
	  float temp = b0*c0*dx*dx+a0*c0*dy*dy+a0*b0*dz*dz;
	  eik->T0[i][j][k] = sqrt(temp/(a0*b0*c0));
	  if(temp>1.0e-10){
	    eik->px0[i][j][k] = sqrt(b0*c0/a0)*dx/sqrt(temp);
	    eik->py0[i][j][k] = sqrt(a0*c0/b0)*dy/sqrt(temp);
	    eik->pz0[i][j][k] = sqrt(a0*b0/c0)*dz/sqrt(temp);
	  }else{
	    eik->px0[i][j][k] = 0;
	    eik->py0[i][j][k] = 0;
	    eik->pz0[i][j][k] = 0;
	  }
	}
      }
    }
    //initialize scaling factor tau: set tau=1 at source and HUGE everywhere else
    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  eik->tau[i][j][k] = HUGE;
    for(i=0; i<2; i++)
      for(j=0; j<2; j++)
	for(k=0; k<2; k++)
	  eik->tau[eik->shotx[i]][eik->shoty[j]][eik->shotz[k]] = 1.;
    fast_sweeping(eik);

    // calculate rhs && T
    compute_ray_parameters(eik->tau, eik->px, eik->py, eik->pz, eik->h1, eik->h2, eik->h3, eik->n1, eik->n2, eik->n3);
    for(i=0; i<n1; i++){
      for(j=0; j<n2; j++){
	for(k=0; k<n3; k++){
	  eik->px[i][j][k] = eik->px0[i][j][k]*eik->tau[i][j][k] + eik->px[i][j][k]*eik->T0[i][j][k]; 
	  eik->py[i][j][k] = eik->py0[i][j][k]*eik->tau[i][j][k] + eik->py[i][j][k]*eik->T0[i][j][k]; 
	  eik->pz[i][j][k] = eik->pz0[i][j][k]*eik->tau[i][j][k] + eik->pz[i][j][k]*eik->T0[i][j][k];
	  float Vnmoc = eik->Vnmo[i][j][k];
	  float V0c = eik->V0[i][j][k];
	  float etac = eik->eta[i][j][k];
	  float Vnmoc2 = Vnmoc*Vnmoc;
	  float V0c2 = V0c*V0c;
	  float pzc = eik->pz[i][j][k];
	  float pxc = eik->px[i][j][k];
	  float pyc = eik->py[i][j][k];
	  eik->rhs[i][j][k] = 1+2*etac*Vnmoc2*V0c2*pzc*pzc*(pxc*pxc+pyc*pyc);
	  eik->T[i][j][k] = eik->T0[i][j][k]*eik->tau[i][j][k];
	}
      }
    }
    maxval = 0.0;
    for(i=0; i<n1; i++)
      for(j=0; j<n2; j++)
	for(k=0; k<n3; k++)
	  maxval = fmax(maxval, fabs(arr[i][j][k]-eik->T[i][j][k]));
    if(maxval<eik->epsilon) break;

  } 

  free3float(eik->T0);
  free3float(eik->px0);
  free3float(eik->py0);
  free3float(eik->pz0);
  free3float(eik->tau);
  free3float(eik->px);
  free3float(eik->py);
  free3float(eik->pz);
  free3float(eik->rhs);
  free3float(arr);
}
