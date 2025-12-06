#include "cstd.h"
#include "eikonal.h"

#define FAR_TIME 1000.0

void solveQuartic(double a, double b, double c, double d, double e, double* roots, int* num);
float shooting_method(float x0, float y0, float z0, float x1, float y1, float z1, float Vnmoc, float V0c, float eta);
void shooting_method_3D(float dx, float dy, float dz, float Vnmoc, float V0c, float etac, float* vel, float* pxc, float* pyc, float* pzc);

bool equals(float a, float b){
  return fabs(a-b) < 1e-12 ? true : false;
}

void Polynomial2rootsolver(float q1,float q0,double* root,int* nd)
{
  //x*x+q1*x+q0=0
  //root:only record the real root
  //nd:the number of real root
  float delta=q1*q1-4*q0;
  root[0]=0;
  root[1]=0;
  if(delta>0){
    root[0]=(-q1+sqrt(delta))/2.;
    root[1]=(-q1-sqrt(delta))/2.;
    *nd=2;
  }else if(fabs(delta)<1.0e-12){
    root[0]=-q1/2.;
    *nd=1;
  }else{
    *nd=0;
  }
}

void calculate_quartic_coefficients(float* a, float* b, float* c, float* d, float* e,
				    float a1, float b1, float a2, float b2, float a3, float b3, float Vnmoc, float V0c, float etac){
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c;
  *a = 2*etac*V0c2*Vnmoc2*a3*a3*(a1*a1+a2*a2);
  *b = 4*etac*V0c2*Vnmoc2*a3*((a1*a1+a2*a2)*b3+(a1*b1+a2*b2)*a3);
  *c = 2*etac*V0c2*Vnmoc2*((a1*a1+a2*a2)*b3*b3 + (b1*b1+b2*b2)*a3*a3 + 4*(a1*b1+a2*b2)*a3*b3)
    -(Vnmoc2*(1+2*etac)*(a1*a1+a2*a2)+V0c2*a3*a3);
  *d = 4*etac*V0c2*Vnmoc2*b3*((a1*b1+a2*b2)*b3 + (b1*b1+b2*b2)*a3)
    -2*(Vnmoc2*(1+2*etac)*(a1*b1+a2*b2)+V0c2*a3*b3);
  *e = 1 + 2*etac*V0c2*Vnmoc2*b3*b3*(b1*b1+b2*b2) - (Vnmoc2*(1+2*etac)*(b1*b1+b2*b2)+V0c2*b3*b3);
}

void get_velocity(eikonal_t *eik){
  for(int i=0; i<eik->n1; i++){
    for(int j=0; j<eik->n2; j++){
      for(int k=0; k<eik->n3; k++){
	float Vnmoc = eik->Vnmo[i][j][k], V0c = eik->V0[i][j][k], etac = eik->eta[i][j][k];
	eik->v_z[i][j][k] = V0c;
	eik->v_xy[i][j][k] = Vnmoc*sqrt(1+2*etac);
	eik->v_xz[i][j][k] = shooting_method(i*eik->h1, j*eik->h2, k*eik->h3, (i+1)*eik->h1, j*eik->h2, (k+1)*eik->h3, Vnmoc, V0c, etac);
	eik->v_yz[i][j][k] = shooting_method(i*eik->h1, j*eik->h2, k*eik->h3, i*eik->h1, (j+1)*eik->h2, (k+1)*eik->h3, Vnmoc, V0c, etac);
	eik->v_xyz[i][j][k] = shooting_method(i*eik->h1, j*eik->h2, k*eik->h3, (i+1)*eik->h1, (j+1)*eik->h2, (k+1)*eik->h3, Vnmoc, V0c, etac);
      }
    }
  }
}

void init_T0(eikonal_t *eik){
  int xc = (int)(eik->x_source/eik->h1);
  int yc = (int)(eik->y_source/eik->h2);
  int zc = (int)(eik->z_source/eik->h3);
  float Vnmoc = eik->Vnmo[xc][yc][zc], V0c = eik->V0[xc][yc][zc], etac = eik->eta[xc][yc][zc];
  float vel, dx, dy, dz, dist, pxc, pyc, pzc;
  for(int i=0; i<eik->n1; i++){
    for(int j=0; j<eik->n2; j++){
      for(int k=0; k<eik->n3; k++){
	dx = i*eik->h1 - eik->x_source;
	dy = j*eik->h2 - eik->y_source;
	dz = k*eik->h3 - eik->z_source;
	dist = sqrt(dx*dx+dy*dy+dz*dz);
	shooting_method_3D(dx,dy,dz,Vnmoc,V0c,etac,&vel,&pxc,&pyc,&pzc);
	eik->T0[i][j][k] = dist/vel;
	eik->px0[i][j][k] = pxc;
	eik->py0[i][j][k] = pyc;
	eik->pz0[i][j][k] = pzc;
      }
    }
  }
}

void init_tau(eikonal_t *eik){
  for(int i=0; i<eik->n1; i++)
    for(int j=0; j<eik->n2; j++)
      for(int k=0; k<eik->n3; k++)
	eik->tau[i][j][k] = FAR_TIME;

  int xc = (int)(eik->x_source/eik->h1);
  int yc = (int)(eik->y_source/eik->h2);
  int zc = (int)(eik->z_source/eik->h3);
  int num = 1;
  for(int i=xc-num+1; i<=xc+num; i++)
    for(int j=yc-num+1; j<=yc+num; j++)
      for(int k=zc-num+1; k<=zc+num; k++)
	eik->tau[i][j][k] = 1.0;
}


bool is_causality_3D(int im, int jm, int km, int ia, int ja, int ka, int ib, int jb, int kb, int ic, int jc, int kc, float taum, eikonal_t *eik)
{
  float T0a = eik->T0[ia][ja][ka], T0b = eik->T0[ib][jb][kb], T0c = eik->T0[ic][jc][kc], T0m = eik->T0[im][jm][km];
  float taua = eik->tau[ia][ja][ka], taub = eik->tau[ib][jb][kb], tauc = eik->tau[ic][jc][kc];
  float Ta = T0a * taua, Tb = T0b * taub, Tc = T0c * tauc, Tm = T0m * taum;
  float Vnmoc = eik->Vnmo[im][jm][km], V0c = eik->V0[im][jm][km], etac = eik->eta[im][jm][km];
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c;
  float h1 = eik->h1, h2 = eik->h2, h3 = eik->h3;
  float dx,dy,dz,px,py,pz,vx,vy,vz;
  if(Tm < Ta || Tm < Tb || Tm < Tc) return false;

  if(ia != im){
    dx = (im-ia)*h1;
    px = (Tm-Ta)/dx;
    if(ka != kb){
      dy = (jc-jb)*h2;
      dz = (kb-ka)*h3;
      py = (Tc-Tb)/dy;
      pz = (Tb-Ta)/dz;
      vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
      int flag = (vx*dx > 0) && (vy*dy < 0) && (vz*dz < 0);
      float a = h2*fabs(vz) - h3*fabs(vy);
      float b = h3*fabs(vx) - h1*fabs(vz);
      if(flag && a>0 && b>0) return true;
      else return false;
    }else if(ka == kb){
      dy = (jb-ja)*h2;
      dz = (kc-kb)*h3;
      py = (Tb-Ta)/dy;
      pz = (Tc-Tb)/dz;
      vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
      int flag = (vx*dx > 0) && (vy*dy < 0) && (vz*dz < 0);
      float a = h3*fabs(vy) - h2*fabs(vz);
      float b = h2*fabs(vx) - h1*fabs(vy);
      if(flag && a>0 && b>0) return true;
      else return false;
    }
  }else if(ja != jm){
    dy = (jm-ja)*h2;
    py = (Tm-Ta)/dy;
    if(ka != kb){
      dx = (ic-ib)*h1;
      dz = (kb-ka)*h3;
      px = (Tc-Tb)/dx;
      pz = (Tb-Ta)/dz;
      vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
      int flag = (vx*dx < 0) && (vy*dy > 0) && (vz*dz < 0);
      float a = h1*fabs(vz) - h3*fabs(vx);
      float b = h3*fabs(vy) - h2*fabs(vz);
      if(flag && a>0 && b>0) return true;
      else return false;
    }else if(ka == kb){
      dx = (ib-ia)*h1;
      dz = (kc-kb)*h3;
      px = (Tb-Ta)/dx;
      pz = (Tc-Tb)/dz;
      vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
      vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
      int flag = (vx*dx < 0) && (vy*dy > 0) && (vz*dz < 0);
      float a = h3*fabs(vx) - h1*fabs(vz);
      float b = h1*fabs(vy) - h2*fabs(vx);
      if(flag && a>0 && b>0) return true;
      else return false;
    }
  }else if(ka != km){
    dz = (km-ka)*h3;
    pz = (Tm-Ta)/dz;
    if(ja != jb){
      dx = (ic-ib)*h1;
      dy = (jb-ja)*h2;
      px = (Tc-Tb)/dx;
      py = (Tb-Ta)/dy;
      vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(py*py+pz*pz));
      vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(px*px+pz*pz));
      vz = 2*pz*V0c2*(1-2*Vnmoc2*etac*(px*px+py*py));
      int flag = (vx*dx < 0) && (vy*dy < 0) && (vz*dz > 0);
      float a = h1*fabs(vy) - h2*fabs(vx);
      float b = h2*fabs(vz) - h3*fabs(vy);
      if(flag && a>0 && b>0) return true;
      else return false;
    }else if(ja == jb){
      dx = (ib-ia)*h1;
      dy = (jc-jb)*h2;
      px = (Tb-Ta)/dx;
      py = (Tc-Tb)/dy;
      vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(py*py+pz*pz));
      vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(px*px+pz*pz));
      vz = 2*pz*V0c2*(1-2*Vnmoc2*etac*(px*px+py*py));
      int flag = (vx*dx < 0) && (vy*dy < 0) && (vz*dz > 0);
      float a = h2*fabs(vx) - h1*fabs(vy);
      float b = h1*fabs(vz) - h3*fabs(vx);
      if(flag && a>0 && b>0) return true;
      else return false;
    }
  }
    
}

float tri_solver_2D(int im, int jm, int km ,int ia, int ja, int ka, int ib, int jb, int kb, eikonal_t *eik)
{
  float h1 = eik->h1, h2 = eik->h2, h3 = eik->h3;
  float taua = eik->tau[ia][ja][ka], taub = eik->tau[ib][jb][kb], taum = eik->tau[im][jm][km];
  float T0a = eik->T0[ia][ja][ka], T0b = eik->T0[ib][jb][kb], T0m = eik->T0[im][jm][km];
  float px0m = eik->px0[im][jm][km], py0m = eik->py0[im][jm][km], pz0m = eik->pz0[im][jm][km];
  float Vnmoc = eik->Vnmo[im][jm][km], V0c = eik->V0[im][jm][km], etac = eik->eta[im][jm][km];
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c;
  float v_z = eik->v_z[im][jm][km], v_xy = eik->v_xy[im][jm][km], v_xz = eik->v_xz[im][jm][km], v_yz = eik->v_yz[im][jm][km], v_xyz = eik->v_xyz[im][jm][km];

  int ie, je, ke, ig, jg, kg, flag=1;
  float taux, tauy, tauz, px, py, pz;   // d(tau)/dx, dT/dx
  float taue, taug, T0e, T0g, dx, dy, dz, vx, vy, vz, result=FAR_TIME;   // vx,vy,vz : characteristic line direction

  double* root = alloc1double(4);
  int num, niter = 25;
  if(taua == FAR_TIME && taub == FAR_TIME) return taum;

  if(km == ka && km == kb){ // xy plane
    if(im != ia && jm != ja){
      ie = ia; je = ja; ke = ka;
      ig = ib; jg = jb; kg = kb;
    }else{
      ie = ib; je = jb; ke = kb;
      ig = ia; jg = ja; kg = ka;
    }
    taue = eik->tau[ie][je][ke]; taug = eik->tau[ig][jg][kg];
    T0e = eik->T0[ie][je][ke]; T0g = eik->T0[ig][jg][kg];
    if(taug == FAR_TIME) return fmin(taum, (taue * T0e + sqrt(h1*h1+h2*h2)/v_xy) / T0m);
    if(im == ig){
      if(taue == FAR_TIME) return fmin(taum, (taug * T0g + h2/v_xy) / T0m);
      dy = (jm-jg)*h2; dx = (ie-ig)*h1;
      taux = (taue-taug)/dx;
      float a1, b1, a2, b2, a3, b3;
      a1 = px0m; b1 = taux*T0m;
      a2 = py0m + T0m/dy; b2 = -T0m*taug/dy;
      a3 = 0; b3 = 0;
      float a,b,c,d,e;
      calculate_quartic_coefficients(&a,&b,&c,&d,&e,a1,b1,a2,b2,a3,b3,Vnmoc,V0c,etac);
      float q1 = d/c, q0 = e/c;
      Polynomial2rootsolver(q1,q0,root,&num);
      for(int i=0; i<num; i++){
	px = (taue*T0e - taug*T0g)/dx; py = (root[i]*T0m - taug*T0g)/dy;
	if(px*dx < 0 && py*dy > 0 && fabs(py)/h2 > fabs(px)/h1 && root[i]>0){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h2*h2)/v_xy) / T0m;
	float t2 = (T0g * taug + h2/v_xy) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }else if(jm == jg){
      if(taue == FAR_TIME) return fmin(taum, (taug * T0g + h1/v_xy) / T0m);
      dx = (im-ig)*h1; dy = (je-jg)*h2;
      tauy = (taue-taug)/dy;
      float a1, b1, a2, b2, a3, b3;
      a1 = px0m + T0m/dx; b1 = -T0m*taug/dx;
      a2 = py0m; b2 = tauy*T0m;
      a3 = 0; b3 = 0;
      float a,b,c,d,e;
      calculate_quartic_coefficients(&a,&b,&c,&d,&e,a1,b1,a2,b2,a3,b3,Vnmoc,V0c,etac);
      float q1 = d/c, q0 = e/c;
      Polynomial2rootsolver(q1,q0,root,&num);
      for(int i=0; i<num; i++){
	px = (root[i]*T0m - taug*T0g)/dx; py = (taue*T0e - taug*T0g)/dy;
	if(px*dx > 0 && py*dy < 0 && fabs(py)/h2 < fabs(px)/h1 && root[i]>0){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }               
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h2*h2)/v_xy) / T0m;
	float t2 = (T0g * taug + h1/(v_xy)) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }
  }else if(im == ia && im == ib){ // yz plane
    if(jm != ja && km != ka){
      ie = ia; je = ja; ke = ka;
      ig = ib; jg = jb; kg = kb;
    }else{
      ie = ib; je = jb; ke = kb;
      ig = ia; jg = ja; kg = ka;
    }
    T0e = eik->T0[ie][je][ke]; T0g = eik->T0[ig][jg][kg];
    taue = eik->tau[ie][je][ke]; taug = eik->tau[ig][jg][kg];
    if(taug == FAR_TIME) return fmin(taum, (taue * T0e + sqrt(h2*h2+h3*h3)/v_yz) / T0m);
    if(jm == jg){
      if(taue == FAR_TIME) return fmin(taum, (T0g * taug + h3/v_z) / T0m);
      dz = (km-kg)*h3; dy = (je-jg)*h2;
      tauy = (taue-taug)/dy;
      float a1, b1, a2, b2, a3, b3;
      a1 = 0; b1 = 0;
      a2 = py0m; b2 = tauy*T0m;
      a3 = pz0m + T0m/dz; b3 = -T0m*taug/dz;
      float a,b,c,d,e;
      calculate_quartic_coefficients(&a,&b,&c,&d,&e,a1,b1,a2,b2,a3,b3,Vnmoc,V0c,etac);
      solveQuartic(a,b,c,d,e,root,&num);
      for(int i=0; i<num; i++){
	pz = (root[i]*T0m - taug*T0g)/dz; py = (taue*T0e - taug*T0g)/dy;
	vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*py*py);
	if(py*dy < 0 && pz*dz > 0 && fabs(vz)/h3 > fabs(vy)/h2 && root[i]>0){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h2*h2+h3*h3)/v_yz) / T0m;
	float t2 = (T0g * taug + h3/v_z) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }else if(km == kg){
      if(taue == FAR_TIME) return fmin(taum, (taug * T0g + h2/v_xy) / T0m);
      dy = (jm-jg)*h2; dz = (ke-kg)*h3;
      tauz = (taue-taug)/dz;
      float a1, b1, a2, b2, a3, b3;
      a1 = 0; b1 = 0;
      a2 = py0m + T0m/dy; b2 = -T0m*taug/dy;
      a3 = pz0m; b3 = tauz*T0m;
      float a,b,c,d,e;
      calculate_quartic_coefficients(&a,&b,&c,&d,&e,a1,b1,a2,b2,a3,b3,Vnmoc,V0c,etac);
      solveQuartic(a,b,c,d,e,root,&num);
      for(int i=0; i<num; i++){
	pz = (taue*T0e - taug*T0g)/dz; py = (root[i]*T0m - taug*T0g)/dy;
	vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*py*py);
	if(py*dy > 0 && pz*dz < 0 && fabs(vz)/h3 < fabs(vy)/h2 && root[i]>0){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h2*h2+h3*h3)/v_yz) / T0m;
	float t2 = (T0g * taug + h2/v_xy) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }
  }else if(jm == ja && jm == jb){ // xz plane
    if(im != ia && km != ka){
      ie = ia; je = ja; ke = ka;
      ig = ib; jg = jb; kg = kb;
    }else{
      ie = ib; je = jb; ke = kb;
      ig = ia; jg = ja; kg = ka;
    }
    T0e = eik->T0[ie][je][ke]; T0g = eik->T0[ig][jg][kg];
    taue = eik->tau[ie][je][ke]; taug = eik->tau[ig][jg][kg];
    if(taug == FAR_TIME) return fmin(taum, (T0e * taue + sqrt(h1*h1+h3*h3)/v_xz) / T0m);
    if(im == ig){
      if(taue == FAR_TIME) return fmin(taum, (T0g * taug + h3/v_z) / T0m);
      dz = (km-kg)*h3; dx = (ie-ig)*h1;
      taux = (taue-taug)/dx;
      float a1, b1, a2, b2, a3, b3;
      a1 = px0m; b1 = taux*T0m;
      a2 = 0; b2 = 0;
      a3 = pz0m + T0m/dz; b3 = -T0m*taug/dz;
      float a,b,c,d,e;
      calculate_quartic_coefficients(&a,&b,&c,&d,&e,a1,b1,a2,b2,a3,b3,Vnmoc,V0c,etac);
      solveQuartic(a,b,c,d,e,root,&num);
      for(int i=0; i<num; i++){
	pz = (root[i]*T0m - taug*T0g)/dz; px = (taue*T0e - taug*T0g)/dx;
	vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*px*px);
	if(px*dx < 0 && pz*dz > 0 && fabs(vz)/h3 > fabs(vx)/h1 && root[i]>0){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h3*h3)/v_xz) / T0m;
	float t2 = (T0g * taug + h3/v_z) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }else if(km == kg){
      if(taue == FAR_TIME) return fmin(taum, (taug * T0g + h1/v_xy) / T0m);
      dx = (im-ig)*h1; dz = (ke-kg)*h3;
      tauz = (taue-taug)/dz;
      float a1, b1, a2, b2, a3, b3;
      a1 = px0m + T0m/dx; b1 = -T0m*taug/dx;
      a2 = 0; b2 = 0;
      a3 = pz0m; b3 = tauz*T0m;
      float a,b,c,d,e;
      calculate_quartic_coefficients(&a,&b,&c,&d,&e,a1,b1,a2,b2,a3,b3,Vnmoc,V0c,etac);
      solveQuartic(a,b,c,d,e,root,&num);
      for(int i=0; i<num; i++){
	pz = (taue*T0e - taug*T0g)/dz; px = (root[i]*T0m - taug*T0g)/dx;
	vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*px*px);
	if(px*dx > 0 && pz*dz < 0 && fabs(vz)/h3 < fabs(vx)/h1 && root[i]>0){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h3*h3)/v_xz) / T0m;
	float t2 = (T0g * taug + h1/v_xy) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }
  }
  else{   // tilting plane
    if(im != ia && jm != ja && km != ka){
      ie = ia; je = ja; ke = ka;
      ig = ib; jg = jb; kg = kb;
    }else if(im != ib && jm != jb && km != kb){
      ie = ib; je = jb; ke = kb;
      ig = ia; jg = ja; kg = ka;
    }
    T0e = eik->T0[ie][je][ke]; T0g = eik->T0[ig][jg][kg];
    taue = eik->tau[ie][je][ke]; taug = eik->tau[ig][jg][kg];
    if(taug == FAR_TIME) return fmin(taum, (taue * T0e + sqrt(h1*h1+h2*h2+h3*h3)/v_xyz) / T0m);
    if(im == ig && jm == jg && km != kg){  // case 1
      // causality : px*dx<0, py*dy<0, pz*dz>0 
      if(taue == FAR_TIME) return fmin(taum, (T0g * taug + h3/v_z) / T0m);
      dz = (km-kg)*h3; dx = (ie-ig)*h1; dy = (je-jg)*h2;
      float k = dy/dx;
      float a1 = (dx*px0m + dy*py0m)/(dx+k*dy), b1 = T0m*(taue-taug)/(dx+k*dy);
      float a2 = k*a1, b2 = k*b2;
      float a3 = pz0m + T0m/dz, b3 = -T0m*taug/dz;
      float a, b, c, d, e;
      calculate_quartic_coefficients(&a, &b, &c, &d, &e, a1, b1, a2, b2, a3, b3, Vnmoc, V0c, etac);
      solveQuartic(a, b, c, d, e, root, &num);
      for(int i=0; i<num; i++){
	px = a1*root[i]+b1;
	py = a2*root[i]+b2;
	pz = a3*root[i]+b3;
	// get root r, check causality
	vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
	float vr = sqrt(vx*vx+vy*vy);
	float Tm= T0m*root[i], Tg = T0g*taug, Te = T0e*taue;
	if((fabs(vz)/h3 - vr/sqrt(h1*h1+h2*h2) > 0) && (pz*dz > 0) && (px*dx <0) && (py*dy <0) && Tm>Tg && Tg>Te){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }                                                         
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h2*h2+h3*h3)/v_xyz) / T0m;
	float t2 = (T0g * taug + h3/v_z) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }
    else if(im != ig && jm != jg && km == kg){  // case 2
      // causality : px*dx>0, py*dy>0, pz*dz<0
      if(taue == FAR_TIME) return fmin(taum, (taug * T0g + sqrt(h1*h1+h2*h2)/v_xy) / T0m);
      dx = (im-ig)*h1; dy = (jm-jg)*h2; dz = (ke-kg)*h3;    
      float k = dy/dx;
      float a1 = (T0m + dx*px0m + dy*py0m)/(dx+k*dy), b1 = -T0m*taug/(dx+k*dy);
      float a2 = k*a1, b2 = k*b1;
      float a3 = pz0m, b3 = T0m*(taue-taug)/dz; 
      float a, b, c, d, e;
      calculate_quartic_coefficients(&a, &b, &c, &d, &e, a1, b1, a2, b2, a3, b3, Vnmoc, V0c, etac);
      solveQuartic(a, b, c, d, e, root, &num);
      for(int i=0; i<num; i++){
	px = a1*root[i]+b1;
	py = a2*root[i]+b2;
	pz = a3*root[i]+b3;
	// get root r, check causality
	vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz);
	vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
	float vr = sqrt(vx*vx+vy*vy);
	float Tm = T0m*root[i], Tg = T0g*taug, Te = T0e*taue;
	if((fabs(vz)/h3 - vr/sqrt(h1*h1+h2*h2) < 0) && (pz*dz < 0) && (px*dx >0) && (py*dy >0) && Tm>Tg && Tg>Te){
	  flag = 0;
	  result = fmin(root[i], result);
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h2*h2+h3*h3)/v_xyz) / T0m;
	float t2 = (T0g * taug + sqrt(h1*h1+h2*h2)/v_xy) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }
    else if(im == ig && jm != jg && km == kg){  // case 3
      // causality : px*dx<0, py*dy>0, pz*dz<0
      if(taue == FAR_TIME) return fmin(taum, (T0g * taug + h2/v_xy) / T0m);
      dy = (jm-jg)*h2; dz = (ke-kg)*h3; dx = (ie-ig)*h1; 
      float alpha = dx*px0m + dz*pz0m, beta;
      if(!equals(alpha,0)){
	float LEFT = 0, RIGHT = 1/V0c, MID, temp, r;
	// temp : pz is known, temp = px*px+py*py by eikonal equation
	// function f(pz) = temp - (px*px+py*py)
	float val_left, val_right, val_mid;
	// calculate f(LEFT)
	pz = 0; temp = 1./(Vnmoc2*(1+2*etac));
	vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	vx = dx*vz/dz;
	px = vx/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	beta = dx*px + dz*pz - T0m*(taue - taug);
	r = beta/alpha;
	py = (T0m*(r-taug) + dy*py0m*r)/dy;
	val_left = temp - (px*px + py*py);
	// calculate f(RIGHT)
	// obviously f(RIGHT) < 0
	val_right = -1;
	if(val_left*val_right < 0){ // there is a root
	  int iter = 0;
	  while(RIGHT - LEFT > 1e-12 && iter < niter){
	    MID = (LEFT + RIGHT)/2;
	    pz = dz>0? -MID : MID;
	    temp = (1-V0c2*pz*pz)/(Vnmoc2*(1+2*etac-2*etac*V0c2*pz*pz));
	    vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	    vx = dx*vz/dz;
	    px = vx/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	    beta = dx*px + dz*pz - T0m*(taue - taug);
	    r = beta/alpha;
	    py = (T0m*(r-taug) + dy*py0m*r)/dy;
	    val_mid = temp - (px*px + py*py);
	    if(equals(val_mid,0)) break;
	    else if(val_left*val_mid < 0){
	      RIGHT = MID;
	      val_right = val_mid;
	    }else{
	      LEFT = MID;
	      val_left = val_mid;
	    }
	    iter++;
	  }
	  // get root r, check causality
	  vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(px*px+pz*pz));
	  vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(py*py+pz*pz));
	  vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
	  float Tm= T0m*r, Tg = T0g*taug, Te = T0e*taue;
	  if((fabs(vy)/h2 - fabs(vz)/h3 > 0) && (pz*dz < 0) && (px*dx <0) && (py*dy >0)  && Tm>Tg && Tg>Te){
	    flag = 0;
	    result = fmin(r, result);
	  }
	}
      }                
      if(flag){
	float t1 = (T0e*taue + sqrt(h1*h1+h2*h2+h3*h3)/v_xyz) / T0m;
	float t2 = (T0g*taug + h2/v_xy) / T0m;
	result = fmin(result, fmin(t1, t2));
      }            
    }else if(im != ig && jm == jg && km != kg){   // case 4
      // causality : px*dx>0, py*dy<0, pz*dz>0
      if(taue == FAR_TIME) return fmin(taum, (T0g * taug + sqrt(h1*h1+h3*h3)/v_xz) / T0m);
      dx = (im-ig)*h1; dz = (km-kg)*h3; dy = (je-jg)*h2; 
      float alpha = T0m + dx*px0m + dz*pz0m, beta;
      if(!equals(alpha,0)){
	float LEFT = 0, RIGHT = 1/V0c, MID, temp, r;
	// function f(pz) = temp - (px*px+py*py)
	float val_left, val_right, val_mid;
	// calculate f(LEFT)
	pz = 0; temp = 1./(Vnmoc2*(1+2*etac));
	vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	vx = dx*vz/dz;
	px = vx/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	beta = dx*px + dz*pz + T0m*taug;
	r = beta/alpha;
	py = (T0m*(taue - taug) + dy*py0m*r)/dy;
	val_left = temp - (px*px + py*py);
	// calculate f(RIGHT)
	// obviously f(RIGHT) < 0
	val_right = -1;
	if(val_left*val_right < 0){ // there is a root
	  int iter = 0;
	  while(RIGHT - LEFT > 1e-12 && iter < niter){
	    MID = (LEFT + RIGHT)/2;
	    pz = dz>0? MID : -MID;
	    temp = (1-V0c2*pz*pz)/(Vnmoc2*(1+2*etac-2*etac*V0c2*pz*pz));
	    vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	    vx = dx*vz/dz;
	    px = vx/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	    beta = dx*px + dz*pz + T0m*taug;
	    r = beta/alpha;
	    py = (T0m*(taue - taug) + dy*py0m*r)/dy;
	    val_mid = temp - (px*px + py*py);
	    if(equals(val_mid,0)) break;
	    else if(val_left*val_mid < 0){
	      RIGHT = MID;
	      val_right = val_mid;
	    }else{
	      LEFT = MID;
	      val_left = val_mid;
	    }
	    iter++;
	  }
	  // get root r, check causality
	  vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(px*px+pz*pz));
	  vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(py*py+pz*pz));
	  vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
	  float Tm= T0m*r, Tg = T0g*taug, Te = T0e*taue;
	  if((fabs(vy)/h2 - fabs(vz)/h3 < 0) && (pz*dz > 0) && (px*dx >0) && (py*dy <0) && Tm>Tg && Tg>Te){
	    flag = 0;
	    result = fmin(r, result);
	  }
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h2*h2+h3*h3)/v_xyz) / T0m;
	float t2 = (T0g * taug + sqrt(h1*h1+h3*h3)/v_xz) / T0m;
	result = fmin(result, fmin(t1, t2));
      }              
    }else if(im != ig && jm == jg && km == kg){   // case 5
      // causality : px*dx>0, py*dy<0, pz*dz<0
      if(taue == FAR_TIME) return fmin(taum, (taug * T0g + h1/v_xy) / T0m);
      dx = (im-ig)*h1; dy = (je-jg)*h2; dz = (ke-kg)*h3;
      float alpha = dy*py0m + dz*pz0m, beta;
      if(!equals(alpha,0)){
	float LEFT = 0, RIGHT = 1/V0c, MID, temp, r;
	// function f(pz) = temp - (px*px+py*py)
	float val_left, val_right, val_mid;
	// calculate f(LEFT)
	pz = 0; temp = 1./(Vnmoc2*(1+2*etac));
	vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	vy = dy*vz/dz;
	py = vy/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	beta = dy*py + dz*pz - T0m*(taue - taug);
	r = beta/alpha;
	px = (T0m*(r - taug) + dx*px0m*r)/dx;
	val_left = temp - (px*px + py*py);
	// calculate f(RIGHT)
	// obviously f(RIGHT) < 0
	val_right = -1;
	if(val_left*val_right < 0){ // there is a root
	  int iter = 0;
	  while(RIGHT - LEFT > 1e-12 && iter < niter){
	    MID = (LEFT + RIGHT)/2;
	    pz = dz>0? -MID : MID;
	    temp = (1-V0c2*pz*pz)/(Vnmoc2*(1+2*etac-2*etac*V0c2*pz*pz));
	    vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	    vy = dy*vz/dz;
	    py = vy/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	    beta = dy*py + dz*pz - T0m*(taue - taug);
	    r = beta/alpha;
	    px = (T0m*(r - taug) + dx*px0m*r)/dx;
	    val_mid = temp - (px*px + py*py);
	    if(equals(val_mid,0)) break;
	    else if(val_left*val_mid < 0){
	      RIGHT = MID;
	      val_right = val_mid;
	    }else{
	      LEFT = MID;
	      val_left = val_mid;
	    }
	    iter++;
	  }
	  // get root r, check causality
	  vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(px*px+pz*pz));
	  vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(py*py+pz*pz));
	  vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
	  float Tm= T0m*r, Tg = T0g*taug, Te = T0e*taue;
	  if((fabs(vx)/h1 - fabs(vz)/h3 > 0) && (pz*dz < 0) && (px*dx >0) && (py*dy <0) && Tm>Tg && Tg>Te){
	    flag = 0;
	    result = fmin(r, result);
	  }
	}
      }
      if(flag){
	float t1 = (T0e*taue + sqrt(h1*h1+h2*h2+h3*h3)/v_xyz) / T0m;
	float t2 = (T0g*taug + h1/v_xy) / T0m;
	result = fmin(result, fmin(t1, t2));
      }          
    }else if(im == ig && jm != jg && km != kg){   // case 6
      // causality : px*dx<0, py*dy>0, pz*dz>0
      if(taue == FAR_TIME) return fmin(taum, (taug * T0g + sqrt(h2*h2+h3*h3)/v_yz) / T0m);
      dy = (jm-jg)*h2; dx = (ie-ig)*h1; dz = (km-kg)*h3;
      float alpha = T0m + dy*py0m + dz*pz0m, beta;
      if(!equals(alpha,0)){
	float LEFT = 0, RIGHT = 1/V0c, MID, temp, r;
	// function f(pz) = temp - (px*px+py*py)
	float val_left, val_right, val_mid;
	// calculate f(LEFT)
	pz = 0; temp = 1./(Vnmoc2*(1+2*etac));
	vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	vy = dy*vz/dz;
	py = vy/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	beta = dy*py + dz*pz + T0m*taug;
	r = beta/alpha;
	px = (T0m*(taue - taug) + dx*px0m*r)/dx;
	val_left = temp - (px*px + py*py);
	// calculate f(RIGHT)
	// obviously f(RIGHT) < 0
	val_right = -1;
	if(val_left*val_right < 0){ // there is a root
	  int iter = 0;
	  while(RIGHT - LEFT > 1e-12 && iter < niter){
	    MID = (LEFT + RIGHT)/2;
	    pz = dz>0? MID : -MID;
	    temp = (1-V0c2*pz*pz)/(Vnmoc2*(1+2*etac-2*etac*V0c2*pz*pz));
	    vz = 2*V0c2*pz*(1 - 2*Vnmoc2*etac*temp);
	    vy = dy*vz/dz;
	    py = vy/(2*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*pz*pz));
	    beta = dy*py + dz*pz + T0m*taug;
	    r = beta/alpha;
	    px = (T0m*(taue - taug) + dx*px0m*r)/dx;
	    val_mid = temp - (px*px + py*py);
	    if(equals(val_mid,0)) break;
	    else if(val_left*val_mid < 0){
	      RIGHT = MID;
	      val_right = val_mid;
	    }else{
	      LEFT = MID;
	      val_left = val_mid;
	    }
	    iter++;
	  }
	  // get root r, check causality
	  vy = 2*py*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(px*px+pz*pz));
	  vx = 2*px*(Vnmoc2*(1+2*etac)-2*V0c2*Vnmoc2*etac*(py*py+pz*pz));
	  vz = 2*V0c2*pz*(1-2*Vnmoc2*etac*(px*px+py*py));
	  float Tm= T0m*r, Tg = T0g*taug, Te = T0e*taue;
	  if((fabs(vx)/h1 - fabs(vz)/h3 < 0) && (pz*dz > 0) && (px*dx <0) && (py*dy >0) && Tm>Tg && Tg>Te){
	    flag = 0;
	    result = fmin(r, result);
	  }
	}
      }
      if(flag){
	float t1 = (T0e * taue + sqrt(h1*h1+h2*h2+h3*h3)/v_xyz) / T0m;
	float t2 = (T0g * taug + sqrt(h2*h2+h3*h3)/v_yz) / T0m;
	result = fmin(result, fmin(t1, t2));
      }
    }     
  }  
    
  free1double(root);
  return fmin(result, taum);      
}
    
float tri_solver_3D(int im, int jm, int km, int ia, int ja, int ka, int ib, int jb, int kb, int ic, int jc, int kc, eikonal_t *eik)
{
  float h1 = eik->h1, h2 = eik->h2, h3 = eik->h3;
  float taua = eik->tau[ia][ja][ka], taub = eik->tau[ib][jb][kb], tauc = eik->tau[ic][jc][kc], taum = eik->tau[im][jm][km];
  float T0a = eik->T0[ia][ja][ka], T0b = eik->T0[ib][jb][kb], T0c = eik->T0[ic][jc][kc], T0m = eik->T0[im][jm][km];
  float px0m = eik->px0[im][jm][km], py0m = eik->py0[im][jm][km], pz0m = eik->pz0[im][jm][km];
  float taux, tauy, tauz, px, py, pz;          // taux -> d(tau)/dx , px -> dT/dx
  float Vnmoc = eik->Vnmo[im][jm][km], V0c = eik->V0[im][jm][km], etac = eik->eta[im][jm][km];
  float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c;

  if(taua == FAR_TIME && taub == FAR_TIME && tauc == FAR_TIME) return taum;
  if(taua == FAR_TIME) return fmin(taum, tri_solver_2D(im, jm, km, ib, jb, kb, ic, jc, kc, eik));
  if(taub == FAR_TIME) return fmin(taum, tri_solver_2D(im, jm, km, ia, ja, ka, ic, jc, kc, eik));
  if(tauc == FAR_TIME) return fmin(taum, tri_solver_2D(im, jm, km, ia, ja, ka, ib, jb, kb, eik));

  double* root = alloc1double(4);
  float dx, dy, dz;
  float result = FAR_TIME;
  int flag = 1, num_root = 0;

  if(ia != im){             // case 1,2
    dx = (im-ia)*h1;
    if(ka != kb){         // case 1
      dy = (jc-jb)*h2;
      dz = (kb-ka)*h3;
      tauy = (tauc-taub)/dy;
      tauz = (taub-taua)/dz;
    }else if(ka == kb){   // case 2
      dy = (jb-ja)*h2;
      dz = (kc-kb)*h3;
      tauy = (taub-taua)/dy;
      tauz = (tauc-taub)/dz;
    }
    float a1 = T0m/dx + px0m, b1 = -T0m*taua/dx;
    float a2 = py0m, b2 = T0m*tauy;
    float a3 = pz0m, b3 = T0m*tauz;
    float a, b, c, d, e;
    calculate_quartic_coefficients(&a, &b, &c, &d, &e, a1, b1, a2, b2, a3, b3, Vnmoc, V0c, etac);
    solveQuartic(a, b, c, d, e, root, &num_root);
    for(int i=0; i<num_root; i++){
      if(is_causality_3D(im, jm, km, ia, ja, ka, ib, jb, kb, ic, jc, kc, root[i], eik)){
	result = fmin(root[i], result);
	flag = 0;
      }
    }
  }else if(ja != jm){       // case 3,4
    dy = (jm-ja)*h2;
    if(ka != kb){         // case 3
      dx = (ic-ib)*h1;
      dz = (kb-ka)*h3;
      taux = (tauc-taub)/dx;
      tauz = (taub-taua)/dz;
    }else if(ka == kb){   // case 4
      dx = (ib-ia)*h1;
      dz = (kc-kb)*h3;
      taux = (taub-taua)/dx;
      tauz = (tauc-taub)/dz;
    }
    float a1 = px0m, b1 = T0m*taux;
    float a2 = T0m/dy + py0m, b2 = -taua*T0m/dy;
    float a3 = pz0m, b3 = T0m*tauz;
    float a, b, c, d, e;
    calculate_quartic_coefficients(&a, &b, &c, &d, &e, a1, b1, a2, b2, a3, b3, Vnmoc, V0c, etac);
    solveQuartic(a, b, c, d, e, root, &num_root);
    for(int i=0; i<num_root; i++){
      if(is_causality_3D(im, jm, km, ia, ja, ka, ib, jb, kb, ic, jc, kc, root[i], eik)){
	result = fmin(root[i], result);
	flag = 0;
      }
    }
  }else if(ka != km){ // case 5,6
    dz = (km-ka)*h3;
    if(ja != jb){         // case 5
      dx = (ic-ib)*h1;
      dy = (jb-ja)*h2;
      taux = (tauc-taub)/dx;
      tauy = (taub-taua)/dy;
    }else if(ja == jb){   // case 6
      dx = (ib-ia)*h1;
      dy = (jc-jb)*h2;
      taux = (taub-taua)/dx;
      tauy = (tauc-taub)/dy;
    }
    float a1 = px0m, b1 = T0m*taux;
    float a2 = py0m, b2 = T0m*tauy;
    float a3 = T0m/dz + pz0m, b3 = -T0m*taua/dz;
    float a, b, c, d, e;
    calculate_quartic_coefficients(&a, &b, &c, &d, &e, a1, b1, a2, b2, a3, b3, Vnmoc, V0c, etac);
    solveQuartic(a, b, c, d, e, root, &num_root);
    for(int i=0; i<num_root; i++){
      if(is_causality_3D(im, jm, km, ia, ja, ka, ib, jb, kb, ic, jc, kc, root[i], eik)){
	result = fmin(root[i], result);
	flag = 0;
      }
    }
  }
  if(flag){
    float t1 = tri_solver_2D(im, jm, km, ia, ja, ka, ib, jb, kb, eik);
    float t2 = tri_solver_2D(im, jm, km, ib, jb, kb, ic, jc, kc, eik);
    float t3 = tri_solver_2D(im, jm, km, ia, ja, ka, ic, jc, kc, eik);
    result = fmin(result, fmin(t1, fmin(t2, t3)));
  }
  free1double(root);   
  return fmin(result, taum);
}

void neighbors(int** arr){
  /*  arr[][0,1,2] -> ia,ja,ka;
      arr[][3,4,5] -> ib,jb,kb;
      arr[][6,7,8] -> ic,jc,kc; */
  // only record the division corresponding to the first octant
  arr[0][0] = 1; arr[0][1] = 0; arr[0][2] = 0;
  arr[0][3] = 1; arr[0][4] = 0; arr[0][5] = 1;
  arr[0][6] = 1; arr[0][7] = 1; arr[0][8] = 1;

  arr[1][0] = 1; arr[1][1] = 0; arr[1][2] = 0;
  arr[1][3] = 1; arr[1][4] = 1; arr[1][5] = 0;
  arr[1][6] = 1; arr[1][7] = 1; arr[1][8] = 1;

  arr[2][0] = 0; arr[2][1] = 1; arr[2][2] = 0;
  arr[2][3] = 0; arr[2][4] = 1; arr[2][5] = 1;
  arr[2][6] = 1; arr[2][7] = 1; arr[2][8] = 1;

  arr[3][0] = 0; arr[3][1] = 1; arr[3][2] = 0;
  arr[3][3] = 1; arr[3][4] = 1; arr[3][5] = 0;
  arr[3][6] = 1; arr[3][7] = 1; arr[3][8] = 1;

  arr[4][0] = 0; arr[4][1] = 0; arr[4][2] = 1;
  arr[4][3] = 0; arr[4][4] = 1; arr[4][5] = 1;
  arr[4][6] = 1; arr[4][7] = 1; arr[4][8] = 1;

  arr[5][0] = 0; arr[5][1] = 0; arr[5][2] = 1;
  arr[5][3] = 1; arr[5][4] = 0; arr[5][5] = 1;
  arr[5][6] = 1; arr[5][7] = 1; arr[5][8] = 1;
}

// determine upwind stencil in x,y,z directions
void determine_upwind_stencil(eikonal_t *eik, int xc, int yc, int zc, int* sx, int* sy, int* sz){
  if(xc == 0) *sx = 1;
  else if(xc == eik->n1-1) *sx = -1;
  else{
    if(eik->T0[xc-1][yc][zc]*eik->tau[xc-1][yc][zc] < eik->T0[xc+1][yc][zc]*eik->tau[xc+1][yc][zc]) *sx = -1;
    else *sx = 1;
  }

  if(yc == 0) *sy = 1;
  else if(yc == eik->n2-1) *sy = -1;
  else{
    if(eik->T0[xc][yc-1][zc]*eik->tau[xc][yc-1][zc] < eik->T0[xc][yc+1][zc]*eik->tau[xc][yc+1][zc]) *sy = -1;
    else *sy = 1;
  }

  if(zc == 0) *sz = 1;
  else if(zc == eik->n3-1) *sz = -1;
  else{
    if(eik->T0[xc][yc][zc-1]*eik->tau[xc][yc][zc-1] < eik->T0[xc][yc][zc+1]*eik->tau[xc][yc][zc+1]) *sz = -1;
    else *sz = 1;
  }
}

void local_solver(int im, int jm, int km, eikonal_t *eik, int** arr)
{
  int xs = (int)(eik->x_source/eik->h1);
  int ys = (int)(eik->y_source/eik->h2);
  int zs = (int)(eik->z_source/eik->h3);
  if((im == xs || im == xs+1) && (jm == ys || jm == ys+1) && (km == zs || km == zs+1)) return;
  int sx,sy,sz,ia,ja,ka,ib,jb,kb,ic,jc,kc;
  float result = FAR_TIME;
  determine_upwind_stencil(eik, im, jm, km, &sx, &sy, &sz);
    
  for(int n=0; n<6; n++){
    ia = im + sx*arr[n][0];
    ja = jm + sy*arr[n][1];
    ka = km + sz*arr[n][2];
    ib = im + sx*arr[n][3];
    jb = jm + sy*arr[n][4];
    kb = km + sz*arr[n][5];
    ic = im + sx*arr[n][6];
    jc = jm + sy*arr[n][7];
    kc = km + sz*arr[n][8];
    float temp = tri_solver_3D(im, jm, km, ia, ja, ka, ib, jb, kb, ic, jc, kc, eik);
    result = fmin(result, temp);
  }
  eik->tau[im][jm][km] = fmin(eik->tau[im][jm][km], result);
}

void fast_sweeping(eikonal_t* eik, int** arr)
{
  int i, j, k;
  int n1 = eik->n1;
  int n2 = eik->n2;
  int n3 = eik->n3;
  float maxval;
  float*** tau_old = alloc3float(n3,n2,n1);
  for(int iter=0; iter<eik->niter; iter++)
    {
      printf("# %d-th loop\n",iter+1);
      memcpy(&tau_old[0][0][0], &eik->tau[0][0][0], n3*n2*n1*sizeof(float));

      for(i=0; i<n1; i++)
	for(j=0; j<n2; j++)
	  for(k=0; k<n3; k++)
	    local_solver(i, j, k, eik, arr);
       
      for(i=n1-1; i>=0; i--)
	for(j=0; j<n2; j++)
	  for(k=0; k<n3; k++)
	    local_solver(i, j, k, eik, arr);
                
      for(i=0; i<n1; i++)
	for(j=n2-1; j>=0; j--)
	  for(k=0; k<n3; k++)
	    local_solver(i, j, k, eik, arr);
               
      for(i=n1-1; i>=0; i--)
	for(j=n2-1; j>=0; j--)
	  for(k=0; k<n3; k++)
	    local_solver(i, j, k, eik, arr);
               
      for(i=0; i<n1; i++)
	for(j=0; j<n2; j++)
	  for(k=n3-1; k>=0; k--)
	    local_solver(i, j, k, eik, arr);
               
      for(i=n1-1; i>=0; i--)
	for(j=0; j<n2; j++)
	  for(k=n3-1; k>=0; k--)
	    local_solver(i, j, k, eik, arr);
                
      for(i=0; i<n1; i++)
	for(j=n2-1; j>=0; j--)
	  for(k=n3-1; k>=0; k--)
	    local_solver(i, j, k, eik, arr);
               
      for(i=n1-1; i>=0; i--)
	for(j=n2-1; j>=0; j--)
	  for(k=n3-1; k>=0; k--)
	    local_solver(i, j, k, eik, arr);
        
      maxval = 0;
      for(i=0; i<n1; i++)
	for(j=0; j<n2; j++)
	  for(k=0; k<n3; k++)
	    maxval = fmax(maxval, fabs(tau_old[i][j][k]-eik->tau[i][j][k]));
      if(maxval < eik->epsilon) break;
    }
  free3float(tau_old);
}

void eikonal_solver(eikonal_t* eik)
{
  int n1 = eik->n1;
  int n2 = eik->n2;
  int n3 = eik->n3;
  eik->v_z = alloc3float(n3,n2,n1);
  eik->v_xy = alloc3float(n3,n2,n1);
  eik->v_xz = alloc3float(n3,n2,n1);
  eik->v_yz = alloc3float(n3,n2,n1);
  eik->v_xyz = alloc3float(n3,n2,n1);
  get_velocity(eik);

  int** arr = alloc2int(9,6);
  neighbors(arr);

  eik->T0 = alloc3float(n3,n2,n1);
  eik->px0 = alloc3float(n3,n2,n1);
  eik->py0 = alloc3float(n3,n2,n1);
  eik->pz0 = alloc3float(n3,n2,n1);
  eik->tau = alloc3float(n3,n2,n1);
  init_T0(eik);
  init_tau(eik);

  fast_sweeping(eik, arr);

  for(int i=0; i<n1; i++)
    for(int j=0; j<n2; j++)
      for(int k=0; k<n3; k++)
	eik->T[i][j][k] = eik->tau[i][j][k] * eik->T0[i][j][k];

  free3float(eik->T0);
  free3float(eik->px0);
  free3float(eik->py0);
  free3float(eik->pz0);
  free3float(eik->tau);
  free3float(eik->v_z);
  free3float(eik->v_xy);
  free3float(eik->v_xz);
  free3float(eik->v_yz);
  free3float(eik->v_xyz);
  free2int(arr);
}


    
