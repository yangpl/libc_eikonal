#include <math.h>

float shooting_method(float x0, float y0, float z0, float x1, float y1, float z1, float Vnmoc, float V0c, float eta)
{
  float V_nmo = Vnmoc*Vnmoc, V0 = V0c*V0c;
  float dx = fabs(x0-x1);
  float dy = fabs(y0-y1);
  float dz = fabs(z0-z1);
  float r = sqrt(dx*dx+dy*dy);
  if(dz<1.0e-10)
    return sqrt((1+2*eta)*V_nmo);
  float tan = r/dz;
  float pa,pb,pc,pz;
  float alpha,beta;
  float vr,vz,v;
  pa=0; pb=1/sqrt((1+2*eta)*V_nmo); pc=(pa+pb)/2.;
  int iter = 0;
  while(fabs(pb-pa)>=1.0e-10){
    iter++;
    pz = (1-V_nmo*(1+2*eta)*pc*pc)/(V0*(1-2*eta*V_nmo*pc*pc));
    pz = sqrt(pz);
    alpha = V_nmo*(1+2*eta)-2*V0*V_nmo*eta*pz*pz;
    beta = V0*(1-2*V_nmo*eta*pc*pc);
    float q = 2*pc*pc*alpha+2*pz*pz*beta;
    vr = 2*pc*alpha/q;
    vz = 2*pz*beta/q;
    v = sqrt(vr*vr+vz*vz);
    float c = vr/vz;
    if(fabs(c-tan)<1.0e-10){
      return v;
    }else if(c>tan){
      pb = pc;
    }else{
      pa = pc;
    }
    pc = (pa+pb)/2.;
  }
  return v;
}
