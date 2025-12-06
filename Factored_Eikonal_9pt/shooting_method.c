#include <math.h>

void shooting_method_2D(float dx, float dz, float Vnmoc, float V0c, float etac, float* vel, float* pxc, float* pzc){
	float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c;
	float sx = -1, sz = -1;
	if(dx>0) sx = 1;
	if(dz>0) sz = 1;
	if(fabs(dz)<1.0e-12){
		*vel = Vnmoc*sqrt(1+2*etac); *pxc = sx*1.0/(*vel); *pzc = 0;
	}else if(fabs(dx)<1.0e-12){
		*vel = V0c; *pxc = 0; *pzc = sx*1.0/V0c;
	}else{
		float tan = fabs(dx/dz);
		float pa = 0, pb = 1/sqrt(Vnmoc2*(1+2*etac)), pc = (pa+pb)/2.;
		float pz, vx, vz, v;
		while(fabs(pb-pa)>=1.0e-10){
			pz = sqrt((1-Vnmoc2*(1+2*etac)*pc*pc)/(V0c2*(1-2*etac*Vnmoc2*pc*pc)));
			float alpha = Vnmoc2*(1+2*etac-2*etac*V0c2*pz*pz);
			float beta = V0c2*(1-2*etac*Vnmoc2*pc*pc);
			float q = 2*pc*pc*alpha+2*pz*pz*beta;
			vx = 2*pc*alpha/q;
			vz = 2*pz*beta/q;
			v = sqrt(vx*vx+vz*vz);
			float c = vx/vz;
			if(fabs(c-tan)<1.0e-15) break;
			else if(c > tan) pb = pc;
			else pa = pc;
			pc = (pa+pb)/2.0;
		}
		*vel = v; *pxc = sx*pc; *pzc = sz*pz;
	}
}

// dx : x_target_point - x_source_point
void shooting_method_3D(float dx, float dy, float dz, float Vnmoc, float V0c, float etac, float* vel, float* pxc, float* pyc, float* pzc){
	float Vnmoc2 = Vnmoc*Vnmoc, V0c2 = V0c*V0c;
	float sx = -1, sy = -1, sz = -1;
	if(dx>0) sx = 1;
	if(dy>0) sy = 1;
	if(dz>0) sz = 1;
	float r = sqrt(dx*dx+dy*dy);
	float cos = fabs(dx)/r, sin = fabs(dy)/r;
	if(fabs(dz)<1.0e-12){
		*vel = Vnmoc*sqrt(1+2*etac); *pxc = sx*cos/(*vel); *pyc = sy*sin/(*vel); *pzc = 0;
	}else if(r<1.0e-12){
		*vel = V0c; *pxc = 0; *pyc = 0; *pzc = sz*1.0/V0c;
	}else{
		float tan = fabs(r/dz);
		float pa = 0, pb = 1/sqrt(Vnmoc2*(1+2*etac)), pc = (pa+pb)/2.;
		float pz, vr, vz, v;
		while(fabs(pb-pa)>=1.0e-10){
			pz = sqrt((1-Vnmoc2*(1+2*etac)*pc*pc)/(V0c2*(1-2*etac*Vnmoc2*pc*pc)));
			float alpha = Vnmoc2*(1+2*etac-2*etac*V0c2*pz*pz);
			float beta = V0c2*(1-2*etac*Vnmoc2*pc*pc);
			float q = 2*pc*pc*alpha+2*pz*pz*beta;
			vr = 2*pc*alpha/q;
			vz = 2*pz*beta/q;
			v = sqrt(vr*vr+vz*vz);
			float c = vr/vz;
			if(fabs(c-tan)<1.0e-15) break;
			else if(c > tan) pb = pc;
			else pa = pc;
			pc = (pa+pb)/2.0;
		}
		*vel = v; *pxc = sx*cos*pc; *pyc = sy*sin*pc; *pzc = sz*pz;
	}
}

// return velocity
float shooting_method(float x0, float y0, float z0, float x1, float y1, float z1, float Vnmoc, float V0c, float eta){
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
		if(fabs(c-tan)<1.0e-15){
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
