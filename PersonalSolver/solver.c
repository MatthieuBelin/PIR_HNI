// #include "grid/linear.h"

void flux (double h, double u, double &flux1, double &flux2) {
	flux1 = h*u;
	flux2 = h*u*u + G*h*h/2.;
}

void HLL (double hL, double uL, double hR, double uR, double &F1, double &F2, double &vmax) {
	
	// wave speeds
	double sL = min(uL - sqrt(G*hL), uR - sqrt(G*hR));
	double sR = max(uL + sqrt(G*hL), uR + sqrt(G*hR));
	vmax = max(abs(sL), abs(sR));
	
	// fluxes
	double f1L, f2L;
	flux(hL, uL, f1L, f2L);
	double f1R, f2R;
	flux(hR, uR, f1R, f2R);
	
	if (0. < sL) {
		F1 = f1L;
		F2 = f2L;
	}
	else if (sR < 0.) {
		F1 = f1R;
		F2 = f2R;
	}
	else {
	        F1 = (sR*f1L - sL*f1R)/(sR-sL) + sL*sR*(hR - hL)/(sR-sL);
        	F2 = (sR*f2L - sL*f2R)/(sR-sL) + sL*sR*(hR*uR - hL*uL)/(sR-sL);
        }	
}
