#include "grid/quadtree.h"
#include "paraview2d.h"
#include "run.h"

#define TMAX 1.

//==============================Constants==============================

scalar h[]; // Height of the water
scalar zb[]; // Bathymetrie
vector u[]; // Speed in x and y direction
face vector fh[];
tensor fu[];

double g = 1.;

double DX;
double dry = 1.e-10;
int level = 6;

//==============================boundary conditions==============================
u.n[left] = neumann(0.); u.n[right] = neumann(0.);
u.n[bottom] = neumann(0.); u.n[top] = neumann(0.);

//==============================Functions for the solver==============================

void flux (double h1, double u1, double *flux1, double *flux2) {
	*flux1 = h1*u1;
	*flux2 = h1*u1*u1 + g*h1*h1/2.;
}

void HLL (double hL, double uL, double hR, double uR, double *F1, double *F2, double *vmax) {
	
	// wave speeds
	double sL = min(uL - sqrt(g*hL), uR - sqrt(g*hR));
	double sR = max(uL + sqrt(g*hL), uR + sqrt(g*hR));
	*vmax = max(abs(sL), abs(sR));
	
	// fluxes
	double f1L, f2L;
	flux(hL, uL, &f1L, &f2L);
	double f1R, f2R;
	flux(hR, uR, &f1R, &f2R);
	
	if (0. < sL) {
		*F1 = f1L;
		*F2 = f2L;
	}
	else if (sR < 0.) {
		*F1 = f1R;
		*F2 = f2R;
	}
	else {
	        *F1 = (sR*f1L - sL*f1R)/(sR-sL) + sL*sR*(hR - hL)/(sR-sL);
        	*F2 = (sR*f2L - sL*f2R)/(sR-sL) + sL*sR*(hR*uR - hL*uL)/(sR-sL);
        }	
}

//==============================Events to solve==============================

event init(t = 0) {
	// Initial conditions

	// Riemann problem
	double dist;	
	foreach() {
		dist = sqrt(x*x + y*y);
		h[] = dist < 1. ? 4. : 1.;
		zb[] = 0.;
		u.x[] = 0.;
		u.y[] = 0.;
	}
	
	DX = (L0/N); // cartesian grid
}

event plot (t += 0.1) {
	fprintf (stdout, "# t = %g\n", t); // time
	// paraview
	output_paraview (slist = {h}, vlist = {u});
}

event solve (i++) {	
	// Compute the solution for this iteration
	
	// Neumann boundary condition
	boundary ({h, u});
	
	// Fluxes at interfaces
	face vector Fh[];	// Fh.x,	Fh.y
	tensor Fhu[];		// Fhu.x.x	Fhu.x.y
				// Fhu.y.x	Fhu.y.y
	
	// Fluxes
	double vmax = 0.;
	// Source term
	face vector s[];
	double zminushalf, zplushalf, hminushalf, hplushalf; // z[i-1/2], z[i+1/2], h[i-1/2+], h[i+1/2-]
	foreach_face() {
		// Face in x AND y directions
		double hL = h[-1], uL = u.x[-1], vL = u.y[-1]; // left value
		double hR = h[],   uR = u.x[],   vR = u.y[];   // right value
		if (hL > dry || hR > dry) {
			// F(W_L,W_R) = (F1, F2)
			double F1, F2, vmax_loc;
			// numerical flux
			HLL (hL, uL, hR, uR, &F1, &F2, &vmax_loc);
			// get fluxes
			Fh.x[] = F1;
			Fhu.x.x[] = F2;
			Fhu.y.x[] = (F1 > 0. ? F1*vL : F1*vR);  // HLLC
			// Source term
            zminushalf = max(zb[-1], zb[]);
            zplushalf = max(zb[], zb[+1]);
            hminushalf = max(0, h[] + zb[] - zminushalf);
            hplushalf = max(0, h[] + zb[] - zplushalf);
            s.x[] = g/2*(hplushalf*hplushalf - hminushalf*hminushalf);

			vmax = max(vmax, vmax_loc);
		} else {
			Fh.x[] = Fhu.x.x[] = Fhu.y.x[] = 0.;
		}
	}
	
	// Timestep under CFL
	DT = CFL*DX/fabs(vmax);
	
	// next timestep
	double dt = dtnext(DT);
	
	// Update the solution
	foreach() {
		// Fh.x,    Fh.y
		// Fhu.x.x, Fhu.x.y
		// Fhu.y.x, Fhu.y.y

		// update h
		double h_old = h[]; // backup h at t = t^n 
		// note that: FL = F[], FR = F[1]
		double FhL = Fh.x[]  + Fh.y[]; // left + bottom
		double FhR = Fh.x[1,0] + Fh.y[0,1]; // right + top 
		h[] += - (dt/DX)*(FhR - FhL); // h at t = t^{n+1}
    
		// update u,v
		if (h[] > dry) {
			foreach_dimension() { // u, v
				double FhuL = Fhu.x.x[]  + Fhu.x.y[]; // left + bottom
				double FhuR = Fhu.x.x[1,0] + Fhu.x.y[0,1]; // right + top 
				double hu = h_old*u.x[] - (dt/DX)*(FhuR - FhuL + s.x[]);
				u.x[] = hu/h[];
			}
		} else { // dry ==> set u = v = 0
			u.x[] = u.y[] = 0.;
		}
	}
}

event end (t = TMAX) {
	printf ("i = %d t = %g\n", i, t);
	return 1;
}

//event adapt (i++) {
//	adapt_wavelet ({h}, (double []){4e-3}, maxlevel = level);
//}

//==============================Main function==============================

int main() {
	// resolution
	// if (argc > 1) level = atoi(argv[1]);

	// domain
	g = 1.; 
	L0 = 6.;   
	X0 = -L0/2.; Y0 = -L0/2.;
	N = 1 << level; // 2^8 = 256
	CFL = 0.4;
	run();
}
