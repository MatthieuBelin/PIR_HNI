// qcc -O2 -Wall solver.c -o solver -lm

#include "grid/quadtree.h"
#include "paraview2d.h"
#include "run.h"

//double TMAX = 5.;

//==============================Constants==============================

scalar h[]; // Height of the water
scalar zb[]; // Bathymetrie
scalar zs[]; // Free surface
vector u[]; // Speed in x and y direction
face vector fh[];
tensor fu[];

double g = 9.81;

double h0, a, r0;
double omega;
double TMAX;

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

void HLL (double hG, double uL, double zL, double hD, double uR, double zR, double *F1, double *F2L, double *F2R, double *vmax) {

    double zi = max(zL, zR);
    double hL = max(0, hG + zL - zi);
    double hR = max(0, hD + zR - zi);

	// wave speeds
	double sL = min(uL - sqrt(g*hL), uR - sqrt(g*hR));
	double sR = max(uL + sqrt(g*hL), uR + sqrt(g*hR));
	*vmax = max(fabs(sL), fabs(sR));

    if (max(hR, hL) < dry) {
        *F1 = 0.;
        *F2L = 0.;
        *F2R = 0.;
        return;
    }

	// fluxes
	double f1L, f2L;
	flux(hL, uL, &f1L, &f2L);
	double f1R, f2R;
	flux(hR, uR, &f1R, &f2R);

	double F2;
	if (0. < sL) {
		*F1 = f1L;
		F2 = f2L;
	}
	else if (sR < 0.) {
		*F1 = f1R;
		F2 = f2R;
	}
	else {
	    *F1 = (sR*f1L - sL*f1R)/(sR-sL) + sL*sR*(hR - hL)/(sR-sL);
        F2 = (sR*f2L - sL*f2R)/(sR-sL) + sL*sR*(hR*uR - hL*uL)/(sR-sL);
	}

	*F2L = F2 + 0.5*g*(hG*hG - hL*hL);
    *F2R = F2 + 0.5*g*(hD*hD - hR*hR);
}

//==============================Events to solve==============================

event init(t = 0) {
	// Initial conditions

	L0 = 4.;
	X0 = 0.;
	Y0 = 0.;
    double A = (a*a - r0*r0)/(a*a + r0*r0);

	// Riemann problem
//	double dist, dist2;
//	foreach() {
//		u.x[] = 0.;
//		u.y[] = 0.;
//
//        dist = sqrt(x*x + y*y); // for the wave
//        dist2 = sqrt((y-2)*(y-2)); // for the wall
//
//        zb[] = dist2 < 0.1 ? 14. : 0.;
//        h[] = dist < 1. ? 4. : 1.;
//	}

    // Radially-symmetrical paraboloid
    double r2;
    foreach() {
        r2 = (x-L0/2.)*(x-L0/2.) + (y-L0/2.)*(y-L0/2.);
        zb[] = -h0 * (1.-(r2/(a*a)));

        h[] = max(0., h0*(sqrt(1.-A*A)/(1.-A) - 1. - r2/(a*a) * ((1.-A*A)/((1.-A)*(1.-A)) - 1.)) - zb[]);
        zs[] = zb[] + h[];

        u.x[] = 0.;
        u.y[] = 0.;
    }

	DX = (L0/N); // cartesian grid
}

event plot (t += 0.1) {
	fprintf (stdout, "# t = %g\n", t); // time
	// paraview
	output_paraview (slist = {h, zb, zs}, vlist = {u});
}

event solve (i++) {	
	// Compute the solution for this iteration

	// Neumann boundary condition
	boundary ({zs, zb, h, u});

	// Fluxes at interfaces
	face vector Fh[];	// Fh.x,	Fh.y
	tensor FhuL[], FhuR[];		// Fhu.x.x	Fhu.x.y
				// Fhu.y.x	Fhu.y.y

	// Fluxes
	double vmax = 0.;
	foreach_face() {
		// Face in x AND y directions
		double hL = h[-1], zL = zb[-1], uL = u.x[-1], vL = u.y[-1]; // left value
		double hR = h[], zR = zb[], uR = u.x[], vR = u.y[];   // right value
		if (hL > dry || hR > dry) {
			// F(W_L,W_R) = (F1, F2)
			double F1, F2L, F2R, vmax_loc;
			// numerical flux
			HLL (hL, uL, zL, hR, uR, zR, &F1, &F2L, &F2R, &vmax_loc);
			// get fluxes
			Fh.x[] = F1;
			FhuL.x.x[] = F2L;
			FhuL.y.x[] = (F1 > 0. ? F1*vL : F1*vR);  // HLLC
            FhuR.x.x[] = F2R;
            FhuR.y.x[] = (F1 > 0. ? F1*vL : F1*vR);  // HLLC

			vmax = max(vmax, vmax_loc);
		} else {
			Fh.x[] = FhuL.x.x[] = FhuL.y.x[] = FhuR.x.x[] = FhuR.y.x[] = 0.;
		}
	}

	boundary_flux ({Fh, FhuL, FhuR});

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
        zs[] = h[] + zb[];

		// update u,v
		if (h[] > dry) {
			foreach_dimension() { // u, v
				double dFhuL = FhuR.x.x[]  + FhuR.x.y[]; // left + bottom
				double dFhuR = FhuL.x.x[1,0] + FhuL.x.y[0,1]; // right + top
				double hu = h_old*u.x[] - (dt/DX)*(dFhuR - dFhuL);
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

	// Default domain
	g = 1.; 
	L0 = 6.;
	X0 = -L0/2.; Y0 = -L0/2.;
	N = 1 << level; // 2^6
	CFL = 0.4;

    // Radially-symmetrical paraboloid
    h0 = 0.1, a = 1, r0 = 0.8;
    omega = sqrt(8.*g*h0)/a;
    TMAX = 3.*2.*M_PI/omega;

	run();
}
