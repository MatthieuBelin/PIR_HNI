// qcc -O2 -Wall solver.c -o solver -lm

#include "grid/quadtree.h"
#include "paraview2d.h"
#include "run.h"


//==============================Variables==============================

scalar h[]; // Height of the water
scalar zb[]; // Bathymetrie
scalar zs[]; // Free surface
vector u[]; // Speed in x and y direction
face vector fh[];
tensor fu[];

double g = 9.81;

double TMAX = 1.;

double DX;
double dry = 1.e-10;
int level = 6;

//==============================Main function==============================

int main(int argc, char * argv[]) {
    // resolution
    if (argc > 1) level = atoi(argv[1]);

//    // Default domain
//    g = 1.;
//    L0 = 6.;
//    X0 = -L0/2.; Y0 = -L0/2.;

    L0 = 4.;
    X0 = 0.;
    Y0 = 0.;

    N = 1 << level; // 2^6
    init_grid(N);
    CFL = 0.4;
    run();
}

//==============================boundary conditions==============================
u.n[left] = neumann(0.); u.n[right] = neumann(0.);
u.n[bottom] = neumann(0.); u.n[top] = neumann(0.);

//==============================Functions for the solver==============================

void flux (double h1, double u1, double *flux1, double *flux2) {
	*flux1 = h1*u1;
	*flux2 = h1*u1*u1 + g*h1*h1/2.;
}

// Roe's state
void WRoe (double hL, double uL, double hR, double uR, // inputs
           double * hRoe, double * uRoe ) {
    // Roe's state
    *hRoe = (hL + hR)/2.;
    *uRoe = (uL*sqrt(hL) + uR*sqrt(hR))/(sqrt(hL) + sqrt(hR));
}

void HLL (double hL, double uL, double hR, double uR, double *F1, double *F2, double *vmax) {
    // Roe's state
    double hRoe, uRoe;
    WRoe (hL, uL, hR, uR, &hRoe, &uRoe);

//    double zi = max(zL, zR);
//    double hL = max(0, hG + zL - zi);
//    double hR = max(0, hD + zR - zi);

	// wave speeds
	double sL = min(uL - sqrt(g*hL), uRoe - sqrt(g*hRoe));
	double sR = max(uR + sqrt(g*hR), uRoe + sqrt(g*hRoe));
	*vmax = max(fabs(sL), fabs(sR));

	// fluxes
	double f1L, f2L;
	flux(hL, uL, &f1L, &f2L);
	double f1R, f2R;
	flux(hR, uR, &f1R, &f2R);

//	double F2;
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
//
//	*F2L = F2 + 0.5*g*(hG*hG - hL*hL);
//    *F2R = F2 + 0.5*g*(hD*hD - hR*hR);
}

//==============================Events to solve==============================

/* FOR MESH ADAPTATION (see http://basilisk.fr/src/saint-venant.h) */
event defaults (i=0) {
    // default initial condition
    foreach () {
        for (scalar s in {h,u,zb}) {
            s[] = 0.; // initial condition
        }
    }

    // for refine/coarsen
    for (scalar s in {h,u,zb}) {
        // see http://basilisk.fr/src/grid/multigrid-common.h
        s.refine = s.prolongation = refine_linear; // refine
        s.restriction = restriction_volume_average; // coarsen
    }
}
/* END MESH ADAPTATION */

event init(t = 0) {
	// Initial conditions

	/* dam-break */
	double dist, dist2;
	foreach() {
		u.x[] = 0.;
		u.y[] = 0.;

        dist = sqrt(x*x + y*y); // for the wave
        dist2 = sqrt((y-2)*(y-2)); // for the wall

        zb[] = dist2 < 0.1 ? 14. : 0.;
        h[] = dist < 1. ? 4. : 1.;
	}

//    // Radially-symmetrical paraboloid
//
//    const double h0 = 0.1, a = 1., r0 = 0.8;
//    double A = (a*a - r0*r0)/(a*a + r0*r0);
//    double omega = sqrt(8.*g*h0)/a;
//    TMAX = 3.*2.*M_PI/omega;
//
//    foreach() {
//        double r2 = (x-L0/2.)*(x-L0/2.) + (y-L0/2.)*(y-L0/2.);
//        zb[] = -h0 * (1.-(r2/(a*a)));
//        double eta = h0*(sqrt(1.-A*A)/(1.-A) - 1. - r2/(a*a) * ((1.-A*A)/((1.-A)*(1.-A)) - 1.));
//        h[] = max(0., eta - zb[]);
//        zs[] = zb[] + h[];
//
//        u.x[] = u.y[] = 0.;
//    }
}

event plot (t <= TMAX; t += 0.1) {
	fprintf (stdout, "# t = %g\n", t); // time
	// paraview
	output_paraview (slist = {h, zs, zb}, vlist = {u});
}

event adapt (i++) {
    astats s = adapt_wavelet ({h}, (double[]){1.e-3}, maxlevel=9, minlevel=7);
    fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

event integration (i++) {
	// Compute the solution for this iteration

	// Neumann boundary condition
	boundary ({zs, zb, h, u});

	// Fluxes at interfaces
	face vector Fh[];	// Fh.x,	Fh.y
	//tensor FhuL[], FhuR[];		// Fhu.x.x	Fhu.x.y
				// Fhu.y.x	Fhu.y.y
	//face vector FhuL[], FhuR[], Fhv[]; // Car l'antidiagonale est identique pour les deux tenseurs
    tensor Fhu[];   // FhuL.x.x, Fhu.x.y
                    // Fhu.y.x,  FhuR.y.y
    face vector Fhv[]; // Fhv = FhuR.x.x, FhuR.y.y

	// Fluxes
	double dtmax = 1.E+10;
	foreach_face(reduction (min:dtmax)) {
		// Face in x AND y directions
		double hL = h[-1], zL = zb[-1], uL = u.x[-1], vL = u.y[-1]; // left value
		double hR = h[], zR = zb[], uR = u.x[], vR = u.y[];   // right value

		// Audusse2004 --> M = Minus ; P = Plus
		double zi = max(zL, zR); // eq 2.9
		double hM = max(0., hL + zL - zi), uM = uL, vM = vL; // eq 2013
		double hP = max(0., hR + zR -zi), uP = uR, vP = vR;
		// source term
		double sM = g/2. * (sq(hM) - sq(hL)); // eq 2.14
		double sP = g/2. * (sq(hR) - sq(hP));
		if (hM > dry || hP > dry) {
			// F(W_L,W_R) = (F1, F2)
			double F1, F2, F3, vmax;
			// numerical flux
			HLL (hM, uM, hP, uP, &F1, &F2, &vmax); //(hL, uL, zL, hR, uR, zR, &F1, &F2L, &F2R, &vmax_loc);
			// get fluxes
            F3 = (F1 > 0. ? F1*vM : F1*vP);  // HLLC
			Fh.x[] = fm.x[] * F1;
			Fhu.x.x[] = fm.x[] * (F2 - sM);
            Fhu.y.x[] = fm.x[] * F3;
			Fhv.x[] = fm.x[] * (F2 + sP);

			// Time step
			double dxloc = Delta*cm[]/fm.x[];
			dtmax = min(dtmax, CFL*dxloc/vmax);
		} else {
			Fh.x[] = Fhu.x.x[] = Fhu.y.x[] = Fhv.x[] = 0.;
		}
	}

	boundary_flux ({Fh, Fhu, Fhv});

	// next timestep
	double dt = dtnext(dtmax);
	
	// Update
	foreach() {
	    double dxloc = cm[]*Delta; // size of the cell
		// Fh.x,    Fh.y
		// Fhu.x.x, Fhu.x.y
		// Fhu.y.x, Fhu.y.y

		// update h
		double h_old = h[]; // backup h at t = t^n 
		// note that: FL = F[], FR = F[1]
		double FhL = Fh.x[]  + Fh.y[]; // left + bottom
		double FhR = Fh.x[1,0] + Fh.y[0,1]; // right + top 
		h[] = h_old - (dt/dxloc)*(FhR - FhL); // h at t = t^{n+1}
        zs[] = h[] + zb[];

		// update u,v
		if (h[] > dry) {
			foreach_dimension() { // u, v
				double FhuL = Fhu.x.y[]  + Fhv.x[]; // left + bottom
				double FhuR = Fhu.x.x[1,0] + Fhu.x.y[0,1]; // right + top
				double hu = h_old*u.x[] - (dt/dxloc)*(FhuR - FhuL);
				u.x[] = hu/h[];
			}
		} else { // dry ==> set u = v = 0
			u.x[] = u.y[] = 0.;
		}
	}
}
