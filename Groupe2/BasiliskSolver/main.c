//
// Created by matthieu on 25/05/2021.
//

#include "spherical.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#define MGD mgp.i
#include "layered/perfs.h"
scalar h;
vector u;
#include "terrain.h"
#include "okada.h"

#include "paraview2d.h"

#define MAXLEVEL 13
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define ETAMAXE  5e-2 // error on maximum free surface elevation (5 cm)

int main()
{
    // Earth radius in metres
    Radius = 6371220.;
    // the domain is 73 degrees squared
    size (73.);
    // centered on 142,38 longitude,latitude
    origin (142. - L0/2., 38. - L0/2.);
    // acceleration of gravity in m/min^2
    G = 9.81*sq(60.);
    // 32^2 grid points to start with
    init_grid (1 << MINLEVEL);

    run();
}

scalar etamax[];

int my_adapt() {
    #if QUADTREE
        scalar eta[];
        foreach()
            eta[] = h[] > dry ? h[] + zb[] : 0;
        boundary ({eta});
        astats s = adapt_wavelet ({eta, etamax}, (double[]){ETAE,ETAMAXE},
			    MAXLEVEL, MINLEVEL);
        fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
        return s.nf;
    #else // Cartesian
        return 0;
    #endif
}

event init (i = 0)
{
    terrain(zb, "~/Documents/terrain/etopo2", "~/Documents/terrain/srtm_japan", NULL);

    if (restore (file = "restart"))
        conserve_elevation();
    else {
        conserve_elevation();

        foreach()
        h[] = max(0., - zb[]);
        #include "tohoku/faults.h"

        boundary ({h});
    }

//    u.n[left]   = - radiation(0);
//    u.n[right]  = + radiation(0);
//    u.t[bottom] = - radiation(0);
//    u.t[top]    = + radiation(0);

    u.n[left] = neumann(0.); u.n[right] = neumann(0.);
    u.n[bottom] = neumann(0.); u.n[top] = neumann(0.);

    boundary({u, zb});
}

event friction (i++)
{
    foreach() {
        double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
        foreach_dimension()
            u.x[] /= a;
        if (h[] > dry && h[] + zb[] > etamax[])
            etamax[] = h[] + zb[];
    }
    boundary ({etamax, u});
}

//event logfile (i++)
//{
//    stats s = statsf (h);
//    norm n = normf (u.x);
//    if (i == 0)
//    fprintf (stderr,
//        "t i h.min h.max h.sum u.x.rms u.x.max dt mgD.i speed tn\n");
//        fprintf (stderr, "%g %d %g %g %g %g %g %g %d %g %ld\n",
//        t, i, s.min, s.max, s.sum, n.rms, n.max, dt, MGD,
//        perf.speed, grid->tn);
//}
//
//event snapshots (t += 15) {
//    char name[80];
//    sprintf (name, "dump-%g", t);
//    dump (name);
//}
//
//event flooding (t = 390)
//{
//    FILE * fp = fopen ("flooding", "w");
//    output_field ({h,zb,etamax}, fp,
//                  box = {{140.4,37.01},{142.44,39.11}},
//                  n = 512, linear = true);
//}
//
//event figures (t = 60; t <= 180; t += 60)
//{
//    scalar m[], etam[];
//    foreach() {
//        etam[] = eta[]*(h[] > dry);
//        m[] = etam[] - zb[];
//    }
//    boundary ({m, etam});
//    char name[80];
//    sprintf (name, "eta-%g.png", t);
//    output_ppm (etam, mask = m, min = -1, max = 2, file = name, n = 1024,
//            linear = true, box = {{123,14},{177,55}},
//            opt = "-fill white -opaque black");
//
//    sprintf (name, "level-%g.png", t);
//    scalar l[];
//    foreach()
//        l[] = level;
//    output_ppm (l, min = 5, max = 13, file = name, n = 1024,
//    linear = false, box = {{123,14},{177,55}});
//
//    if (t == 60)
//        output_ppm (etam, mask = m, min = -1, max = 2, file = "zoom-60.png",
//        n = 1024,
//        linear = true, box = {{140.6,30.8},{154.6,42}},
//        opt = "-fill white -opaque black");
//    else if (t == 120)
//        output_ppm (etam, mask = m, min = -1, max = 2, file = "zoom-120.png",
//        n = 1024,
//        linear = true, box = {{151.8,27.56},{165.8,38.68}},
//        opt = "-fill white -opaque black");
//    else if (t == 180)
//        output_ppm (etam, mask = m, min = -1, max = 2, file = "zoom-180.png",
//        n = 1024,
//        linear = true, box = {{158.8,24.25},{172.7,35.3}},
//        opt = "-fill white -opaque black");
//}
//
//event movies (t += 0.5)
//{
//    scalar m[], etam[];
//    foreach() {
//        m[] = eta[]*(h[] > dry) - zb[];
//        etam[] = h[] < dry ? 0. : (zb[] > 0. ? eta[] - zb[] : eta[]);
//    }
//    boundary ({m, etam});
//    output_ppm (etam, mask = m, min = -1, max = 2,
//            n = 1024, linear = true, file = "eta.mp4");
//
//    output_ppm (etam, mask = m, min = -5, max = 5,
//        n = 1024, linear = true,
//        box = {{140.4,37.51},{142.44,38.61}},
//        file = "eta-sendai.mp4");
//
//    output_ppm (etam, mask = m, min = -5, max = 5,
//        n = 1024, linear = true,
//        box = {{140.5,36.53},{142.52,37.62}},
//        file = "eta-fukushima.mp4");
//
//    output_ppm (etam, mask = m, min = -5, max = 5,
//        n = 1024, linear = true,
//        box = {{141.32,39.42},{143.44,40.50}},
//        file = "eta-miyako.mp4");
//
//    output_ppm (etam, mask = m, min = -5, max = 5,
//        n = 1024, linear = true,
//        box = {{141.11,38.51},{143.19,39.60}},
//        file = "eta-ofunato.mp4");
//
//    scalar l = etam;
//    foreach()
//        l[] = level;
//    output_ppm (l, min = MINLEVEL, max = MAXLEVEL, n = 1024,
//        file = "level.mp4");
//
//    #if _OPENMP
//        foreach()
//            etam[] = pid();
//        double tmax = omp_get_max_threads() - 1;
//        output_ppm (etam, max = tmax, n = 512, file = "pid.mp4");
//    #endif // _OPENMP
//}

#include "tohoku/gauges.h"

event adapt (i++) my_adapt();

event plot (t <= 10 ; t += 0.5)
{
    fprintf (stdout, "# t = %g\n", t); // time
    // paraview
    scalar zs[];
    foreach() zs[] = h[] + zb[];
        output_paraview (slist = {h, zs, zb}, vlist = {u});
}
