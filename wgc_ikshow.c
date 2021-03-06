/**
# Gravito-Capillary Waves

A similar test to the [capillary wave](capwave.c) but for a pure
gravity wave, using the [reduced gravity](/src/reduced.h) approach.

We use a constant-resolution grid, the Navier--Stokes solver for
two-phase flows and reduced gravity. */

//#include "grid/cartesian1D.h"
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
//#include "test/prosperetti-gravity.h"
#include "iforce.h"
#include "curvature.h"
//#include "curvaturebump.h"


/**
We make sure that the boundary conditions for the face-centered
velocity field are consistent with the centered velocity field (this
affects the advection term). */

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

/**
We will store the accumulated error in *se* and the number of samples
in *ne*. */

int main() {
  /**
  The domain is 2x2 to minimise finite-size effects. The viscosity is
  constant. The acceleration of gravity is 50. */
  
  L0 = 2.;
  Y0 = -L0/2.;
  G.y = 50.; //'G' has an attribute 'y' to specify direction
  f.sigma = G.y*(L0/10.)*(L0/10.); //'f' interface has an attribute 'sigma' for surface tension
  rho1 = 1, rho2 = 0.1225;
  mu1 = mu2 = 0.000; //Inviscid Analysis
  TOLERANCE = 1e-6;
  /**
  We vary the resolution to check for convergence. */
  //  for (N = 16; N <= 128; N *= 2) {
  N = 128;
  //for(f.sigma = 0; f.sigma <=  5; f.sigma = f.sigma + 0.25)
  //{  
  run();
  //}
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity. */
event init (t = 0) {
  //fraction (f, y - 0.1*sin(1.*pi*x));
  double a = 0.4;
  fraction(f,y - 0.1*exp(-(x-1)*(x-1)/(2*a*a)))
}

/**
By default tracers are defined at $t-\Delta t/2$. We use the *first*
keyword to move VOF advection before the *amplitude* output i.e. at
$t+\Delta/2$. This improves the results. */

event vof (i++, first);

/**
We output the amplitude at times matching exactly those in the
reference file. */

event amplitude (t += 0.005; t <= 1.) {
  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  output_ppm (pos,linear = false);
  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points *N*. */
  char name[100];
  double st;
  st = f.sigma;
  sprintf (name, "Kakade/wgc_2d/NewStuff/wave_nu0_new-%d-%g.txt", N,st);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t*16.032448313657, max);
  fflush (fp);

  if (N == 64)
    output_facets (f, stdout);
    sprintf (name, "Kakade/wgc_2d/wave-%d.gfs", i);
    output_gfs (file = name);
}

/**
At the end of the simulation, we output on standard error the
resolution (number of grid points per wavelength) and the relative RMS
error. */

event error (t = end)
  fprintf (stderr, "%g %g\n", N/L0,f.sigma); //Number of grid points per wavelength

#if 0
event gfsview (i += 1) {
  static FILE * fp = popen ("gfsview2D -s gravity.gfv", "w");
  output_gfs (fp);
}
#endif

event images (t += 0.005) {
  scalar pos[];
  position (f, pos, {0,1});
  output_ppm (pos);
}

event acceleration (t++)
{
  scalar phi = f.phi;
  coord G1;
  G1.y = -G.y;
  positionbump (f, phi, G1, Z, add = false);
  f.phi = phi;
}



