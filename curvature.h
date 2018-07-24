/**
# Curvature of an interface

The curvature field is defined only in interfacial cells. In all the
other cells it takes the value *nodata*. 

On trees, we need to redefine the restriction function to take
this into account i.e. the curvature of the parent cell is the average
of the curvatures in the interfacial child cells. */

#if TREE
static void curvature_restriction(Point point, scalar kappa)
{
  double k = 0., s = 0.;
  foreach_child()
    if (kappa[] != nodata)
      k += kappa[], s++;
  kappa[] = s ? k/s : nodata;
}

/**
The prolongation function performs a similar averaging, but using the
same stencil as that used for bilinear interpolation, so that the
symmetries of the volume fraction field and curvature field are
preserved. */

static void curvature_prolongation (Point point, scalar kappa)
{
  foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)
    #if dimension > 1
      for (int j = 0; j <= 1; j++)
    #endif
      #if dimension > 2
	for (int k = 0; k <= 1; k++)
      #endif
	  if (coarse(kappa,child.x*i,child.y*j,child.z*k) != nodata)
	    sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    kappa[] = s ? sk/s : nodata;
  }
}
#endif // TREE

/**
## Height-function curvature

To compute the curvature, we estimate the derivatives of the height
functions in a given direction (*x*, *y* or *z*). We first check that
all the heights are defined and that their orientations are the
same. We then compute the curvature as
$$
\kappa = \frac{h_{xx}}{(1 + h_x^2)^{3/2}}
$$
in two dimensions, or
$$
\kappa = \frac{h_{xx}(1 + h_y^2) + h_{yy}(1 + h_x^2) - 2h_{xy}h_xh_y}
{(1 + h_x^2 + h_y^2)^{3/2}}
$$
in three dimensions. */

#include "heights.h"

#if dimension == 2
foreach_dimension()
static double kappa_y (Point point, vector h) {
  int ori = orientation(h.y[]);
  for (int i = -1; i <= 1; i++)
    if (h.y[i] == nodata || orientation(h.y[i]) != ori)
      return nodata;
  double hx = (h.y[1] - h.y[-1])/2.;
  double hxx = (h.y[1] + h.y[-1] - 2.*h.y[])/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}
#else // dimension == 3
foreach_dimension()
static double kappa_z(Point point, vector h) {
  int ori = orientation(h.z[]);
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (h.z[i,j] == nodata || orientation(h.z[i,j]) != ori)
	return nodata;
  double hx = (h.z[1] - h.z[-1])/2.;
  double hy = (h.z[0,1] - h.z[0,-1])/2.;
  double hxx = (h.z[1] + h.z[-1] - 2.*h.z[])/Delta;
  double hyy = (h.z[0,1] + h.z[0,-1] - 2.*h.z[])/Delta;
  double hxy = (h.z[1,1] + h.z[-1,-1] - h.z[1,-1] - h.z[-1,1])/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#endif

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the curvature. This is done by the function below which
returns the HF curvature given a volume fraction field *c* and a
height function field *h*. */

static double height_curvature (Point point, scalar c, vector h)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *c*) and corresponding HF curvature function *kappa*
  (defined above). */

  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;  
  foreach_dimension()
    n.x.n = c[1] - c[-1], n.x.kappa = kappa_x;
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);
  
  /**
  We sort these pairs in decreasing order of $|n|$. */
  
  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#if dimension == 3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormKappa, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormKappa, n.y, n.z);
#endif

  /**
  We try each curvature function in turn. */

  double kappa = nodata;
  foreach_dimension()
    if (kappa == nodata) {
      kappa = n.x.kappa (point, h);
      if (kappa != nodata) {
	kappaf = n.x.kappa;
	if (n.x.n < 0.)
	  kappa = - kappa;
      }
    }

  if (kappa != nodata) {
    
    /**
     We limit the maximum curvature to $1/\Delta$. */
	
    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
    
    /**
     We add the axisymmetric curvature if necessary. */
      
#if AXI
    double nr, r = y, hx;
    if (kappaf == kappa_x) {
      hx = (height(h.x[0,1]) - height(h.x[0,-1]))/2.;
      nr = hx*(orientation(h.x[]) ? 1 : -1);
    }
    else {
      r += height(h.y[])*Delta;
      hx = (height(h.y[1,0]) - height(h.y[-1,0]))/2.;
      nr = orientation(h.y[]) ? -1 : 1;
    }
    /* limit the minimum radius to half the grid size */
    kappa += nr/max (sqrt(1. + sq(hx))*r, Delta/2.);
#endif
  }
  
  return kappa;
}

/**
## Parabolic fit of "mixed" height-functions

When the standard height function curvature calculation is not
possible (for example because not enough heights are available in any
given direction), one can try to combine all the available heights
(thus using "mixed" directions) to obtain points on the
interface. These point locations can then be fitted with a parabola
(using least-mean-square optimisation) and the resulting curvature can
be computed. The fitting functions are defined in the file included
below. */

#include "parabola.h"

/**
Given *n* (interface) point coordinates, this function returns the
number of "independent" points i.e. points which are more than
half-a-cell distant from all the other points. */

static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      foreach_dimension()
	d2 += sq(p[i].x - p[j].x);
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}

/**
Given a volume fraction field *c* and a height function field *h*,
this function returns the "mixed heights" parabola-fitted curvature
(or *nodata* if the curvature cannot be computed). */

static double height_curvature_fit (Point point, scalar c, vector h)
{

  /**
  The coordinates of the interface points and the number of
  interface points. */
  
  coord ip[dimension == 2 ? 6 : 27];
  int n = 0;

  /**
  We collect the points along all directions. */
  
  foreach_dimension() {

    /**
    We don't want to mix heights with different orientations. We first
    find the "dominant" orientation *ori*. */
    
    int n1 = 0, n2 = 0;
#if dimension == 2
    for (int i = -1; i <= 1; i++)
      if (h.y[i] != nodata) {
	if (orientation(h.y[i])) n1++; else n2++;
      }
#else // dimension == 3
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h.z[i,j] != nodata) {
	  if (orientation(h.z[i,j])) n1++; else n2++;
	}
#endif
    int ori = (n1 > n2);

    /**
    We look for height-functions with the dominant orientation and
    store the corresponding interface coordinates (relative to the
    center of the cell and normalised by the cell size). */

#if dimension == 2
    for (int i = -1; i <= 1; i++)
      if (h.y[i] != nodata && orientation(h.y[i]) == ori)
	ip[n].x = i, ip[n++].y = height(h.y[i]);
#else // dimension == 3
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h.z[i,j] != nodata && orientation(h.z[i,j]) == ori)
	  ip[n].x = i, ip[n].y = j, ip[n++].z = height(h.z[i,j]);
#endif
  }

  /**
  If we don't have enough independent points, we cannot do the
  parabolic fit. */
  
  if (independents (ip, n) < (dimension == 2 ? 3 : 9))
    return nodata;

  /**
  We recover the interface normal and the centroid of the interface
  fragment and initialize the parabolic fit. */
  
  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (c[], m);
  double area = plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);
#if dimension == 2
  NOT_UNUSED(area);
  parabola_fit_add (&fit, fc, PARABOLA_FIT_CENTER_WEIGHT);
#else // dimension == 3
  parabola_fit_add (&fit, fc, area*100.);
#endif
  
  /**
  We add the collected interface positions and compute the
  curvature. */

  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}

/**
## Parabolic fit of centroids

If all else fails, we try a parabolic fit of interface centroids. */

static double centroids_curvature_fit (Point point, scalar c)
{

  /**
  We recover the interface normal and the centroid of the interface
  fragment and initialize the parabolic fit. */
  
  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (c[], m);
  plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  /**
  We add the interface centroids in a $3^d$ neighborhood and compute
  the curvature. */

  coord r = {x,y,z};
  foreach_neighbor(1)
    if (c[] > 0. && c[] < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = plane_alpha (c[], m);
      double area = plane_area_center (m, alpha, &fc);
      coord rn = {x,y,z};
      foreach_dimension()
	fc.x += (rn.x - r.x)/Delta;
      parabola_fit_add (&fit, fc, area);
    }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}

/**
## General curvature computation

We first need to define "interfacial cells" i.e. cells which contain
an interface. A simple test would just be that the volume fraction is
neither zero nor one. As usual things are more complicated because of
round-off errors. They can cause the interface to be exactly aligned
with cell boundaries, so that cells on either side of this interface
have fractions exactly equal to zero or one. The function below takes
this into account. */

static inline bool interfacial (Point point, scalar c)
{
  if (c[] >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	if (c[i] <= 0.)
	  return true;
  }
  else if (c[] <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	if (c[i] >= 1.)
	  return true;
  }
  else // c[] > 0. && c[] < 1.
    return true;
  return false;
}

/**
The function below computes the mean curvature *kappa* of the
interface defined by the volume fraction *c*. It uses a combination of
the methods above: statistics on the number of curvatures computed
which each method is returned in a *cstats* data structure. 

If *sigma* is different from zero the curvature is multiplied by *sigma*.

If *add* is *true*, the curvature (optionally multiplied by *sigma*)
is added to field *kappa*. */

typedef struct {
  int h; // number of standard HF curvatures
  int f; // number of parabolic fit HF curvatures
  int a; // number of averaged curvatures
  int c; // number of centroids fit curvatures
} cstats;

struct Curvature {
  scalar c, kappa;
  double sigma;
  bool add;
};

trace
cstats curvature(struct Curvature p)
{
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  cstats s = {0,0,0,0};
  vector h[];
  heights (c, h);

  /**
  On trees we set the prolongation and restriction functions for
  the curvature. */
  
#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

  /**
  We first compute a temporary curvature *k*: a "clone" of
  $\kappa$. */
  
  scalar k[];
  scalar_clone (k, kappa);

  foreach() {

    /**
    If we are not in an interfacial cell, we set $\kappa$ to *nodata*. */

    if (!interfacial (point, c))
      k[] = nodata;

    /**
    Otherwise we try the standard HF curvature calculation first, and
    the "mixed heights" HF curvature second. */ 
    
    else if ((k[] = height_curvature (point, c, h)) != nodata)
      s.h++;
    else if ((k[] = height_curvature_fit (point, c, h)) != nodata)
      s.f++;
  }
  boundary ({k});

  foreach() {
    
    /**
    We then construct the final curvature field using either the
    computed temporary curvature... */

    double kf;
    if (k[] < nodata)
      kf = k[];
    else if (interfacial (point, c)) {

      /**
      ...or the average of the curvatures in the $3^{d}$ neighborhood
      of interfacial cells. */
      
      double sk = 0., a = 0.;
      foreach_neighbor(1)
	if (k[] < nodata)
	  sk += k[], a++;
      if (a > 0.)
	kf = sigma*sk/a, s.a++;
      else

	/**
	Empty neighborhood: we try centroids as a last resort. */
	
	kf = sigma*centroids_curvature_fit (point, c), s.c++;
    }
    else
      kf = nodata;

    /**
    We add or set *kappa*. */
    
    if (kf == nodata)
      kappa[] = nodata;
    else if (p.add)
      kappa[] += sigma*kf;
    else
      kappa[] = sigma*kf;      
  }
  boundary ({kappa});

  return s;
}

/**
# Position of an interface

This is similar to curvature but this time for the position of the
interface, defined as
$$
pos = \mathbf{G}\cdot(\mathbf{x} - \mathbf{Z})
$$
with $\mathbf{G}$ and $\mathbf{Z}$ two vectors and $\mathbf{x}$ the
coordinates of the interface.

This is defined only in interfacial cells. In all the other cells it
takes the value *nodata*.

We first need a function to compute the position $\mathbf{x}$ of an
interface. For accuracy, we first try to use height functions. */

foreach_dimension()
static double pos_x (Point point, vector h, coord * G, coord * Z) {
  if (h.x[] == nodata)
    return nodata;
  coord o = {x, y, z};
  o.x += height(h.x[])*Delta;
  double pos = 0.;
  foreach_dimension()
    pos += (o.x - Z->x)*G->x;
  return pos;
}

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the position. This is done by the function below which returns
the HF position given a volume fraction field *f*, a height function
field *h* and vectors *G* and *Z*. */

static double height_position (Point point, scalar f, vector h,
			       coord * G, coord * Z)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *f*) and corresponding HF position function *pos*
  (defined above). */

  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  foreach_dimension()
    n.x.n = f[1] - f[-1], n.x.pos = pos_x;
  
  /**
  We sort these pairs in decreasing order of $|n|$. */
  
  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormPos, n.x, n.y);
#if dimension == 3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormPos, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormPos, n.y, n.z);
#endif

  /**
  We try each position function in turn. */

  double pos = nodata;
  foreach_dimension()
    if (pos == nodata)
      pos = n.x.pos (point, h, G, Z);

  return pos;
}

/**
The position() function fills field *pos* with
$$
\mathbf{G}\cdot(\mathbf{x} - \mathbf{Z})
$$
with $\mathbf{x}$ the position of the interface defined by $f$.

If *add* is *true*, the position is added to *pos*. */

struct Position {
  scalar f, pos;
  coord G, Z;
  bool add;
};

void position(struct Position p)
{
  scalar f = p.f, pos = p.pos;
  coord * G = &p.G, * Z = &p.Z;

  /**
  On trees we set the prolongation and restriction functions for
  the position. */
  
#if TREE
  pos.refine = pos.prolongation = curvature_prolongation;
  pos.restriction = curvature_restriction;
#endif
  
  vector h[];
  heights (f, h);
  foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z);
      if (hp == nodata) {

	/**
	If the height function is not defined, we use the centroid of
	the reconstructed VOF interface. */
	
	coord n = mycs (point, f), o = {x,y,z}, c;
	double alpha = plane_alpha (f[], n);
	plane_area_center (n, alpha, &c);
	hp = 0.;
	foreach_dimension()
	  hp += (o.x + Delta*c.x - Z->x)*G->x;
      }
      if (p.add)
	pos[] += hp;
      else
	pos[] = hp;
    }
    else
      pos[] = nodata;
  }
  boundary ({pos});  
}

struct Positionbumptrans {
  scalar f, pos; 
  double a, Po, t, s, b;
  coord G, Z;
  bool add;
};

void positionbumptrans(struct Positionbumptrans p)
{
  scalar f = p.f, pos = p.pos;
  double a = p.a, Po = p.Po, t = p.t, s = p.s , b = p.b; 
  coord * G = &p.G, * Z = &p.Z;


  /**
  On trees we set the prolongation and restriction functions for
  the position. */
  
#if TREE
  pos.refine = pos.prolongation = curvature_prolongation;
  pos.restriction = curvature_restriction;
#endif
  vector h[];
  heights (f, h);
  foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z)+Po*exp(-(x-(b+t*s))*(x-(b+t*s))/(2*a*a));
      if (hp == nodata) {
	/**
	If the height function is not defined, we use the centroid of
	the reconstructed VOF interface. */
	
        coord n = mycs (point, f), o = {x,y,z}, c;
	double alpha = plane_alpha (f[], n);
	plane_area_center (n, alpha, &c);
	hp = 0.;
	foreach_dimension()
	  hp += ((o.x + Delta*c.x - Z->x)*G->x)+Po*exp(-(x-(b+t*s))*(x-(b+t*s))/(2*a*a));
      }
      if (p.add)
	pos[] += hp;
      else
	pos[] = hp;
    }
    else
      pos[] = nodata;
  }
  boundary ({pos});  
}

struct Positionbumpstat {
  scalar f, pos; 
  double a, Po, b;
  coord G, Z;
  bool add;
};

void positionbumpstat(struct Positionbumpstat p)
{
  scalar f = p.f, pos = p.pos;
  double a = p.a, Po = p.Po, b = p.b; 
  coord * G = &p.G, * Z = &p.Z;


  /**
  On trees we set the prolongation and restriction functions for
  the position. */
  
#if TREE
  pos.refine = pos.prolongation = curvature_prolongation;
  pos.restriction = curvature_restriction;
#endif
  vector h[];
  heights (f, h);
  foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z)+Po*exp(-(x-b)*(x-b)/(2*a*a));
      if (hp == nodata) {
	/**
	If the height function is not defined, we use the centroid of
	the reconstructed VOF interface. */
	
        coord n = mycs (point, f), o = {x,y,z}, c;
	double alpha = plane_alpha (f[], n);
	plane_area_center (n, alpha, &c);
	hp = 0.;
	foreach_dimension()
	  hp += ((o.x + Delta*c.x - Z->x)*G->x)+Po*exp(-(x-b)*(x-b)/(2*a*a));
      }
      if (p.add)
	pos[] += hp;
      else
	pos[] = hp;
    }
    else
      pos[] = nodata;
  }
  boundary ({pos});  
}


struct Positionbumpstat2D {
  scalar f, pos; 
  double a, Po, b;
  coord G, Z;
  bool add;
};

void positionbumpstat2D(struct Positionbumpstat2D p)
{
  scalar f = p.f, pos = p.pos;
  double a = p.a, Po = p.Po, b = p.b; 
  coord * G = &p.G, * Z = &p.Z;

  /**
  On trees we set the prolongation and restriction functions for
  the position. */
  
#if TREE
  pos.refine = pos.prolongation = curvature_prolongation;
  pos.restriction = curvature_restriction;
#endif
  vector h[];
  heights (f, h);
  foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z)+Po*exp(-((x-b)*(x-b)+(y-b)*(y-b))/(2*a*a));
      if (hp == nodata) {
	/**
	If the height function is not defined, we use the centroid of
	the reconstructed VOF interface. */
	
        coord n = mycs (point, f), o = {x,y,z}, c;
	double alpha = plane_alpha (f[], n);
	plane_area_center (n, alpha, &c);
	hp = 0.;
	foreach_dimension()
	  hp += ((o.x + Delta*c.x - Z->x)*G->x)+Po*exp(-((x-b)*(x-b)+(y-b)*(y-b))/(2*a*a));
      }
      if (p.add)
	pos[] += hp;
      else
	pos[] = hp;
    }
    else
      pos[] = nodata;
  }
  boundary ({pos});  
}

struct Positionbumptrans2D {
  scalar f, pos; 
  double a, Po, b, t, s;
  coord G, Z;
  bool add;
};

void positionbumptrans2D(struct Positionbumptrans2D p)
{
  scalar f = p.f, pos = p.pos;
  double a = p.a, Po = p.Po, b = p.b, t =p.t, s = p.s; 
  coord * G = &p.G, * Z = &p.Z;

  /**
  On trees we set the prolongation and restriction functions for
  the position. */
  
#if TREE
  pos.refine = pos.prolongation = curvature_prolongation;
  pos.restriction = curvature_restriction;
#endif
  vector h[];
  heights (f, h);
  foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z)+Po*exp(-((x-cos(s*pi*t)*b)*(x-cos(s*pi*t)*b)+(y-sin(s*pi*t)*b)*(y-sin(s*pi*t)*b))/(2*a*a));
      if (hp == nodata) {
	/**
	If the height function is not defined, we use the centroid of
	the reconstructed VOF interface. */
	
        coord n = mycs (point, f), o = {x,y,z}, c;
	double alpha = plane_alpha (f[], n);
	plane_area_center (n, alpha, &c);
	hp = 0.;
	foreach_dimension()
	  hp += ((o.x + Delta*c.x - Z->x)*G->x)+Po*exp(-((x-cos(s*pi*t)*b)*(x-cos(s*pi*t)*b)+(y-sin(s*pi*t)*b)*(y-sin(s*pi*t)*b))/(2*a*a));
      }
      if (p.add)
	pos[] += hp;
      else
	pos[] = hp;
    }
    else
      pos[] = nodata;
  }
  boundary ({pos});  
}
