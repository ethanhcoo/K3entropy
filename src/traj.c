
/* Traj:  draws trajectories of K3 dynamical systems */
/* Writes to stdout a Postscript file */

#include "k3.h"
#include "cx.h" 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Colors */
#define BLACK 1
#define GRAY 2
#define Color 150 /*Ethan added this*/

/* Defaults */
#define TOPONLY 1

static complex center = CENTER;
static double radius = RADIUS;
static int idim = IDIM, jdim= JDIM;
static int iterlim = ITERLIM;
static int onlytop = TOPONLY;

int main(argc,argv)
	int argc;
	char *argv[];

{
	int i,n;
	double a,b,x,y;
	point p;

	init_color_map();
	scan_window(center,radius);
	if(argc < 2) random_orbits(50);

	for(i=1; i<argc; i++)
	{	if(argv[i][0] != '-') usage();
		switch(argv[i][1])
	{
		case 'c':
		if(argc <= (i+2)) usage();
		sscanf(argv[++i],"%lf",&center.x);
		sscanf(argv[++i],"%lf",&center.y);
		scan_window(center,radius);
		break;

		case 'i':
		if(argc <= (i+1)) usage();
		sscanf(argv[++i],"%d",&iterlim);
		break;

		case 'd':
		if(argc <= (i+2)) usage();
		sscanf(argv[++i],"%d",&idim);
		sscanf(argv[++i],"%d",&jdim);
		scan_dim(idim,jdim);
		break;

		case 'o':
		if(argc <= (i+2)) usage();
		sscanf(argv[++i],"%lf",&x);
		sscanf(argv[++i],"%lf",&y);
		p = lift(x,y);
		orbit(p);
		break;

		case 'p':
		if(argc <= (i+2)) usage();
		sscanf(argv[++i],"%lf",&a);
		sscanf(argv[++i],"%lf",&b);
		setk3(a,b);
		break;

		case 'r':
		if(argc <= (i+1)) usage();
		sscanf(argv[++i],"%lf",&radius);
		scan_window(center,radius);
		break;

		case 'R':
		if(argc <= (i+1)) usage();
		sscanf(argv[++i],"%d",&n);
		random_orbits(n);
		break;

		case 's':
		init_surface();
		break;

		case 't':
		onlytop = !onlytop;
		break;

		default:
		usage();
	}
	}
	scan_end(argc,argv);
	exit(0);
}

void usage()
{
	fprintf(stderr,"Usage:  traj [options]\n");
	fprintf(stderr,"  -c [x y]   window center\n");
	fprintf(stderr,"  -d [n m]   nxm raster\n");
	fprintf(stderr,"  -i [lim]   iteration limit\n");
	fprintf(stderr,"  -o [x y]   draw real orbit\n");
	fprintf(stderr,"  -p [a b] K3 parameters\n");
	fprintf(stderr,"  -r [r]     window radius\n");
	fprintf(stderr,"  -R [n]     random orbits\n");
	fprintf(stderr,"  -s         color surface light gray\n");
	fprintf(stderr,"  -t         toggle top sheet filter\n");
	fprintf(stderr,"Postscript file written to stdout\n");
	exit(1);
}

/* Color map */
void init_color_map()
{
	scan_gray(0,1.0);
	scan_gray(BLACK,0.0);
	scan_gray(GRAY,0.99);
}

/* Initialize surface */
/* Turn points inside surface gray */
void init_surface()
{
	complex p;
	double dx, dy;

	dx = radius/idim;
	dy = radius/jdim;

	for(p.x=center.x-radius; p.x<center.x+radius; p.x += dx)
	for(p.y=center.y-radius; p.y<center.y+radius; p.y += dy)
		if(!unreal(p.x,p.y)) scan_set(p,Color); /*Ethan changed GRAY to COLOR*/
}

/* Draw orbit of a point */
/* Draw only upper sheet if point is real and onlytop is set */
void orbit(p)
	point p;
{
	int i;
	complex proj();

	for(i=0; i<iterlim; i++)
	{	if(!onlytop || top(p))
			scan_set(proj(p),Color); /*Ethan changed BLACK to COLOR*/
		p = f(p);
	}
}

/* Draw orbits of n random points in window */
/* Keeps looking until a real point is found */
void random_orbits(n)
	int n;
{
	point p;
	double random_real(), x, y;

	int i,j;
	for(i=0; i<n; i++)
	{	for(j=0; j<1000; j++)
		{	x = center.x + radius*(2*random_real() -1);
			y = center.y + radius*(2*random_real() -1);
			if(!unreal(x,y)) break;
		}
		if(unreal(x,y))
		{	fprintf(stderr,"Couldn't find any real points!\n");
			fprintf(stderr,"Check window\n");
			exit(2);
		}
	p = lift(x,y);
	orbit(p);
	}
}


/* Random real in [0,1] */
#define NORM 2147483647.0
double random_real()
{
	double x;
	long random();
	x = random()/NORM;
	return(x);
}

/* Project to complex x-y plane */
complex proj(p)
	point p;
{
	complex z;

	z.x = p.x;
	z.y = p.y;
	return(z);
}
