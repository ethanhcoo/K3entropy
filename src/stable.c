
/* Stable:  draws stable manifold of K3 dynamical systems */
/* Writes to stdout a Postscript file */

#include "k3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


static double centerx = CENTERX, centery = CENTERY, radius = RADIUS;
static double length=LENGTH, sep = SEP;
static int drawtop = 1, drawbot = 1, verbose=VERBOSE, spin=0, unravel = 0;

/* Very large array to dynamically store paths*/
point ps[PSMAX];

int main(int argc, char *argv[]) {
	int i,n;
	double a,b;

	//read in values from stable.run
	for(i=1; i<argc; i++)
	{	if(argv[i][0] != '-') usage();
		switch(argv[i][1])
	{
		case 'b':
		drawtop = 0;
		break;

		case 'c':
		if(argc <= (i+2)) usage();
		sscanf(argv[++i],"%lf",&centerx);
		sscanf(argv[++i],"%lf",&centery);
		break;

		case 'e':
		if(argc <= (i+1)) usage();
		sscanf(argv[++i],"%lf",&sep);
		setsep(sep);
		break;

		case 'l':
		if(argc <= (i+1)) usage();
		sscanf(argv[++i],"%lf",&length);
		break;

		case 'p':
		if(argc <= (i+2)) usage();
		sscanf(argv[++i],"%lf",&a);
		sscanf(argv[++i],"%lf",&b);
		setk3(a,b);
		break;

		case 'q':
		verbose = 0;
		break;

		case 'r':
		if(argc <= (i+1)) usage();
		sscanf(argv[++i],"%lf",&radius);
		break;

		case 's':
		spin = 1;
		break;

		case 't':
		drawbot = 0;
		break;

		case 'u': /*ethan added unravel*/
		unravel = 1;
		break;

		default:
		usage();
	}
	}

	//recurrent search corner 
	
	/*
	point p = lift(0.25900, 1.32300);
	for(int i = 0; i <= 1000; i++){
		if(i % 23 == 0){
			pprint(p);
		}
		p = f(p);
	}
	recurrent_search();
	abort(); */

	//open postscript file, declare window
	ps_open(argc,argv);
	ps_window(centerx,centery,radius);
	
	// draws background
	point starter = lift(.1,.1);
	//draw_background(starter,20000);

	//draw stable manifold: this will no longer work since I changed f from sym(x(y(z))) to x(y(z)), and thus f no longer has a fixed point at (0,0,1).
	//draw_manifold(); 

	//draw a small ball around periodic point, and plot its image
	point q = lift(0.89560, 0.00490); //0.89460, 0.00480
	double epsilon = .0001;
	double delta = .00001;
	point s;
	point l;
	s.x = q.x - epsilon;
	for(int i = 0; i < 2*epsilon/delta; i++){
		s.y = q.y - epsilon;
		for(int j = 0; j < 2*epsilon/delta; j++){
			s = lift(s.x, s.y);
			ps_dot_transparent(s);
			l.x = s.x;
			l.y = s.y;
			l.z = s.z;
			for(int m = 0; m < 11; m++){
				l = f(l);
			}
			//ps_dot(l, 0);
			s.y += delta;
		}
		s.x += delta;
	}
	s.x = q.x - epsilon;
	s.y = q.y - epsilon;
	for(int j = 0; j < 2*epsilon/delta; j++){
		s = lift(s.x, s.y);
		//ps_dot_transparent(s);
		l.x = s.x;
		l.y = s.y;
		l.z = s.z;
		for(int m = 0; m < 11; m++){
			l = f(l);
		}
		ps_dot(l, 0);
		s.y += delta;
	}
	
	
	
	//close postscript file
	ps_close();

	exit(0);
}

void usage()
{
	fprintf(stderr,"Usage:  stable [options]\n");
	fprintf(stderr,"  -b          draw only bottom\n");
	fprintf(stderr,"  -c [x y]    window center\n");
	fprintf(stderr,"  -e [eps]    accuracy\n");
	fprintf(stderr,"  -l [length] length of stable manifold drawn\n");
	fprintf(stderr,"  -p [a b]    K3 parameters\n");
	fprintf(stderr,"  -q          quiet\n");
	fprintf(stderr,"  -r [r]      window radius\n");
	fprintf(stderr,"  -s          spin picture\n");
	fprintf(stderr,"  -t          draw only top\n");
	fprintf(stderr,"Postscript file written to stdout\n");
	exit(1);
}

point initial_point()
{
	point p;
	p = lift(sep/10,-sep/10);
	return(p);
}
void draw_background(point z, int n)/*Ethan added this*/
{
	for (int i = 0; i<n; i++){
		ps_dot_transparent(z);
		z = f(z);
	}
}
void draw_manifold() //this will no longer work since I changed f from sym(x(y(z))) to x(y(z)), and thus f no longer has a fixed point at (0,0,1).
{
	int i, n;
	double tot, t;
	n=1;
	ps[0] = initial_point();
	ps[1] = f(ps[0]);
	draw_path(n);
	tot=pathlength(n);
	while(length > tot)
	{	n=fpath(ps,n);
		tot = tot + pathlength(n);
		di(ps,n);
		draw_path(n);
	}
	if(verbose) 
		fprintf(stderr,"Total path length:  %.2lf\n",tot);
}


void draw_path(int n)
{
	int i;
	point p, q;

	for(i=0; i<n; i++)
	{	p = ps[i]; q = ps[i+1];
/*		if(rotate) {p=rotate(p); q=rotate(q);} */
		if(drawtop &&  top(p)) ps_line(p,q, unravel);
		if(drawbot && !top(p)) ps_line_color(p,q, unravel); /*Ethan changed ps_line to ps_line_color to draw the bottom orange*/
	}
}

double pathlength(int n)
{
	int i;
	double d;

	d=0;
	for(i=0; i<n; i++) d = d + dist(ps[i],ps[i+1]);
	return(d);
}

void di(point *ps,int n)
{
	double d;

	if(verbose)
	{	d=pathlength(n);
		fprintf(stderr,"Points: %8d; Length %.2lf\n",
			n,d);
	}
}
