// Postscript plotting routines

#include "k3.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>

#define PAGEX 540       /* Page size: 540 default */
#define PAGEY 720       /* Offset of lower left corner: 720 default */
#define CORNX 36
#define CORNY 36
#define LINEWIDTH 1 //.3	 /* Line width in device units */
#define BORDER 0.10     /* Default border percentage */

//defaults
static double linewidth=LINEWIDTH;
double size, xmax, xmin, ymax, ymin;

//homeomorphism from half-plane onto unit disk
point normalize(point p){
	double temp = sqrt(norm(p));
	p.x = p.x/(1+fabs(p.x));
	p.y = p.y/(1+fabs(p.y));
	p.z = p.z/(1+fabs(p.z));

	temp = sqrt(norm(p));
	p.x = p.x*pow(temp, 0);
	p.y = p.y*pow(temp, 0);
	p.z = p.z*pow(temp, 0);

	p.x = p.x* 2;
	p.y = p.y*.5;
	p.z = p.z*1;

	return p;
}


//starts plot
void ps_open(argc,argv)
	int argc;
	char *argv[];
{
	int i;

	printf("%%!\n");
	printf("%%%% Traj Version 2.0\n");
	printf("%%%% Command line:\n");
        printf("%%%% ");
        for(i=0; i<argc; i++)
                printf("%s ",argv[i]);
	printf("\n");
}

// Scale and box
void ps_window(centerx, centery, radius)
	double centerx, centery;
	double radius;
{
	int marginx, marginy;
	double s, sx, sy, tx, ty;

	xmin = centerx - radius;
	xmax = centerx + radius;
	ymin = centery - radius;
	ymax = centery + radius;

// Print header 
	sx = PAGEX/(xmax-xmin);
	sy = PAGEY/(ymax-ymin);
	s  = (sx < sy) ? sx : sy;
	marginx = (PAGEX-s*(xmax-xmin))/2;
	marginy = (PAGEY-s*(ymax-ymin))/2;

	printf("%%%% Plot of stable manifolds\n");
	printf("%%%%\n");
	printf("%%%%BoundingBox:  %d %d %d %d\n",
		CORNX+marginx,
		CORNY+marginy,
		CORNX+PAGEX-marginx,
		CORNY+PAGEY-marginy);
	printf("%d %d translate\n",
		CORNX+marginx,
		CORNY+marginy);
	printf("%.2lf setlinewidth\n",linewidth);
	printf("%%%%\n");
	printf("%%%% Change to complex coordinates\n");
	printf("%.5lf %.5lf scale\n",s,s);
	printf("currentlinewidth %.5lf div setlinewidth\n",s);
	printf("%.5lf %.5lf translate\n",-xmin,-ymin); 
	printf("%%%% Clipping box %%%%\n");
	printf(
"newpath\n %.5lf %.5lf moveto %.5lf %.5lf lineto\n %.5lf %.5lf lineto %.5lf %.5lf lineto\nclosepath clip\n",
		xmin,ymin,xmin,ymax,xmax,ymax,xmax,ymin);
	
	printf("0.9 setgray\n");
	printf("/l {newpath moveto lineto 0 0 0 setrgbcolor stroke 3} def\n"); //black paths 
	printf("/o {newpath moveto lineto 1 .1 0 setrgbcolor stroke} def\n"); //red paths 
	printf("/d {newpath 0 .5 1 setrgbcolor arc fill} def\n"); //blue dots
	printf("/b {newpath .5 .5 .5 setrgbcolor arc fill} def\n"); //grey dots
	printf("/p {newpath 1 .5 0 setrgbcolor arc fill} def\n"); //orange dots

}

// Determine if point is inside window
int inside(p)
	point p;
{
	if (p.x < xmin) return(0);
	if (p.x > xmax) return(0);
	if (p.y < ymin) return(0);
	if (p.y > ymax) return(0);
	return(1);
}

//plot line between p and q
void ps_line(point p,point q, int unravel)
{
	
	if(unravel){
		p = normalize(stereo_proj(rescale(p)));
		q = normalize(stereo_proj(rescale(q)));
		if(inside(p) || inside(q)){
			printf("%.5lf %.5lf %.5lf %.5lf l\n",
			p.x,p.y,q.x,q.y);
		}
		return;
	} else{
		//use to plot z-x axis
		p = rotate(p);
		p = rotate(p);
		q = rotate(q);
		q = rotate(q);
	}

	if(top(p)){
		if(inside(p) || inside(q)){
			//p = rescale(p); //Ethan added for sphere
			//q = rescale(q);
			printf("%.5lf %.5lf %.5lf %.5lf l\n",
			p.x,p.y,q.x,q.y);
		}
	} else if(inside(p) || inside(q)){
		printf("%.5lf %.5lf %.5lf %.5lf o\n",
		p.x,p.y,q.x,q.y);
	}

		
}


void ps_dot(point p, int unravel)
{

	if(unravel){
		p = normalize(stereo_proj(rescale(p)));
		printf("%.5lf %.5lf .01 0 360 d\n", //.015 controls dot size
		p.x,p.y);
		return;
	} else{
		//p = rescale(p); //Ethan added for sphere
		//use to plot z-x axis
		p = rotate(p);
		p = rotate(p);
	//
	}
	if(inside(p) && top(p)){
		printf("%.5lf %.5lf .02 0 360 d\n", 
		p.x,p.y);
	}
	if(inside(p) && !top(p)){
		printf("%.5lf %.5lf .02 0 360 p\n", 
		p.x,p.y);
	}
}

void ps_dot_transparent(point p)/*Ethan added this*/
{
	//use to plot z-x axis
	p = rotate(p);
	p = rotate(p);

	if(inside(p))
		printf("%.5lf %.5lf .003 0 360 b\n", //originally .003
		p.x,p.y);
}

void ps_close()
{
	printf("showpage\n");
}
