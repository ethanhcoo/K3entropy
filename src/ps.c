
/* Postscript plotting routines */

#include "k3.h"
#include "math.h"

#include <stdio.h>

#define PAGEX 540       /* Page size: was 540 */
#define PAGEY 720       /* Offset of lower left corner */
#define CORNX 36
#define CORNY 36
#define LINEWIDTH .02// was 0.20  /* Line width in device units */
#define BORDER 0.10     /* Default border percentage */

/* Defaults */
static double linewidth=LINEWIDTH;
double size, xmax, xmin, ymax, ymin;

point normalize(point p){
	double temp = sqrt(norm(p));
	p.x = p.x/(1+temp);
	p.y = p.y/(1+temp);
	p.z = p.z/(1 + temp);

	temp = sqrt(norm(p));
	p.x = p.x*pow(temp, 0);
	p.y = p.y*pow(temp, 0);
	p.z = p.z*pow(temp, 0);

	p.x = p.x*2;
	p.y = p.y*.5;
	p.z = p.z*1;

	return p;
}


/* Start plot */
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

/* Scale and box */
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

/* Print header */
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
	printf("/l {newpath moveto lineto 0 0 0 setrgbcolor stroke 3} def\n"); //stroke changes line thickness
	printf("/o {newpath moveto lineto 1 0 0 setrgbcolor stroke} def\n"); /*Ethan added this line for orange paths*/
	printf("/d {newpath 0 0 1 setrgbcolor arc fill} def\n"); /*Ethan added this line for orange paths*/
	printf("/b {newpath .5 .5 .5 setrgbcolor arc fill} def\n"); /*Ethan added this line for backgroung surface*/
	printf("/p {newpath 0 1 0 setrgbcolor arc fill} def\n"); /*Ethan added this line for orange paths*/

}

/* Determine if point is inside window */
int inside(p)
	point p;
{
	if (p.x < xmin) return(0);
	if (p.x > xmax) return(0);
	if (p.y < ymin) return(0);
	if (p.y > ymax) return(0);
	return(1);
}
	
void ps_line(point p,point q, int unravel)
{
	if(unravel){
		p = normalize(stereo_proj(rescale(p)));
		q = normalize(stereo_proj(rescale(q)));
		if(inside(p) || inside(q)){
			printf("%.5lf %.5lf %.5lf %.5lf o\n",
			p.x,p.y,q.x,q.y);
		}
		return;
	}
	

	if(inside(p) || inside(q)){
		//p = rescale(p); //Ethan added for sphere
		//q = rescale(q);
		printf("%.5lf %.5lf %.5lf %.5lf l\n",
		p.x,p.y,q.x,q.y);
	}
}

void ps_line_color(point p, point q, int unravel) /*Ethan added this. It draws a colored segment*/
{
	if(unravel){
		p = normalize(stereo_proj(rescale(p)));
		q = normalize(stereo_proj(rescale(q)));
	} else{
		//p = rescale(p); //Ethan added for sphere
		//q = rescale(q); //Ethan added for sphere
	}
	
	if(inside(p) || inside(q)){
		
			printf("%.5lf %.5lf %.5lf %.5lf o\n",
			p.x,p.y,q.x,q.y);
		}
}

void ps_dot(point p, int unravel)/*Ethan added this. It draws a dot*/
{
	if(unravel){
		p = normalize(stereo_proj(rescale(p)));
		printf("%.5lf %.5lf .005 0 360 d\n", //ethan added
		p.x,p.y);
		return;
	} else{
		//p = rescale(p); //Ethan added for sphere
	}
	if(inside(p) && top(p)){
		printf("%.5lf %.5lf .03 0 360 d\n", //originally .005, changed to 5e-7
		p.x,p.y);
	}
	if(inside(p) && !top(p)){
		printf("%.5lf %.5lf .03 0 360 p\n", /*fist number controls dot size*/
		p.x,p.y);
	}
}

void ps_dot_transparent(point p)/*Ethan added this*/
{
	if(inside(p))
		printf("%.5lf %.5lf .003 0 360 b\n", //originally .003
		p.x,p.y);
}

void ps_close()
{
	printf("showpage\n");
}
