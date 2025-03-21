
/* For finding recurrent points with controllable precision */


#include <math.h>
#include "k3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h> //for precision orbits


static double centerx = CENTERX, centery = CENTERY, radius = RADIUS;
static double length=LENGTH, sep = SEP;
static int drawtop = 1, drawbot = 1, verbose=VERBOSE, spin=0, unravel = 0;

/* Very large array to dynamically store paths*/
point ps[PSMAX];

//very important for determine_arrays
int pos_array[PSMAX]; //HOW BIG SHOULD I MAKE THESE??
int above_array[PSMAX];

int orbit_length = 23;
const char *word[PSMAX];
int above = -100000;
int current_pos = -10000;

//new point structure
typedef struct {
    mpf_t x, y, z;
} mpf_point;

void init_point(mpf_point *p) {
    mpf_init(p->x);
    mpf_init(p->y);
    mpf_init(p->z);
}

//function declarations
void involute_exact(mpf_t x, mpf_t y, mpf_t z);
void ix_exact(mpf_point *p);
void iy_exact(mpf_point *p);
void iz_exact(mpf_point *p);
void f_exact(mpf_point *p);
void lift_exact(mpf_point *p,  mpf_t x,  mpf_t y);
void lift_exact_under(mpf_point *p,  mpf_t x,  mpf_t y);
void lift_exact_charts(mpf_point *p,  mpf_t x,  mpf_t y, int i); //lifts with the ith chart
bool unreal_exact( mpf_t x,  mpf_t y);
void lift_exact_y(mpf_point *p,  mpf_t x,  mpf_t z);
void f_c(mpf_t a,  mpf_t b,  int i, int j); //i = incoming chart, j = outgoing chart
void projection_chart(mpf_point p, mpf_t a, mpf_t b, int i);
void dist_2plane(mpf_t result, mpf_t a1, mpf_t b1, mpf_t a2, mpf_t b2);

mpf_t ae, be;

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

	

	mpf_set_default_prec(500); 

    mpf_t result, epsilon, delta, temp;
    mpf_inits(ae, be, result, epsilon, delta, temp, NULL);
	
    mpf_t a_start, b_start;
    mpf_inits(a_start,b_start,NULL);

    mpf_t a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, a9, b9, a10, b10;
    mpf_inits(a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, a9, b9, a10, b10, NULL);


    mpf_set_str(delta, "1e-29", 10); //recurrence requirement //1e-7

    mpf_set_str(ae, "10", 10);
	mpf_set_str(be, "2", 10);

    mpf_set_str(a1, "1.041643093944314148360673792017", 10); 
    mpf_set_str(b1, "1.726895448754858426328854724474", 10); 
    
    mpf_set_str(a2, "-0.439586738044637984442175311821", 10); 
    mpf_set_str(b2, "0.555943953085459715621476373770", 10);
    
    mpf_set_str(a3, "1.111402054756051352317454938205", 10); 
    mpf_set_str(b3, "-0.926435350008842162121508383319", 10);
    
    mpf_set_str(a4, "-0.328869789067645570794396144391", 10); 
    mpf_set_str(b4, "0.867950394543647373540816310310", 10);
    
    mpf_set_str(a5, "1.488818954806569814700993326668", 10); 
    mpf_set_str(b5, "1.004464450964796276608907444033", 10);
    
    mpf_set_str(a6, "0.533829900932504729554816817729", 10); 
    mpf_set_str(b6, "-0.533829900932504729554816817729", 10);
    
    mpf_set_str(a7, "-1.004464450964796276608907444033", 10); 
    mpf_set_str(b7, "-1.488818954806569814700993326668", 10);
    
    mpf_set_str(a8, "-0.867950394543647373540816310310", 10); 
    mpf_set_str(b8, "-0.328869789067645570794396144391", 10);
    
    mpf_set_str(a9, "0.926435350008842162121508383319", 10); 
    mpf_set_str(b9, "-1.111402054756051352317454938205", 10);
    
    mpf_set_str(a10, "0.555943953085459715621476373770", 10); 
    mpf_set_str(b10, "0.439586738044637984442175311821", 10);

    mpf_set(a_start, a1);
    mpf_set(b_start, b1);

    f_c(a1, b1, 6, 4);

    dist_2plane(temp, a1, b1, a2, b2);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success \n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a2, b2, 4, 5);
    
    dist_2plane(temp, a2, b2, a3, b3);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success \n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a3, b3, 5, 6);
    
    dist_2plane(temp, a3, b3, a4, b4);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a4, b4, 6, 5);
    
    dist_2plane(temp, a4, b4, a5, b5);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a5, b5, 5, 5);
    
    dist_2plane(temp, a5, b5, a6, b6);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a6, b6, 5, 5);
    
    dist_2plane(temp, a6, b6, a7, b7);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a7, b7, 5, 1);
    
    dist_2plane(temp, a7, b7, a8, b8);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a8, b8, 1, 5);
    
    dist_2plane(temp, a8, b8, a9, b9);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a9, b9, 5, 3);
    
    dist_2plane(temp, a9, b9, a10, b10);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    }else{
        fprintf(stderr, "fail");
    }
    
    f_c(a10, b10, 3, 6);

    dist_2plane(temp, a10, b10, a_start,b_start);
    if(mpf_cmp(temp, delta) < 0) {
        fprintf(stderr, "success\n");
    } else{
        fprintf(stderr, "fail");
    }
    
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

void involute_exact(mpf_t x, mpf_t y, mpf_t z) {
    mpf_t temp1, temp2, temp3;

    mpf_inits(temp1, temp2, temp3, NULL);

    // temp1 = x * y
    mpf_mul(temp1, x, y);

    // temp2 = 1 + x * x
    mpf_mul(temp2, x, x);
    mpf_add_ui(temp2, temp2, 1);

    // temp3 = 1 + y * y
    mpf_mul(temp3, y, y);
    mpf_add_ui(temp3, temp3, 1);

    // temp2 = temp2 * temp3
    mpf_mul(temp2, temp2, temp3);

    // temp1 = -a * temp1 / temp2
    mpf_mul(temp1, temp1, ae);
    mpf_neg(temp1, temp1);
    mpf_div(temp1, temp1, temp2);

    // result = temp1 - z
    mpf_sub(z, temp1, z);

    mpf_clears(temp1, temp2, temp3, NULL);
}

void ix_exact(mpf_point *p)
{
	involute_exact(p->z,p->y,p->x);
}
void iy_exact(mpf_point *p)
{
	involute_exact(p->z,p->x,p->y);
}
void iz_exact(mpf_point *p)
{
	involute_exact(p->x,p->y,p->z);
}
void f_exact(mpf_point *p)
{
	ix_exact(p);
	iy_exact(p);
	iz_exact(p);
}
void f_c(mpf_t s,  mpf_t t, int i, int j) {
	mpf_point temp;
	init_point(&temp);
	lift_exact_charts(&temp, s, t, i);
	//gmp_fprintf(stderr, "(%.5Ff, %.5Ff, %.5Ff)\n", temp.x, temp.y, temp.z);
	f_exact(&temp);
	//gmp_fprintf(stderr, "after (%.5Ff, %.5Ff, %.5Ff)\n", temp.x, temp.y, temp.z);

	if(j == 1 || j== 4){
		mpf_set(s, temp.x);
		mpf_set(t, temp.y);
	}
	if(j == 2 || j== 5){
		mpf_set(s, temp.x);
		mpf_set(t, temp.z);
	}
	if(j == 3 || j== 6){
		mpf_set(s, temp.y);
		mpf_set(t, temp.z);
	}
	//gmp_fprintf(stderr, "inside (%.5Ff, %.5Ff)\n", s,t);
	mpf_clears(temp.x, temp.y, temp.z, NULL); //this right?
}



void dist_2plane(mpf_t result, mpf_t a1, mpf_t b1, mpf_t a2, mpf_t b2) {
    mpf_t dx, dy, sum;

    mpf_inits(dx, dy, sum, NULL);

    mpf_sub(dx, a1, a2);
    mpf_sub(dy, b1, b2);

    mpf_mul(dx, dx, dx);
    mpf_mul(dy, dy, dy);

    mpf_add(sum, dx, dy);

    mpf_sqrt(result, sum);

    mpf_clears(dx, dy, sum, NULL);
}

void projection_chart(mpf_point p, mpf_t a, mpf_t b, int i) {
   if (i == 1 || i == 4) {
		mpf_set(a, p.x);
		mpf_set(b, p.y);
   }
   if (i == 2 || i == 5) {
		mpf_set(a, p.x);
		mpf_set(b, p.z);
   }
   if (i == 3 || i == 6) {
		mpf_set(a, p.y);
		mpf_set(b, p.z);
   }
}

bool is_recurrent_exact_specific(mpf_t a, mpf_t b, mpf_t delta) /*returns true of p is delta, k-recurrent for k <= 22*/
{
	//fprintf(stderr, "ummm");
	mpf_t temp2;
	mpf_t a_start;
	mpf_t b_start;
	mpf_init(a_start);
	mpf_init(b_start);
	mpf_init(temp2);

	mpf_set(a_start, a);
	mpf_set(b_start, b);
	
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	f_c(a, b, 6, 4);	
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	f_c(a, b, 4, 5);	
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);

	f_c(a, b, 5, 6);
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	
	f_c(a, b, 6, 5);
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	
	f_c(a, b, 5, 5);
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	
	f_c(a, b, 5, 5);
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	
	f_c(a, b, 5, 1);	
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);

	f_c(a, b, 1, 5);
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	
	f_c(a, b, 5, 3);
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	
	f_c(a, b, 3, 6);
	//gmp_fprintf(stderr, "(%.20Ff, %.20Ff)\n", a,b);
	
	
	mpf_clears(temp2, NULL);
    return false;
}


bool unreal_exact( mpf_t x,  mpf_t y) {
	mpf_t qa, qb, qc, temp1, temp2, sqrt_term;

    mpf_inits(qa, qb, qc, temp1, temp2, sqrt_term, NULL);

    // qa = (1 + x*x) * (1 + y*y)
    mpf_mul(temp1, x, x);
    mpf_add_ui(temp1, temp1, 1);
    mpf_mul(temp2, y, y);
    mpf_add_ui(temp2, temp2, 1);
    mpf_mul(qa, temp1, temp2);

    // qb = a * x * y
    mpf_mul(qb, ae, x);
    mpf_mul(qb, qb, y);

    // qc = qa - b
    mpf_sub(qc, qa, be);


    // p.z = (-qb + sqrt(qb*qb - 4*qa*qc)) / (2*qa)
    mpf_mul(temp1, qb, qb); // temp1 = qb * qb
    mpf_mul_ui(temp2, qa, 4); // temp2 = 4 * qa
    mpf_mul(temp2, temp2, qc); // temp2 = 4 * qa * qc
    mpf_sub(sqrt_term, temp1, temp2); // sqrt_term = qb*qb - 4*qa*qc
	if(mpf_cmp_ui(sqrt_term, 0) < 0){
		return true;
	} else {
		return false;
	}
}
void lift_exact(mpf_point *p,  mpf_t x,  mpf_t y) { //MAKE SURE TO THROW ERROR IF NOT REAL!

    mpf_t qa, qb, qc, temp1, temp2, sqrt_term;

    mpf_inits(qa, qb, qc, temp1, temp2, sqrt_term, NULL);

    // qa = (1 + x*x) * (1 + y*y)
    mpf_mul(temp1, x, x);
    mpf_add_ui(temp1, temp1, 1);
    mpf_mul(temp2, y, y);
    mpf_add_ui(temp2, temp2, 1);
    mpf_mul(qa, temp1, temp2);

    // qb = a * x * y
    mpf_mul(qb, ae, x);
    mpf_mul(qb, qb, y);

    // qc = qa - b
    mpf_sub(qc, qa, be);


    // p.x = x
    mpf_set(p->x, x);

    // p.y = y
    mpf_set(p->y, y);

    // p.z = (-qb + sqrt(qb*qb - 4*qa*qc)) / (2*qa)
    mpf_mul(temp1, qb, qb); // temp1 = qb * qb
    mpf_mul_ui(temp2, qa, 4); // temp2 = 4 * qa

    mpf_mul(temp2, temp2, qc); // temp2 = 4 * qa * qc
    mpf_sub(sqrt_term, temp1, temp2); // sqrt_term = qb*qb - 4*qa*qc

    mpf_sqrt(sqrt_term, sqrt_term); // sqrt_term = sqrt(qb*qb - 4*qa*qc)

    mpf_neg(qb, qb); // qb = -qb
    mpf_add(temp1, qb, sqrt_term); // temp1 = -qb + sqrt(qb*qb - 4*qa*qc)
    mpf_mul_ui(temp2, qa, 2); // temp2 = 2 * qa
    mpf_div(p->z, temp1, temp2); // p.z = temp1 / temp2


    mpf_clears(qa, qb, qc, temp1, temp2, sqrt_term, NULL);
}

void lift_exact_under(mpf_point *p,  mpf_t x,  mpf_t y) { //MAKE SURE TO THROW ERROR IF NOT REAL!

    mpf_t qa, qb, qc, temp1, temp2, sqrt_term;

    mpf_inits(qa, qb, qc, temp1, temp2, sqrt_term, NULL);

    // qa = (1 + x*x) * (1 + y*y)
    mpf_mul(temp1, x, x);
    mpf_add_ui(temp1, temp1, 1);
    mpf_mul(temp2, y, y);
    mpf_add_ui(temp2, temp2, 1);
    mpf_mul(qa, temp1, temp2);

    // qb = a * x * y
    mpf_mul(qb, ae, x);
    mpf_mul(qb, qb, y);

    // qc = qa - b
    mpf_sub(qc, qa, be);


    // p.x = x
    mpf_set(p->x, x);

    // p.y = y
    mpf_set(p->y, y);

    // p.z = (-qb + sqrt(qb*qb - 4*qa*qc)) / (2*qa)
    mpf_mul(temp1, qb, qb); // temp1 = qb * qb
    mpf_mul_ui(temp2, qa, 4); // temp2 = 4 * qa

    mpf_mul(temp2, temp2, qc); // temp2 = 4 * qa * qc
    mpf_sub(sqrt_term, temp1, temp2); // sqrt_term = qb*qb - 4*qa*qc

    mpf_sqrt(sqrt_term, sqrt_term); // sqrt_term = sqrt(qb*qb - 4*qa*qc)

    mpf_neg(qb, qb); // qb = -qb
    mpf_sub(temp1, qb, sqrt_term); // temp1 = -qb - sqrt(qb*qb - 4*qa*qc)
    mpf_mul_ui(temp2, qa, 2); // temp2 = 2 * qa
    mpf_div(p->z, temp1, temp2); // p.z = temp1 / temp2
    mpf_clears(qa, qb, qc, temp1, temp2, sqrt_term, NULL);
}


void lift_exact_charts(mpf_point *p,  mpf_t a,  mpf_t b, int i) {
	mpf_t temp0;
	mpf_init(temp0);

	if(i == 1){
		lift_exact(p, a, b);
	}
	if(i == 2){
		lift_exact(p, a, b);
		mpf_set(temp0, p->y);
		mpf_set(p->y, p->z);
		mpf_set(p->z, temp0);
	}
	if(i == 3){
		lift_exact(p, a, b);
		mpf_set(temp0, p->z);
		mpf_set(p->z, p->y);
		mpf_set(p->y, p->x);
		mpf_set(p->x, temp0);
	}
	if(i == 4){
		lift_exact_under(p, a, b);
	}
	if(i == 5){
		lift_exact_under(p, a, b);
		mpf_set(temp0, p->y);
		mpf_set(p->y, p->z);
		mpf_set(p->z, temp0);
	}
	if(i == 6){
		lift_exact_under(p, a, b);
		mpf_set(temp0, p->z);
		mpf_set(p->z, p->y);
		mpf_set(p->y, p->x);
		mpf_set(p->x, temp0);
	}
	mpf_clears(temp0, NULL);
}





//searches along line x = -z
void search_near_exact_line(mpf_point *q, mpf_t epsilon, mpf_t delta, int N) {
    mpf_t x, z, temp; //x_t & y_t used for substitution
	mpf_t a,b;
	mpf_inits(a,b, NULL);

    mpf_point p;
    int length = 2 * N;

    mpf_inits(x, z, temp, NULL);
    init_point(&p);

    // x = q.x - N * epsilon
    mpf_mul_ui(temp, epsilon, N);
    mpf_sub(x, q->x, temp);
    // y = q.y - N * epsilon
    mpf_add(z, q->z, temp);
    for (int j = 0; j <= length; j++) {
        if (!unreal_exact(x, z)) { //check unreal
            
            lift_exact_charts(&p, x, z, 2); //add lift
			projection_chart(p, a, b, 6);
			/*using f, not f_c
			if (is_recurrent_exact(p, delta)) {
                gmp_fprintf(stderr, "Success at (%.20Ff, %.20Ff, %.20Ff)\n", p.x, p.y, p.z);
            } */

			if (is_recurrent_exact_specific(a,b, delta)) {
                gmp_fprintf(stderr, "Success at (%.50Ff, %.50Ff, %.50Ff)\n", p.x, p.y, p.z);
            }
			projection_chart(p, a, b, 6);

			/*
			if (is_recurrent_exact(p, delta)) {			//this bit checks bottoms as well:
                gmp_fprintf(stderr, "Success at (%.20Ff, %.20Ff, %.20Ff)\n", p.x, p.y, p.z);
            } */

        }
            // x = x + epsilon
        mpf_add(x, x, epsilon);
		mpf_sub(z, z, epsilon);
	}
        // x = q.x - N * epsilon
    mpf_clears(x, z, temp, NULL);
    //clear_point(&p);
}

