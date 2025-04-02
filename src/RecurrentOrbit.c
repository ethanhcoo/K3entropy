
/* A particular pseudo-orbit */
#include <math.h>
#include "k3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
//#include <gmp.h> //for precision orbits
#include <stdarg.h>
#include <mpfr.h> 
//#include <mpf2mpfr.h>


//new point structure
typedef struct {
    mpfr_t x, y, z;
} mpfr_point;

void init_point(mpfr_point *p) {
    mpfr_init(p->x);
    mpfr_init(p->y);
    mpfr_init(p->z);
}

//function declarations
void involute_exact(mpfr_t x, mpfr_t y, mpfr_t z);
void ix_exact(mpfr_point *p);
void iy_exact(mpfr_point *p);
void iz_exact(mpfr_point *p);
void f_exact(mpfr_point *p);
void lift_exact(mpfr_point *p,  mpfr_t x,  mpfr_t y);
void lift_exact_under(mpfr_point *p,  mpfr_t x,  mpfr_t y);
void lift_exact_charts(mpfr_point *p,  mpfr_t x,  mpfr_t y, int i); //lifts with the ith chart
bool unreal_exact( mpfr_t x,  mpfr_t y);
void lift_exact_y(mpfr_point *p,  mpfr_t x,  mpfr_t z);
void f_c(mpfr_t a,  mpfr_t b,  int i, int j); //i = incoming chart, j = outgoing chart
void projection_chart(mpfr_point p, mpfr_t a, mpfr_t b, int i);
void dist_2plane(mpfr_t result, mpfr_t a1, mpfr_t b1, mpfr_t a2, mpfr_t b2);

void multiplyMatricesExact(mpfr_t **firstMatrix, mpfr_t **secondMatrix, mpfr_t **result, int ROWS1, int COLS1, int ROWS2, int COLS2);
void multiplyMatricesAuxExact(mpfr_t **firstMatrix, mpfr_t **secondMatrix, mpfr_t **result, int ROWS1, int COLS1, int ROWS2, int COLS2);
void printMatrixExact(mpfr_t **matrix, int ROWS1, int COLS1);
void DixExact(mpfr_point p, mpfr_t **matrix);
void DiyExact(mpfr_point p, mpfr_t **matrix);
void DizExact(mpfr_point p, mpfr_t **matrix);
void DExact(mpfr_t result, mpfr_t x, mpfr_t y);
void DxExact(mpfr_t result, mpfr_t x, mpfr_t y);
void DyExact(mpfr_t result, mpfr_t x, mpfr_t y);
void DxxExact(mpfr_t result, mpfr_t x, mpfr_t y);
void DyyExact(mpfr_t result, mpfr_t x, mpfr_t y);
void DxyExact(mpfr_t result, mpfr_t x, mpfr_t y);
void pxExact(mpfr_t result, mpfr_t s, mpfr_t t);
void pyExact(mpfr_t result, mpfr_t s, mpfr_t t);
double Dx(double x, double y) {
    double result = 200 * x * y * y + 16 * x * (1 + y * y) - 16 * x * (1 + x * x) * (1 + y * y) * (1 + y * y);
    return result;
}

void alpha_xExact(mpfr_t result, const mpfr_t x, const mpfr_t y);
void alpha_yExact(mpfr_t result, const mpfr_t x, const mpfr_t y);

mpfr_t ae, be;

int main(int argc, char *argv[]) {

    mpfr_set_default_rounding_mode(MPFR_RNDZ);

    mpfr_set_default_prec(5000); // Set precision to 500 bits
    mpfr_set_emin (-1073); mpfr_set_emax (1024); //Sets the exponent range to [-1073,1024]

    mpfr_t** matrix = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** coordmatrix = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** matrix_temp = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** temp1 = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** temp2 = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** temp3 = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** temp4 = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** temp5 = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** x_mat = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** y_mat = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));
    mpfr_t** z_mat = (mpfr_t**)malloc(3 * sizeof(mpfr_t*));

    for (int i = 0; i < 3; ++i) {
        matrix[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        coordmatrix[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        matrix_temp[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        temp1[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        temp2[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        temp3[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        temp4[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        temp5[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        x_mat[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        y_mat[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));
        z_mat[i] = (mpfr_t*)malloc(3 * sizeof(mpfr_t));

        for (int j = 0; j < 3; ++j) {
            mpfr_init(matrix[i][j]);
            mpfr_init(coordmatrix[i][j]);
            mpfr_init(matrix_temp[i][j]);
            mpfr_init(temp1[i][j]);
            mpfr_init(temp2[i][j]);
            mpfr_init(temp3[i][j]);
            mpfr_init(temp4[i][j]);
            mpfr_init(temp5[i][j]);
            mpfr_init(x_mat[i][j]);
            mpfr_init(y_mat[i][j]);
            mpfr_init(z_mat[i][j]);

            if (i == j) {
                mpfr_set_ui(matrix[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(coordmatrix[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(matrix_temp[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(temp1[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(temp2[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(temp3[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(temp4[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(temp5[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(x_mat[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(y_mat[i][j], 1, MPFR_RNDZ);
                mpfr_set_ui(z_mat[i][j], 1, MPFR_RNDZ);
            } else {
                mpfr_set_ui(matrix[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(coordmatrix[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(matrix_temp[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(temp1[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(temp2[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(temp3[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(temp4[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(temp5[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(x_mat[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(y_mat[i][j], 0, MPFR_RNDZ);
                mpfr_set_ui(z_mat[i][j], 0, MPFR_RNDZ);
            }
        }
    }

    

    mpfr_t result, epsilon, delta, temp;
    mpfr_inits(ae, be, result, epsilon, delta, temp, NULL);
	
   
    mpfr_t a_start, b_start;
    mpfr_t a_temp, b_temp;
    mpfr_inits(a_start,b_start,NULL);
    mpfr_inits(a_temp,b_temp,NULL);


    //mpfr_t a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, a9, b9, a10, b10;
   // mpfr_inits(a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, a9, b9, a10, b10, NULL);


    mpfr_set_str(delta, "1e-20", 10, MPFR_RNDZ); //recurrence requirement //1e-7

    mpfr_set_str(ae, "10", 10, MPFR_RNDZ);
	mpfr_set_str(be, "2", 10, MPFR_RNDZ);

    mpfr_t a[10];
    for (int i = 0; i < 10; ++i) {
        mpfr_init(a[i]);
    }

    mpfr_set_str(a[0], "1.041643093944314148360673792017", 10, MPFR_RNDZ);
    mpfr_set_str(a[1], "-0.439586738044637984442175311821", 10, MPFR_RNDZ);
    mpfr_set_str(a[2], "1.111402054756051352317454938205", 10, MPFR_RNDZ);
    mpfr_set_str(a[3], "-0.328869789067645570794396144391", 10, MPFR_RNDZ);
    mpfr_set_str(a[4], "1.488818954806569814700993326668", 10, MPFR_RNDZ);
    mpfr_set_str(a[5], "0.533829900932504729554816817729", 10, MPFR_RNDZ);
    mpfr_set_str(a[6], "-1.004464450964796276608907444033", 10, MPFR_RNDZ);
    mpfr_set_str(a[7], "-0.867950394543647373540816310310", 10, MPFR_RNDZ);
    mpfr_set_str(a[8], "0.926435350008842162121508383319", 10, MPFR_RNDZ);
    mpfr_set_str(a[9], "0.555943953085459715621476373770", 10, MPFR_RNDZ);

    mpfr_t b[10];
    for (int i = 0; i < 10; ++i) {
        mpfr_init(b[i]);
    }

    mpfr_set_str(b[0], "1.726895448754858426328854724474", 10, MPFR_RNDZ);
    mpfr_set_str(b[1], "0.555943953085459715621476373770", 10, MPFR_RNDZ);
    mpfr_set_str(b[2], "-0.926435350008842162121508383319", 10, MPFR_RNDZ);
    mpfr_set_str(b[3], "0.867950394543647373540816310310", 10, MPFR_RNDZ);
    mpfr_set_str(b[4], "1.004464450964796276608907444033", 10, MPFR_RNDZ);
    mpfr_set_str(b[5], "-0.533829900932504729554816817729", 10, MPFR_RNDZ);
    mpfr_set_str(b[6], "-1.488818954806569814700993326668", 10, MPFR_RNDZ);
    mpfr_set_str(b[7], "-0.328869789067645570794396144391", 10, MPFR_RNDZ);
    mpfr_set_str(b[8], "-1.111402054756051352317454938205", 10, MPFR_RNDZ);
    mpfr_set_str(b[9], "0.439586738044637984442175311821", 10, MPFR_RNDZ);


    int charts[10] = {6, 4, 5, 6, 5, 5, 5, 1, 5, 3};


    mpfr_set(a_start, a[0], MPFR_RNDZ);
    mpfr_set(b_start, b[0], MPFR_RNDZ);

    for (int i = 0; i < 10; i++) {
        mpfr_set (a_temp, a[i], MPFR_RNDZ);
        mpfr_set (b_temp, b[i], MPFR_RNDZ);
        f_c(a_temp, b_temp, charts[i], charts[(i + 1) % 10]);
        dist_2plane(temp, a_temp, b_temp, a[(i + 1) % 10], b[(i + 1) % 10]);
        if (mpfr_cmp(temp, delta) < 0) {
            fprintf(stderr, "success \n");
        } else {
            fprintf(stderr, "fail \n");
        }
    
        char a_str[64], b_str[64];
        mpfr_sprintf(a_str, "%.15Rf", a[i]);
        mpfr_sprintf(b_str, "%.15Rf", b[i]);
    
        fprintf(stderr, "(x,y) is (%s, %s) \n", a_str, b_str);
    
        mpfr_t result;
        mpfr_init(result);
    
        DExact(result, a[i], b[i]);
        mpfr_sprintf(a_str, "%.15Rf", result);
        fprintf(stderr, "DExact(x,y) is %s \n", a_str);
    
        DxExact(result, a[i], b[i]);
        mpfr_sprintf(a_str, "%.15Rf", result);
        fprintf(stderr, "DxExact(x,y) is %s \n", a_str);
    
        DyExact(result, a[i], b[i]);
        mpfr_sprintf(a_str, "%.15Rf", result);
        fprintf(stderr, "DyExact(x,y) is %s \n", a_str);
    
        DxxExact(result, a[i], b[i]);
        mpfr_sprintf(a_str, "%.15Rf", result);
        fprintf(stderr, "DxxExact(x,y) is %s \n", a_str);
    
        DyyExact(result, a[i], b[i]);
        mpfr_sprintf(a_str, "%.15Rf", result);
        fprintf(stderr, "DyyExact(x,y) is %s \n", a_str);
    
        DxyExact(result, a[i], b[i]);
        mpfr_sprintf(a_str, "%.15Rf", result);
        fprintf(stderr, "DxyExact(x,y) is %s \n", a_str);
    
        mpfr_clear(result);
    
        int k = 0;
        
        for (int j = 0; j < 3; j++) {
            if (j == 3 - charts[i]) {
                pxExact(temp4[j][0], a[i], b[i]);
                pyExact(temp4[j][1], a[i], b[i]);
            } else if (j == 6 - charts[i]) {
                mpfr_neg(temp4[j][0], a[i], MPFR_RNDZ);
                pxExact(temp4[j][0], temp4[j][0], b[i]);
                mpfr_neg(temp4[j][1], b[i], MPFR_RNDZ);
                pyExact(temp4[j][1], temp4[j][1], b[i]);
            } else {
                if (k == 0) {
                    mpfr_set_ui(temp4[j][0], 1, MPFR_RNDZ);
                    mpfr_set_ui(temp4[j][1], 0, MPFR_RNDZ);
                }
                if (k == 1) {
                    mpfr_set_ui(temp4[j][0], 0, MPFR_RNDZ);
                    mpfr_set_ui(temp4[j][1], 1, MPFR_RNDZ);
                }
                if (k == 2) {
                    fprintf(stderr, "error");
                    abort();
                }
                k++;
            }
        }
        k = 0;
    
        
        for (int j = 0; j < 3; j++) {
            if (j == 3 - charts[(i+1)%10]) {
                mpfr_set_ui(temp5[0][j], 0, MPFR_RNDZ);
                mpfr_set_ui(temp5[1][j], 0, MPFR_RNDZ);
            } else if (j == 6 - charts[(i+1)%10]) {
                mpfr_set_ui(temp5[0][j], 0, MPFR_RNDZ);
                mpfr_set_ui(temp5[1][j], 0, MPFR_RNDZ);
            } else {
                if (k == 0) {
                    mpfr_set_ui(temp5[0][j], 1, MPFR_RNDZ);
                    mpfr_set_ui(temp5[1][j], 0, MPFR_RNDZ);
                }
                if (k == 1) {
                    mpfr_set_ui(temp5[0][j], 0, MPFR_RNDZ);
                    mpfr_set_ui(temp5[1][j], 1, MPFR_RNDZ);
                }
                if (k == 2) {
                    fprintf(stderr, "error");
                    abort();
                }
                k++;
            }
        }
    
        fprintf(stderr, "chart is %d \n", charts[i]);
        
        mpfr_point p;
	    init_point(&p);
	    lift_exact_charts(&p, a[i], b[i], charts[i]);
        mpfr_sprintf(a_str, "%.15Rf", p.y);
        mpfr_sprintf(b_str, "%.15Rf", p.z);
    
        fprintf(stderr, "(x,y) is (%s, %s) \n", a_str, b_str);

        DixExact(p, x_mat);
        ix_exact(&p);
        DiyExact(p, y_mat);

        iy_exact(&p);
        DizExact(p, z_mat);

        iz_exact(&p);
    
        multiplyMatricesExact(y_mat, x_mat, temp2, 3, 3, 3, 3);
        multiplyMatricesExact(z_mat, temp2, temp3, 3, 3, 3, 3);
        printMatrixExact(temp3, 3, 3);


    
        printMatrixExact(temp4, 3, 2);

        multiplyMatricesAuxExact(temp5, temp3, temp1, 2, 3, 3, 3);
        multiplyMatricesAuxExact(temp3, temp4, temp1, 2, 3, 3, 2);
        fprintf(stderr, "(Df^c}_{x_%d^c} =  \n", i);
        printMatrixExact(temp4, 2, 2);
        multiplyMatricesAuxExact(temp4, coordmatrix, temp1, 2, 2, 2, 2);
    }
	exit(0);

}


void involute_exact(mpfr_t x, mpfr_t y, mpfr_t z) {
    mpfr_t temp1, temp2, temp3;

    mpfr_inits(temp1, temp2, temp3, NULL);

    // temp1 = x * y
    mpfr_mul(temp1, x, y, MPFR_RNDZ);

    // temp2 = 1 + x * x
    mpfr_mul(temp2, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);

    // temp3 = 1 + y * y
    mpfr_mul(temp3, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp3, temp3, 1, MPFR_RNDZ);

    // temp2 = temp2 * temp3
    mpfr_mul(temp2, temp2, temp3, MPFR_RNDZ);

    // temp1 = -ae * temp1 / temp2
    mpfr_mul(temp1, temp1, ae, MPFR_RNDZ);
    mpfr_neg(temp1, temp1, MPFR_RNDZ);
    mpfr_div(temp1, temp1, temp2, MPFR_RNDZ);

    // result = temp1 - z
    mpfr_sub(z, temp1, z, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, NULL);
}

void ix_exact(mpfr_point *p)
{
	involute_exact(p->z,p->y,p->x);
}
void iy_exact(mpfr_point *p)
{
	involute_exact(p->z,p->x,p->y);
}
void iz_exact(mpfr_point *p)
{
	involute_exact(p->x,p->y,p->z);
}
void f_exact(mpfr_point *p)
{
	ix_exact(p);
	iy_exact(p);
	iz_exact(p);
}

void f_c(mpfr_t s,  mpfr_t t, int i, int j) {
	mpfr_point temp;
	init_point(&temp);
	lift_exact_charts(&temp, s, t, i);
	//gmp_fprintf(stderr, "(%.5Ff, %.5Ff, %.5Ff)\n", temp.x, temp.y, temp.z);
	f_exact(&temp);
	//gmp_fprintf(stderr, "after (%.5Ff, %.5Ff, %.5Ff)\n", temp.x, temp.y, temp.z);

	if(j == 1 || j== 4){
		mpfr_set(s, temp.x,MPFR_RNDZ);
		mpfr_set(t, temp.y,MPFR_RNDZ);
	}
	if(j == 2 || j== 5){
		mpfr_set(s, temp.x,MPFR_RNDZ);
		mpfr_set(t, temp.z,MPFR_RNDZ);
	}
	if(j == 3 || j== 6){
		mpfr_set(s, temp.y,MPFR_RNDZ);
		mpfr_set(t, temp.z,MPFR_RNDZ);
	}
	//gmp_fprintf(stderr, "inside (%.5Ff, %.5Ff)\n", s,t);
	mpfr_clears(temp.x, temp.y, temp.z, NULL); //this right?
}



void dist_2plane(mpfr_t result, mpfr_t a1, mpfr_t b1, mpfr_t a2, mpfr_t b2) {
    mpfr_t dx, dy, sum;

    mpfr_inits(dx, dy, sum, NULL);

    mpfr_sub(dx, a1, a2, MPFR_RNDZ);
    mpfr_sub(dy, b1, b2, MPFR_RNDZ);

    mpfr_mul(dx, dx, dx, MPFR_RNDZ);
    mpfr_mul(dy, dy, dy, MPFR_RNDZ);

    mpfr_add(sum, dx, dy, MPFR_RNDZ);

    mpfr_sqrt(result, sum, MPFR_RNDZ);

    mpfr_clears(dx, dy, sum, NULL);
}



bool unreal_exact( mpfr_t x,  mpfr_t y) {
	mpfr_t qa, qb, qc, temp1, temp2, sqrt_term;

    mpfr_inits(qa, qb, qc, temp1, temp2, sqrt_term, NULL);

    // qa = (1 + x*x) * (1 + y*y)
    mpfr_mul(temp1, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp1, temp1, 1, MPFR_RNDZ);
    mpfr_mul(temp2, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul(qa, temp1, temp2, MPFR_RNDZ);

    // qb = a * x * y
    mpfr_mul(qb, ae, x, MPFR_RNDZ);
    mpfr_mul(qb, qb, y, MPFR_RNDZ);

    // qc = qa - b
    mpfr_sub(qc, qa, be, MPFR_RNDZ);


    // p.z = (-qb + sqrt(qb*qb - 4*qa*qc)) / (2*qa)
    mpfr_mul(temp1, qb, qb, MPFR_RNDZ); // temp1 = qb * qb
    mpfr_mul_ui(temp2, qa, 4, MPFR_RNDZ); // temp2 = 4 * qa
    mpfr_mul(temp2, temp2, qc, MPFR_RNDZ); // temp2 = 4 * qa * qc
    mpfr_sub(sqrt_term, temp1, temp2, MPFR_RNDZ); // sqrt_term = qb*qb - 4*qa*qc
	if(mpfr_cmp_ui(sqrt_term, 0) < 0){
		return true;
	} else {
		return false;
	}
}

void lift_exact(mpfr_point *p,  mpfr_t x,  mpfr_t y) { //MAKE SURE TO THROW ERROR IF NOT REAL!

    mpfr_t qa, qb, qc, temp1, temp2, sqrt_term;

    mpfr_inits(qa, qb, qc, temp1, temp2, sqrt_term, NULL);

    // qa = (1 + x*x) * (1 + y*y)
    mpfr_mul(temp1, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp1, temp1, 1, MPFR_RNDZ);
    mpfr_mul(temp2, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul(qa, temp1, temp2, MPFR_RNDZ);

    // qb = a * x * y
    mpfr_mul(qb, ae, x, MPFR_RNDZ);
    mpfr_mul(qb, qb, y, MPFR_RNDZ);

    // qc = qa - b
    mpfr_sub(qc, qa, be, MPFR_RNDZ);


    // p.x = x
    mpfr_set(p->x, x, MPFR_RNDZ);

    // p.y = y
    mpfr_set(p->y, y, MPFR_RNDZ);

    // p.z = (-qb + sqrt(qb*qb - 4*qa*qc)) / (2*qa)
    mpfr_mul(temp1, qb, qb, MPFR_RNDZ); // temp1 = qb * qb
    mpfr_mul_ui(temp2, qa, 4, MPFR_RNDZ); // temp2 = 4 * qa

    mpfr_mul(temp2, temp2, qc, MPFR_RNDZ); // temp2 = 4 * qa * qc
    mpfr_sub(sqrt_term, temp1, temp2, MPFR_RNDZ); // sqrt_term = qb*qb - 4*qa*qc

    mpfr_sqrt(sqrt_term, sqrt_term, MPFR_RNDZ); // sqrt_term = sqrt(qb*qb - 4*qa*qc)

    mpfr_neg(qb, qb, MPFR_RNDZ); // qb = -qb
    mpfr_add(temp1, qb, sqrt_term, MPFR_RNDZ); // temp1 = -qb + sqrt(qb*qb - 4*qa*qc)
    mpfr_mul_ui(temp2, qa, 2, MPFR_RNDZ); // temp2 = 2 * qa
    mpfr_div(p->z, temp1, temp2, MPFR_RNDZ); // p.z = temp1 / temp2


    mpfr_clears(qa, qb, qc, temp1, temp2, sqrt_term, NULL);
}

void lift_exact_under(mpfr_point *p,  mpfr_t x,  mpfr_t y) { //MAKE SURE TO THROW ERROR IF NOT REAL!

    mpfr_t qa, qb, qc, temp1, temp2, sqrt_term;

    mpfr_inits(qa, qb, qc, temp1, temp2, sqrt_term, NULL);

    // qa = (1 + x*x) * (1 + y*y)
    mpfr_mul(temp1, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp1, temp1, 1, MPFR_RNDZ);
    mpfr_mul(temp2, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul(qa, temp1, temp2, MPFR_RNDZ);

    // qb = a * x * y
    mpfr_mul(qb, ae, x, MPFR_RNDZ);
    mpfr_mul(qb, qb, y, MPFR_RNDZ);

    // qc = qa - b
    mpfr_sub(qc, qa, be, MPFR_RNDZ);


    // p.x = x
    mpfr_set(p->x, x, MPFR_RNDZ);

    // p.y = y
    mpfr_set(p->y, y, MPFR_RNDZ);

    // p.z = (-qb + sqrt(qb*qb - 4*qa*qc)) / (2*qa)
    mpfr_mul(temp1, qb, qb, MPFR_RNDZ); // temp1 = qb * qb
    mpfr_mul_ui(temp2, qa, 4, MPFR_RNDZ); // temp2 = 4 * qa

    mpfr_mul(temp2, temp2, qc, MPFR_RNDZ); // temp2 = 4 * qa * qc
    mpfr_sub(sqrt_term, temp1, temp2, MPFR_RNDZ); // sqrt_term = qb*qb - 4*qa*qc

    mpfr_sqrt(sqrt_term, sqrt_term, MPFR_RNDZ); // sqrt_term = sqrt(qb*qb - 4*qa*qc)

    mpfr_neg(qb, qb, MPFR_RNDZ); // qb = -qb
    mpfr_sub(temp1, qb, sqrt_term, MPFR_RNDZ); // temp1 = -qb - sqrt(qb*qb - 4*qa*qc)
    mpfr_mul_ui(temp2, qa, 2, MPFR_RNDZ); // temp2 = 2 * qa
    mpfr_div(p->z, temp1, temp2, MPFR_RNDZ); // p.z = temp1 / temp2
    mpfr_clears(qa, qb, qc, temp1, temp2, sqrt_term, NULL);
}


void lift_exact_charts(mpfr_point *p,  mpfr_t a,  mpfr_t b, int i) {
	mpfr_t temp0;
	mpfr_init(temp0);

	if(i == 1){
		lift_exact(p, a, b);
	}
	if(i == 2){
		lift_exact(p, a, b);
		mpfr_set(temp0, p->y, MPFR_RNDZ);
		mpfr_set(p->y, p->z, MPFR_RNDZ);
		mpfr_set(p->z, temp0, MPFR_RNDZ);
	}
	if(i == 3){
		lift_exact(p, a, b);
		mpfr_set(temp0, p->z, MPFR_RNDZ);
		mpfr_set(p->z, p->y, MPFR_RNDZ);
		mpfr_set(p->y, p->x, MPFR_RNDZ);
		mpfr_set(p->x, temp0, MPFR_RNDZ);
	}
	if(i == 4){
		lift_exact_under(p, a, b);
	}
	if(i == 5){
		lift_exact_under(p, a, b);
		mpfr_set(temp0, p->y, MPFR_RNDZ);
		mpfr_set(p->y, p->z, MPFR_RNDZ);
		mpfr_set(p->z, temp0, MPFR_RNDZ);
	}
	if(i == 6){
		lift_exact_under(p, a, b);
		mpfr_set(temp0, p->z, MPFR_RNDZ);
		mpfr_set(p->z, p->y, MPFR_RNDZ);
		mpfr_set(p->y, p->x, MPFR_RNDZ);
		mpfr_set(p->x, temp0, MPFR_RNDZ);
	}
	mpfr_clears(temp0, NULL);
}

void multiplyMatricesExact(mpfr_t **firstMatrix, mpfr_t **secondMatrix, mpfr_t **result, int ROWS1, int COLS1, int ROWS2, int COLS2) {
    for (int i = 0; i < ROWS1; ++i) {
        for (int j = 0; j < COLS2; ++j) {
            mpfr_set_ui(result[i][j], 0, MPFR_RNDZ);
            for (int k = 0; k < COLS1; ++k) {
                mpfr_t temp;
                mpfr_init(temp);
                mpfr_mul(temp, firstMatrix[i][k], secondMatrix[k][j], MPFR_RNDZ);
                mpfr_add(result[i][j], result[i][j], temp, MPFR_RNDZ);
                mpfr_clear(temp);
            }
        }
    }
}

void multiplyMatricesAuxExact(mpfr_t **firstMatrix, mpfr_t **secondMatrix, mpfr_t **result, int ROWS1, int COLS1, int ROWS2, int COLS2) {
    for (int i = 0; i < ROWS1; ++i) {
        for (int j = 0; j < COLS2; ++j) {
            mpfr_set_ui(result[i][j], 0, MPFR_RNDZ);
            for (int k = 0; k < COLS1; ++k) {
                mpfr_t temp;
                mpfr_init(temp);
                mpfr_mul(temp, firstMatrix[i][k], secondMatrix[k][j], MPFR_RNDZ);
                mpfr_add(result[i][j], result[i][j], temp, MPFR_RNDZ);
                mpfr_clear(temp);
            }
        }
    }
    for (int i = 0; i < ROWS1; ++i) {
        for (int j = 0; j < COLS2; ++j) {
            mpfr_set(secondMatrix[i][j], result[i][j], MPFR_RNDZ);
        }
    }
}

void printMatrixExact(mpfr_t **matrix, int ROWS1, int COLS1) {
    fprintf(stderr, "{");
    for (int i = 0; i < ROWS1; ++i) {
        fprintf(stderr, "{");
        for (int j = 0; j < COLS1; ++j) {
            char buffer[64];
            mpfr_sprintf(buffer, "%.15Rf", matrix[i][j]);
            fprintf(stderr, "%s", buffer);
            if (j < COLS1 - 1) {
                fprintf(stderr, ",");
            }
        }
        if (i == ROWS1 - 1) {
            fprintf(stderr, "}}\n");
        } else {
            fprintf(stderr, "},\n");
        }
    }
}
void alpha_xExact(mpfr_t result, const mpfr_t x, const mpfr_t y) {
    mpfr_t temp1, temp2, temp3;
    mpfr_inits(temp1, temp2, temp3, NULL);

    // y * (1 - x^2) / ((1 + x^2)^2 * (1+y^2))

    //mpfr_printf("x,y is %.10Rf, %.10Rf\n", x,y);

    // Calculate temp1 = y^2 + 1
    mpfr_mul(temp1, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp1, temp1, 1, MPFR_RNDZ);
   // mpfr_printf("temp1 (y^2 + 1): %.10Rf\n", temp1);

    // Calculate temp2 = (1 + x^2)^2
    mpfr_mul(temp2, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul(temp2, temp2, temp2, MPFR_RNDZ);
   // mpfr_printf("temp2 ((1 + x^2)^2): %.10Rf\n", temp2);

    // Calculate temp1 = (y^2 + 1) * (1 + x^2)^2
    mpfr_mul(temp1, temp1, temp2, MPFR_RNDZ);
   // mpfr_printf("temp1 ((y^2 + 1) * (1 + x^2)^2): %.10Rf\n", temp1);

    // Calculate temp2 = 1 - x^2
    mpfr_mul(temp2, x, x, MPFR_RNDZ);
    mpfr_ui_sub(temp2, 1, temp2, MPFR_RNDZ);
   // mpfr_printf("temp2 (1 - x^2): %.10Rf\n", temp2);

    // Calculate temp2 = (1 - x^2) * y
    mpfr_mul(temp2, temp2, y, MPFR_RNDZ);
   // mpfr_printf("temp2 ((1 - x^2) * y): %.10Rf\n", temp2);

    // Calculate result = temp2 / temp1
    mpfr_div(result, temp2, temp1, MPFR_RNDZ);
   // mpfr_printf("result: %.10Rf\n", result);

    // Clear MPFR variables
    mpfr_clears(temp1, temp2, temp3, NULL);
}
void alpha_yExact(mpfr_t result, const mpfr_t x, const mpfr_t y) {
    alpha_xExact(result, y, x);
}

void DixExact(mpfr_point p, mpfr_t **matrix) {
    mpfr_set_si(matrix[0][0], -1, MPFR_RNDZ);
    mpfr_t temp1, temp2, temp3;
    mpfr_inits(temp1, temp2, temp3, NULL);

    // matrix[0][1] = -a  * p.z * (1 - p.y * p.y) / ((1 + p.z * p.z) * (1 + p.y * p.y) * (1 + p.y * p.y))
    
    alpha_xExact(temp1, p.y, p.z);
    mpfr_mul_ui(matrix[0][1], temp1, 10, MPFR_RNDZ);
    mpfr_neg(matrix[0][1], matrix[0][1], MPFR_RNDZ);



    // matrix[0][2] = -a * p.y * (1 - p.z * p.z) / ((1 + p.y * p.y) * (1 + p.z * p.z) * (1 + p.z * p.z))
    alpha_yExact(temp2, p.y, p.z);
    mpfr_mul_ui(matrix[0][2], temp2, 10, MPFR_RNDZ);
    mpfr_neg(matrix[0][2], matrix[0][2], MPFR_RNDZ);


    mpfr_set_ui(matrix[1][0], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[1][1], 1, MPFR_RNDZ);
    mpfr_set_ui(matrix[1][2], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[2][0], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[2][1], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[2][2], 1, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, NULL);
}

void DiyExact(mpfr_point p, mpfr_t **matrix) {
    mpfr_set_ui(matrix[0][0], 1, MPFR_RNDZ);
    mpfr_set_ui(matrix[0][1], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[0][2], 0, MPFR_RNDZ);

    mpfr_t temp1, temp2, temp3;
    mpfr_inits(temp1, temp2, temp3, NULL);

    // matrix[1][0] = -a * p.z * (1 - p.x * p.x) / ((1 + p.z * p.z) * (1 + p.x * p.x) * (1 + p.x * p.x))
    alpha_xExact(temp1, p.x, p.z);
    mpfr_mul_ui(matrix[1][0], temp1, 10, MPFR_RNDZ);
    mpfr_neg(matrix[1][0], matrix[1][0], MPFR_RNDZ);

    mpfr_set_si(matrix[1][1], -1, MPFR_RNDZ);

    // matrix[1][2] = -a * p.x * (1 - p.z * p.z) / ((1 + p.x * p.x) * (1 + p.z * p.z) * (1 + p.z * p.z))
    alpha_yExact(temp2, p.x, p.z);
    mpfr_mul_ui(matrix[1][2], temp2, 10, MPFR_RNDZ);
    mpfr_neg(matrix[1][2], matrix[1][2], MPFR_RNDZ);


    mpfr_set_ui(matrix[2][0], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[2][1], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[2][2], 1, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, NULL);
}

void DizExact(mpfr_point p, mpfr_t **matrix) {
    mpfr_set_ui(matrix[0][0], 1, MPFR_RNDZ);
    mpfr_set_ui(matrix[0][1], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[0][2], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[1][0], 0, MPFR_RNDZ);
    mpfr_set_ui(matrix[1][1], 1, MPFR_RNDZ);
    mpfr_set_ui(matrix[1][2], 0, MPFR_RNDZ);
    mpfr_t temp1, temp2, temp3;
    mpfr_inits(temp1, temp2, temp3, NULL);

    // matrix[2][0] = -a * p.y * (1 - p.x * p.x) / ((1 + p.y * p.y) * (1 + p.x * p.x) * (1 + p.x * p.x))
    alpha_xExact(temp1, p.x, p.y);
    mpfr_mul_ui(matrix[2][0], temp1, 10, MPFR_RNDZ);
    mpfr_neg(matrix[2][0], matrix[2][0], MPFR_RNDZ);

    // matrix[2][1] = -a * p.x * (1 - p.y * p.y) / ((1 + p.x * p.x) * (1 + p.y * p.y) * (1 + p.y * p.y))
    alpha_yExact(temp1, p.x, p.y);
    mpfr_mul_ui(matrix[2][1], temp1, 10, MPFR_RNDZ);
    mpfr_neg(matrix[2][1], matrix[2][1], MPFR_RNDZ);

    mpfr_set_si(matrix[2][2], -1, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, NULL);
}

void DExact(mpfr_t result, mpfr_t x, mpfr_t y) {
    mpfr_t temp1, temp2, temp3, temp4;
    mpfr_inits(temp1, temp2, temp3, temp4, NULL);

    // temp1 = 100 * x * x * y * y
    mpfr_mul(temp1, x, x, MPFR_RNDZ);
    mpfr_mul(temp1, temp1, y, MPFR_RNDZ);
    mpfr_mul(temp1, temp1, y, MPFR_RNDZ);
    mpfr_mul_ui(temp1, temp1, 100, MPFR_RNDZ);

    // temp2 = 8 * (1 + x * x) * (1 + y * y)
    mpfr_mul(temp2, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul(temp3, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp3, temp3, 1, MPFR_RNDZ);
    mpfr_mul(temp2, temp2, temp3, MPFR_RNDZ);
    mpfr_mul_ui(temp2, temp2, 8, MPFR_RNDZ);

    // temp3 = 4 * (1 + x * x) * (1 + x * x) * (1 + y * y) * (1 + y * y)
    mpfr_mul(temp3, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp3, temp3, 1, MPFR_RNDZ);
    mpfr_mul(temp3, temp3, temp3, MPFR_RNDZ);
    mpfr_mul(temp4, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp4, temp4, 1, MPFR_RNDZ);
    mpfr_mul(temp4, temp4, temp4, MPFR_RNDZ);
    mpfr_mul(temp3, temp3, temp4, MPFR_RNDZ);
    mpfr_mul_ui(temp3, temp3, 4, MPFR_RNDZ);

    // result = temp1 + temp2 - temp3
    mpfr_add(result, temp1, temp2, MPFR_RNDZ);
    mpfr_sub(result, result, temp3, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, temp4, NULL);
}

void DxExact(mpfr_t result, mpfr_t x, mpfr_t y) {
    mpfr_t temp1, temp2, temp3, temp4;
    mpfr_inits(temp1, temp2, temp3, temp4, NULL);

    // temp1 = 200 * x * y * y
    mpfr_mul(temp1, x, y, MPFR_RNDZ);
    mpfr_mul(temp1, temp1, y, MPFR_RNDZ);
    mpfr_mul_ui(temp1, temp1, 200, MPFR_RNDZ);

    // temp2 = 16 * x * (1 + y * y)
    mpfr_mul(temp2, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul(temp2, temp2, x, MPFR_RNDZ);
    mpfr_mul_ui(temp2, temp2, 16, MPFR_RNDZ);

    // temp3 = 16 * x * (1 + x * x) * (1 + y * y) * (1 + y * y)
    mpfr_mul(temp3, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp3, temp3, 1, MPFR_RNDZ);
    mpfr_mul(temp4, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp4, temp4, 1, MPFR_RNDZ);
    mpfr_mul(temp4, temp4, temp4, MPFR_RNDZ);
    mpfr_mul(temp3, temp3, temp4, MPFR_RNDZ);
    mpfr_mul(temp3, temp3, x, MPFR_RNDZ);
    mpfr_mul_ui(temp3, temp3, 16, MPFR_RNDZ);

    // result = temp1 + temp2 - temp3
    mpfr_add(result, temp1, temp2, MPFR_RNDZ);
    mpfr_sub(result, result, temp3, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, temp4, NULL);
}
void DyExact(mpfr_t result, mpfr_t x, mpfr_t y) {
    DxExact(result, y, x);
}

void DxxExact(mpfr_t result, mpfr_t x, mpfr_t y) {
    mpfr_t temp1, temp2, temp3, temp4, temp5;
    mpfr_inits(temp1, temp2, temp3, temp4, temp5, NULL);

    // temp1 = 200 * y * y
    mpfr_mul(temp1, y, y, MPFR_RNDZ);
    mpfr_mul_ui(temp1, temp1, 200, MPFR_RNDZ);

    // temp2 = 16 * (1 + y * y)
    mpfr_mul(temp2, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul_ui(temp2, temp2, 16, MPFR_RNDZ);

    // temp3 = 32 * x * x * (1 + y * y) * (1 + y * y)
    mpfr_mul(temp3, x, x, MPFR_RNDZ);
    mpfr_mul(temp4, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp4, temp4, 1, MPFR_RNDZ);
    mpfr_mul(temp4, temp4, temp4, MPFR_RNDZ);
    mpfr_mul(temp3, temp3, temp4, MPFR_RNDZ);
    mpfr_mul_ui(temp3, temp3, 32, MPFR_RNDZ);

    // temp4 = 16 * (1 + x * x) * (1 + y * y) * (1 + y * y)
    mpfr_mul(temp4, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp4, temp4, 1, MPFR_RNDZ);
    mpfr_mul(temp5, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp5, temp5, 1, MPFR_RNDZ);
    mpfr_mul(temp5, temp5, temp5, MPFR_RNDZ);
    mpfr_mul(temp4, temp4, temp5, MPFR_RNDZ);
    mpfr_mul_ui(temp4, temp4, 16, MPFR_RNDZ);

    // result = temp1 + temp2 - temp3 - temp4
    mpfr_add(result, temp1, temp2, MPFR_RNDZ);
    mpfr_sub(result, result, temp3, MPFR_RNDZ);
    mpfr_sub(result, result, temp4, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, temp4, temp5, NULL);
}
void DyyExact(mpfr_t result, mpfr_t x, mpfr_t y) {
    DxxExact(result, y, x);
}
void DxyExact(mpfr_t result, mpfr_t x, mpfr_t y) {
    mpfr_t temp1, temp2, temp3;
    mpfr_inits(temp1, temp2, temp3, NULL);

    // temp1 = 432 * x * y
    mpfr_mul(temp1, x, y, MPFR_RNDZ);
    mpfr_mul_ui(temp1, temp1, 432, MPFR_RNDZ);

    // temp2 = 64 * x * (1 + x * x) * y * (1 + y * y)
    mpfr_mul(temp2, x, x, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);
    mpfr_mul(temp2, temp2, x, MPFR_RNDZ);
    mpfr_mul(temp2, temp2, y, MPFR_RNDZ);
    mpfr_mul(temp3, y, y, MPFR_RNDZ);
    mpfr_add_ui(temp3, temp3, 1, MPFR_RNDZ);
    mpfr_mul(temp2, temp2, temp3, MPFR_RNDZ);
    mpfr_mul_ui(temp2, temp2, 64, MPFR_RNDZ);

    // result = temp1 - temp2
    mpfr_sub(result, temp1, temp2, MPFR_RNDZ);

    mpfr_clears(temp1, temp2, temp3, NULL);
}

void pxExact(mpfr_t result, mpfr_t s, mpfr_t t) {
    mpfr_t numerator, denominator, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    mpfr_inits(numerator, denominator, temp1, temp2, temp3, temp4, temp5, temp6, temp7, NULL);

    // numerator = s * (-2 + 23 * t * t) - s * s * s * (2 + 27 * t * t)
    mpfr_mul(temp1, t, t, MPFR_RNDZ);
    mpfr_mul_ui(temp1, temp1, 23, MPFR_RNDZ);
    mpfr_add_si(temp1, temp1, -2, MPFR_RNDZ);
    mpfr_mul(temp1, temp1, s, MPFR_RNDZ);

    mpfr_mul(temp2, s, s, MPFR_RNDZ);
    mpfr_mul(temp2, temp2, s, MPFR_RNDZ);

    mpfr_mul_ui(temp3, t, 27, MPFR_RNDZ);
    mpfr_mul(temp3, t, temp3, MPFR_RNDZ);
    mpfr_add_ui(temp3, temp3, 2, MPFR_RNDZ);

    mpfr_mul(temp2, temp2, temp3, MPFR_RNDZ);
    mpfr_sub(temp1, temp1, temp2, MPFR_RNDZ);

    //temp4 = 1 - t^4 - s^4 * (1 + t^2)^2
    mpfr_pow_ui(temp4, t, 4, MPFR_RNDZ);
    mpfr_pow_ui(temp5, s, 4, MPFR_RNDZ);
    mpfr_mul(temp6, t, t, MPFR_RNDZ);
    mpfr_add_ui(temp6, temp6, 1, MPFR_RNDZ);
    mpfr_mul(temp6, temp6, temp6, MPFR_RNDZ);
    mpfr_mul(temp5, temp5, temp6, MPFR_RNDZ);
    mpfr_ui_sub(temp4, 1, temp4, MPFR_RNDZ); //ethan changed
    mpfr_sub(temp4, temp4, temp5, MPFR_RNDZ);

    //temp4 = 1 - t^4 - s^4 * (1 + t^2)^2 + s^2 * (23 * t^2 - 2 * t^4)
    mpfr_mul(temp5, s, s, MPFR_RNDZ);
    mpfr_mul_ui(temp6, t, 23, MPFR_RNDZ);
    mpfr_mul(temp6, temp6, t, MPFR_RNDZ);
    mpfr_pow_ui(temp7, t, 4, MPFR_RNDZ);
    mpfr_mul_ui(temp7, temp7, 2, MPFR_RNDZ);
    mpfr_sub(temp6, temp6, temp7, MPFR_RNDZ);
    mpfr_mul(temp6, temp6, s, MPFR_RNDZ);
    mpfr_mul(temp6, temp6, s, MPFR_RNDZ);
    mpfr_add(temp4, temp4, temp6, MPFR_RNDZ);

    //temp4 = sqrt(1 - t^4 - s^4 * (1 + t^2)^2 + s^2 * (23 * t^2 - 2 * t^4))
    mpfr_sqrt(temp4, temp4, MPFR_RNDZ);


    // numerator += -5 * t * temp4 + 5 * s^2 * t * temp4
    mpfr_mul(temp5, t, temp4, MPFR_RNDZ);
    mpfr_mul_si(temp5, temp5, -5, MPFR_RNDZ);
    mpfr_add(numerator, temp1, temp5, MPFR_RNDZ);

    mpfr_mul(temp5, s, s, MPFR_RNDZ);
    mpfr_mul(temp5, temp5, t, MPFR_RNDZ);
    mpfr_mul(temp5, temp5, temp4, MPFR_RNDZ);
    mpfr_mul_ui(temp5, temp5, 5, MPFR_RNDZ);
    mpfr_add(numerator, numerator, temp5, MPFR_RNDZ);

    // denominator = (1 + s^2)^2 * (1 + t^2) * sqrt(1 - t^4 - s^4 * (1 + t^2)^2 + s^2 * (23 * t^2 - 2 * t^4))
    mpfr_mul(temp1, s, s, MPFR_RNDZ);
    mpfr_add_ui(temp1, temp1, 1, MPFR_RNDZ);
    mpfr_mul(temp1, temp1, temp1, MPFR_RNDZ);

    mpfr_mul(temp2, t, t, MPFR_RNDZ);
    mpfr_add_ui(temp2, temp2, 1, MPFR_RNDZ);

    mpfr_mul(denominator, temp1, temp2, MPFR_RNDZ);
    mpfr_mul(denominator, denominator, temp4, MPFR_RNDZ);

    // result = numerator / denominator
    mpfr_div(result, numerator, denominator, MPFR_RNDZ);

    mpfr_clears(numerator, denominator, temp1, temp2, temp3, temp4, temp5, temp6, temp7, NULL);
}

void pyExact(mpfr_t result, mpfr_t s, mpfr_t t) {
    pxExact(result, t, s);
}