
/* Calculates the derivative of f along an orbit */

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

//very important for determine_arrays
int pos_array[PSMAX]; //HOW BIG SHOULD I MAKE THESE??
int above_array[PSMAX];

int orbit_length = 23;
const char *word[PSMAX];
int above = -100000;
int current_pos = -10000;

int alpha = 10;
int beta = 10;

double px(double s, double t);
double py(double s, double t);
double pxx(double s, double t);
double pxy(double s, double t);
double pyy(double s, double t);
int chart(point p); //outputs 1,2,3,4,5,6 depending on chart p belongs to (might belong to two)

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
	
	int iters = 11;
	point p = lift(0.89460, 0.00480);
    point q = vtan(p);
    // Allocate memory for the array of pointers
    double** tangent = (double**)malloc(3 * sizeof(double*));
    double** tangent2 = (double**)malloc(3 * sizeof(double*));

    // Allocate memory for each row
    for (int i = 0; i < 3; ++i) {
        tangent[i] = (double*)malloc(3 * sizeof(double));
        tangent2[i] = (double*)malloc(3 * sizeof(double));
    }
    tangent[0][0] = q.x;
    tangent[1][0] = q.y;
    tangent[2][0] = q.z;

    tangent2[0][0] = 0;
    tangent2[1][0] = 0;
    tangent2[2][0] = 0;

   // fprintf(stderr, "MADE IT");

    double** matrix = (double**)malloc(3 * sizeof(double*));
    double** matrix_temp = (double**)malloc(3 * sizeof(double*));
    double** temp1 = (double**)malloc(3 * sizeof(double*));
    double** temp2 = (double**)malloc(3 * sizeof(double*));
    double** temp3 = (double**)malloc(3 * sizeof(double*));
    double** temp4 = (double**)malloc(3 * sizeof(double*));
    double** temp5 = (double**)malloc(3 * sizeof(double*));
    double** x_mat = (double**)malloc(3 * sizeof(double*));
    double** y_mat = (double**)malloc(3 * sizeof(double*));
    double** z_mat = (double**)malloc(3 * sizeof(double*));

    for (int i = 0; i < 3; ++i) {
        matrix[i] = (double*)malloc(3 * sizeof(double));
        matrix_temp[i] = (double*)malloc(3 * sizeof(double));
        temp1[i] = (double*)malloc(3 * sizeof(double));
        temp2[i] = (double*)malloc(3 * sizeof(double));
        temp3[i] = (double*)malloc(3 * sizeof(double));
        temp4[i] = (double*)malloc(3 * sizeof(double));
        temp5[i] = (double*)malloc(3 * sizeof(double));
        x_mat[i] = (double*)malloc(3 * sizeof(double));
        y_mat[i] = (double*)malloc(3 * sizeof(double));
        z_mat[i] = (double*)malloc(3 * sizeof(double));
    }
    for(int i = 0; i < 3; i++){
       // fprintf(stderr, "inside");
        for(int j = 0; j < 3; j++){
           // fprintf(stderr, "inside2");
            if(i==j){
                matrix_temp[i][j] = 1;
                //fprintf(stderr, "MADE IT");
                matrix[i][j] = 1;
                temp1[i][j] = 1;
                temp2[i][j] = 1;
                temp3[i][j] = 1;
                temp4[i][j] = 1;
                temp5[i][j] = 1;
                x_mat[i][j] = 1;
                y_mat[i][j] = 1;
                z_mat[i][j] = 1;

            } else{
                matrix_temp[i][j] = 0;
                matrix[i][j] = 0;
                temp1[i][j] = 0;
                temp2[i][j] = 0;
                temp3[i][j] = 0;
                temp4[i][j] = 0;
                temp5[i][j] = 0;
                x_mat[i][j] = 0;
                y_mat[i][j] = 0;
                z_mat[i][j] = 0;        
            }
        }
    }

    
    for(int i = 0; i < 22; i++){
        Dix(p, x_mat); 
        //fprintf(stderr, "%f", x_mat[0][0]);
        p = ix(p);
        Diy(p, y_mat); 
        p = iy(p);
        Diz(p, z_mat);
        p = iz(p); 
        multiplyMatrices(y_mat, x_mat, temp2, 3, 3, 3, 3);
        multiplyMatrices(z_mat, temp2, temp3, 3, 3, 3, 3);
        //pprint(p);
        //fprintf(stderr, "chart is %d\n", chart(p));
        //fprintf(stderr, "chart is %d\n", chart(f(p)));

        //fprintf(stderr, "matrix is \n");
        //printMatrix(temp3, 3, 3);
        int k = 0;
        for(int i = 0; i < 3; i++){
            if(i == 3-chart(p)){
                temp4[i][0] = px(p.x, p.y);
                temp4[i][1] = py(p.x, p.y);
            } else if(i == 6-chart(p)){
                temp4[i][0] = -px(-p.x, p.y);
                temp4[i][1] = -py(-p.x, p.y);
            } else{
                if(k == 0){
                    temp4[i][0] = 1;
                    temp4[i][1] = 0;
                }
                if(k == 1){
                    temp4[i][0] = 0;
                    temp4[i][1] = 1;
                }
                if(k == 2){
                    fprintf(stderr, "oops");
                }
                k++;
            }
            
        }
        k = 0;
        for(int j = 0; j < 3; j++){
            if(j == 3-chart(f(p))){
                temp5[0][j] = 0;
                temp5[1][j] = 0;
            } else if(j == 6-chart(f(p))){
                temp5[0][j] = 0;
                temp5[1][j] = 0;
            } else{
                if(k == 0){
                    temp5[0][j] = 1;
                    temp5[1][j] = 0;
                }
                if(k == 1){
                    temp5[0][j] = 0;
                    temp5[1][j] = 1;
                }
                if(k == 2){
                    fprintf(stderr, "oops");
                }
                k++;
            }
            
        }
        //fprintf(stderr, "trans 1 is \n");
        //printMatrix(temp4, 3, 2);
        //fprintf(stderr, "trans 2 is \n");
        //printMatrix(temp5, 2, 3);
        

        multiplyMatricesAux(temp5, temp3, temp1, 2, 3, 3, 3);
        multiplyMatricesAux(temp3, temp4, temp1, 2, 3, 3, 2);
        fprintf(stderr, "final matrix is \n");
        printMatrix(temp4, 2, 2);


        multiplyMatricesAux(x_mat, matrix, temp1, 3, 3, 3, 3);
        multiplyMatricesAux(y_mat, matrix, temp1, 3, 3, 3, 3);
        multiplyMatricesAux(z_mat, matrix, temp1, 3, 3, 3, 3);

        //printMatrix(matrix, 3, 3);
        //fprintf(sterr, "deep");
        if(i % 22 == 21){
            //fprintf(stderr, "argh");
            //multiplyMatricesAux(matrix, tangent, tangent2, 3, 3, 3, 1);
            //printMatrix(tangent, 3, 1);
            //fprintf(stderr, "%f \n", log(mat_norm(matrix, 3, 1))/(i/11));
            //pprint(p);
        }
    }
    //fprintf(stderr, "MADE IT");
    //printMatrix(tangent, 3, 1);
    //multiplyMatricesAux(matrix, tangent, tangent2, 3, 3, 3, 1);
    //printMatrix(tangent, 3, 1);
   //fprintf(stderr, "matrix is \n");
    //printMatrix(matrix, 3, 3);


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


double px(double s, double t) {
    double numerator = s * (-2 + 23 * t * t) - s * s * s * (2 + 27 * t * t) - 
                       5 * t * sqrt(1 - t * t * t * t - s * s * s * s * (1 + t * t) * (1 + t * t) + s * s * (23 * t * t - 2 * t * t * t * t)) + 
                       5 * s * s * t * sqrt(1 - t * t * t * t - s * s * s * s * (1 + t * t) * (1 + t * t) + s * s * (23 * t * t - 2 * t * t * t * t));
    
    double denominator = (1 + s * s) * (1 + s * s) * (1 + t * t) * 
                         sqrt(1 - t * t * t * t - s * s * s * s * (1 + t * t) * (1 + t * t) + s * s * (23 * t * t - 2 * t * t * t * t));
    
    return numerator / denominator;
}

double py(double s, double t) {
    return px(t,s);
}

double pxx(double s, double t) {
    double term1 = pow((200 * s * t * t + 16 * s * (1 + t * t) - 16 * s * (1 + s * s) * pow(1 + t * t, 2)), 2) / 
                   (4 * pow((100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2)), 1.5));
    
    double term2 = (200 * t * t + 16 * (1 + t * t) - 32 * s * s * pow((1 + t * t), 2) - 16 * (1 + s * s) * pow((1 + t * t), 2)) / 
                   (2 * sqrt(100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2)));
    
    double term3 = 2 * s * (-10 * t + (200 * s * t * t + 16 * s * (1 + t * t) - 16 * s * (1 + s * s) * pow((1 + t * t), 2)) / 
                            (2 * sqrt(100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2))));
    
    double term4 = 4 * s * s * (-10 * s * t + sqrt(100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2)));
    
    double term5 = -(-10 * s * t + sqrt(100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2)));
    
    double denominator = 2 * (1 + s * s) * (1 + t * t);
    
    return (-(term1) + term2) / denominator - term3 / pow((1 + s * s), 2) / (1 + t * t) + term4 / pow((1 + s * s), 3) / (1 + t * t) - term5 / pow((1 + s * s), 2) / (1 + t * t);
}

double pyy(double s, double t) {
    return pxx(t,s);
}

double pxy(double s, double t) {
    double term1 = (200 * s * s * t + 16 * (1 + s * s) * t - 16 * pow((1 + s * s), 2) * t * (1 + t * t)) * 
                   (200 * s * t * t + 16 * s * (1 + t * t) - 16 * s * pow((1 + s * s), 2) * pow((1 + t * t), 2));
    
    double term2 = 4 * pow((100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2)), 1.5);
    
    double term3 = 432 * s * t - 64 * s * (1 + s * s) * t * (1 + t * t);
    
    double term4 = 2 * sqrt(100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2));
    
    double term5 = s * (-10 * s + (200 * s * s * t + 16 * (1 + s * s) * t - 16 * pow((1 + s * s), 2) * t * (1 + t * t)) / term4);
    
    double term6 = t * (-10 * t + (200 * s * t * t + 16 * s * (1 + t * t) - 16 * s * pow((1 + s * s), 2) * pow((1 + t * t), 2)) / term4);
    
    double term7 = 2 * s * t * (-10 * s * t + sqrt(100 * s * s * t * t + 8 * (1 + s * s) * (1 + t * t) - 4 * pow((1 + s * s), 2) * pow((1 + t * t), 2)));
    
    double denominator = 2 * (1 + s * s) * (1 + t * t);
    
    return (-10 - term1 / term2 + term3 / term4) / denominator - term5 / pow((1 + s * s), 2) / (1 + t * t) - term6 / (1 + s * s) / pow((1 + t * t), 2) + term7 / pow((1 + s * s), 2) / pow((1 + t * t), 2);
}

int chart(point p){
    if(top(p)){
        if(fabs(px(p.x, p.y)) < alpha && fabs(py(p.x, p.y)) < alpha && fabs(pxx(p.x, p.y)) < beta && fabs(pxy(p.x, p.y)) < beta && fabs(pyy(p.x, p.y)) < beta) {
            return 1;
        }
        if(fabs(px(p.x, p.z)) < alpha && fabs(py(p.x, p.z)) < alpha && fabs(pxx(p.x, p.z)) < beta && fabs(pxy(p.x, p.z)) < beta && fabs(pyy(p.x, p.z)) < beta) {
            return 2;
        }
        if(fabs(px(p.y, p.z)) < alpha && fabs(py(p.y, p.z)) < alpha && fabs(pxx(p.y, p.z)) < beta && fabs(pxy(p.y, p.z)) < beta && fabs(pyy(p.y, p.z)) < beta) {
            return 3;
        }
    } else{
        //fprintf(stderr, "made it");
        if(fabs(px(-p.x, p.y)) < alpha && fabs(py(-p.x, p.y)) < alpha && fabs(pxx(-p.x, p.y)) < beta && fabs(-pxy(p.x, p.y)) < beta && fabs(-pyy(p.x, p.y)) < beta) {
            return 4;
        }
        if(fabs(px(-p.x, p.z)) < alpha && fabs(py(-p.x, p.z)) < alpha && fabs(pxx(-p.x, p.z)) < beta && fabs(-pxy(p.x, p.z)) < beta && fabs(-pyy(p.x, p.z)) < beta) {
            return 5;
        }
        if(fabs(px(-p.y, p.z)) < alpha && fabs(py(-p.y, p.z)) < alpha && fabs(pxx(-p.y, p.z)) < beta && fabs(-pxy(p.y, p.z)) < beta && fabs(-pyy(p.y, p.z)) < beta) {
            return 6;
        }
    }
    return 0;
}



//notice: p_+(-x, y) = -p_-(x,y)