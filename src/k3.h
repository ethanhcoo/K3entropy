/* K3 SURFACE TRAJECTORIES */
#include <stdbool.h>

/* Point on surface */
typedef struct {double x, y, z;} point;

/* Function types */
point ix(point p), iy(point p), iz(point p), f(point p), lift(double x, double y), rotate(point p), sym(point p);
double dist(point p, point q), involute(double x, double y, double z), norm(point p), pathlength(int n), surf();
void init_color_map(), init_surface(), orbit(), random_orbits(), usage(), 
	scan_window(), scan_dim(), setk3(double a, double b), scan_end(), scan_gray(), 
	scan_set(), ps_colorimage(), printhex(), hlsrgb(), ps_image(),
	patherr(), suberr();

void ps_line(point p, point q, int unravel), di(point *ps, int n), setsep(double s), draw_path(int n), draw_manifold(), ps_close();

/*Ethan's additions start*/
point grad(point p);
bool is_recurrent(point p);
void recurrent_search();
void search_near(point p, double epsilon, int N);
int connect_path(point p, point w, int n);
void pprint(point p); /*prints point*/
void ps_line_color(point p, point q, int unravel); /*colored segment in ps*/
void plot_orbit(point marked_orbit[], point p, int n); /*plots orbit in stable ps*/
void ps_dot(point p, int unravel);/*dot in ps*/
void ps_dot_transparent(point p);
void draw_background(point z, int n);
void subdivide(point p,point q,point ps[],int n);
void subdivide_plane(point p,point q,point ps[],int n);
int compare(const void *a, const void *b);
void star_algorithm(int *arr1, int *arr2, int *size, int *twists, int *len, int k);
void determine_position(double *x_vals, double *y_vals, point p);
point rescale (point p); /*rescales a point to the unit ball. */
point stereo_proj(point p); /*projects points on sphere with respect to (-1,0,0)*/
void ps_open(int argc,char *argv[]);
void ps_window(double centerx, double centery, double radius);
int fpath(point ps[],int np);
void fully_simplify_arrays(int *arr1, int *arr2, int *size);
bool simplify_arrays(int *arr1, int *arr2, int *size);
void s(int *arr1, int *arr2, int *size, int i);
void S(int *arr1, int *arr2, int *size, int i);
bool simplifying_move(int *arr1, int *arr2, int *size, int *twists, int *len);
void determine_arrays(double *x_vals, double *y_vals, int *pos, int *height, int *m, point *marked_orbit, int *b);
void delete_indices(int *arr1, int *arr2, int *size, int *arr1_temp, int *arr2_temp, int k, int *len);
void process_array(int *arr, int *size, int k);
void simplify_final_array(int *arr, int *size);
point transform(point p);
void reverse_array(int *arr, int size);
bool bool_simplify_final_array(int *arr, int *size);

void Dix(point p, double **matrix);
void Diy(point p,double **matrix);
void Diz(point p, double **matrix);
point vtan(point p);
void multiplyMatrices(double **firstMatrix, double **secondMatrix, double **result, int ROWS1, int COLS1, int ROWS2, int COLS2);
void printMatrix(double **matrix, int ROWS1, int COLS1);
double mat_norm(double **matrix, int ROWS1, int COLS1);
void multiplyMatricesAux(double **firstMatrix, double **secondMatrix, double **result, int ROWS1, int COLS1, int ROWS2, int COLS2);



point inv_transform(point p);
point inv_proj(point p);
point newton_plane(point p);
void subdivide_plane(point p,point q,point ps[],int n);

point newton(point p);


/*Ethan's additions end*/


int top(point p), unreal(double x, double y);
int topx(point p);
int topy(point p);

/* K3 surface parameters */
#define A 10.0
#define B 2.0

/* Window */
#define CENTERX 0.0
#define CENTERY 0.0
#define CENTER {0.0,0.0}
#define IDIM 900
#define JDIM 900
#define RADIUS 1.2

/* Length of stable manifold to draw */
#define LENGTH 100

/* Iteration limit */
#define ITERLIM 10000

/* Epsilon for Newton convergence */
#define NEWTONEPS 1e-12
#define NEWTONMAX 10000		/* Max number of steps */ //ETHAN CHANGED FROM 100
#define NEWTONCAUTION 1.0	/* Multiple of grad to follow */

/* Subdivisions and point on stable manifold &  */
#define SUBMAX 100000 //Ethan changed from 1000
#define PSMAX 3500000 //Ethan changed from 3500000

/* Max separation of points in plot */
#define SEP 0.002 /*Ethan changed it from 0.02*/

/* Display extra info */
#define VERBOSE 1
