#define PI 3.1415926535897932384  
#define TWO30  1073741824

typedef struct {double x,y;} complex;

complex add(), sub(), mult(), divide(), recip(), cx_conj();
complex cx_sqrt(), contsqrt(), polar();
complex cx_exp(), cx_log(), cx_sin(), cx_cos(), cx_sinh(), cx_cosh(), power();
complex disk_to_sphere(), mobius();

double arg(), cx_abs(), infnorm();
