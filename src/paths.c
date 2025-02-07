
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



//very important for determine_arrays
int pos_array[PSMAX]; //HOW BIG SHOULD I MAKE THESE??
int above_array[PSMAX];

int orbit_length = 23;
const char *word[PSMAX];
int above = -100000;
int current_pos = -10000;

point initial_point()
{
	point p;
	p = lift(sep/10,-sep/10);
	return(p);
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

int connect_path(point z, point w, int k) /*connects z and w by a straightline and gradient-flows to surface. Iterates path by f^k times. Returns array length( not, MINUS ONE).*/
{
	int n;
	point z_plane = transform(z);
	point w_plane = transform(w);
	ps[0] = z;
	ps[1] = w;
	n = dist(ps[0], ps[1])/sep;
	ps[n] = ps[1]; 
	subdivide_plane(ps[0], ps[1], ps, n); //edited so that initial path is just a path in the plane
	for (int i = 0; i <k; ++i){
		n=fpath(ps,n);
	}
	//draw_path(n); 
	return n+1; //stuff Kurt coded runs on size - 1. this returns size.
}

void plot_orbit(point marked_orbit[], point z, int n)/*Ethan added this*/
{
	marked_orbit[0] = z;
	for (int i = 0; i<n; i++){
		marked_orbit[i] = z;
		ps_dot(z, unravel);
		z = f(z);
	}
}

void draw_background(point z, int n)/*Ethan added this*/
{
	for (int i = 0; i<n; i++){
		ps_dot_transparent(z);
		z = f(z);
	}
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

bool is_recurrent(point p) /*returns true of p is delta, k-recurrent for some k<=22*/
{
    double delta = .00002;
	point q = f(p);
    for (int k = 1; k <= 50; k++) {
		if(dist(p, q) < delta) {
			fprintf(stderr, "(%d)\n", k);
			return true;
    	}
		q = f(q);
    };
    return false;
}

void recurrent_search() /*searches in epsilon-steps for points that are delta, k-recurrent for some k<=22*/
{
	double epsilon = .001;
    double x = .1;
	double y = .1;    
    int length = radius/epsilon;
    point p = lift(x,y);
    for (int j = 1; j <= length; j++) {
		for(int i = 1; i <= length; i++){
		    if(!unreal(x,y)){
                point p = lift(x,y);
                if(is_recurrent(p)) {
                    fprintf(stderr, "Success at (%.5lf, %.5lf, %.5lf)\n", p.x, p.y, p.z);
                };
            };
        	x = x + epsilon;
		};
		y = y + epsilon;
		x = .1; //ethan changed
    };
}

void determine_arrays(double *x_vals, double *y_vals, int *pos, int *height, int *m, point *marked_orbit, int *b) //input pointers to arrays and common array length
{	
	//records position of path
	pos_array[0] = -1000;
	above_array[0] = -1000; //path_origin

	//fprintf(stderr, "Determine_arrays - initialized\n");

	for(int t = 1; t < *m; t++){ 
		determine_position(x_vals, y_vals, ps[t]);
		//fprintf(stderr, "first position determined\n");
		//fprintf(stderr, "value of b is %d\n", *b);
		if(*b==1){
			pos_array[*b] = current_pos;
			//fprintf(stderr, "%f\n", transform(ps[t]).x);
			above_array[*b] = 0; //this does not matter
			(*b)++;
			//fprintf(stderr, "case 1\n");
			continue;
		}
		//fprintf(stderr, "made it");
		if(current_pos - pos_array[(*b)-1] == -1){ //moving left
			//fprintf(stderr, "case 2\n");
			pos_array[*b] = current_pos;
			//fprintf(stderr, " %d\n", current_pos);
			if(transform(ps[t]).y > y_vals[current_pos]){
				above_array[*b] = 1;
				//fprintf(stderr, " not here\n");
			} else {
				above_array[*b] = 0;
			}
			//fprintf(stderr, "(%d, %d)\n", pos_array[b], above_array[b]);
			//determine_position(ps[t+200]);
			//fprintf(stderr, "HERE %d\n", current_pos); //added for testing
			//pprint(transform(ps[t+200]));
			//fprintf(stderr, "%f\n", transform(ps[t+50]).x); //something off here: (1,0) not being counted
			(*b)++;
			continue;
		}
		//fprintf(stderr, "made it");
		if(current_pos - pos_array[(*b)-1] == 1){ //moving right
			//fprintf(stderr, "case 3\n");
			pos_array[*b] = current_pos;
			//fprintf(stderr, "%d", current_pos);
			if(transform(ps[t]).y > y_vals[pos_array[(*b)-1]]){
				above_array[*b] = 1;
			} else {
				above_array[*b] = 0;
			}
			(*b)++;
			continue;
		}
		int last = pos_array[(*b)-1];
		if (pos_array[(*b)-1] < current_pos) { //moving right
			//fprintf(stderr, "%d", current_pos);
			for (int i = last + 1; i <= current_pos; i++) {
				pos_array[*b] = i;
				//fprintf(stderr, "%d", i);
				if((transform(ps[t-1]).y > y_vals[i]) && (transform(ps[t]).y > y_vals[i])){
					above_array[*b] = 1;
				} else if((transform(ps[t-1]).y < y_vals[i]) && (transform(ps[t]).y < y_vals[i])){
					above_array[*b] = 0;
				} else{
					fprintf(stderr, "issue with over/under\n");
					abort();
					//pprint(ps[t-1]);
					//pprint(ps[t]);
					//fprintf(stderr, "%f",  y_vals[i]); 
				}
				(*b)++;
			}
    	} else if (pos_array[(*b)-1] > current_pos) { //moving left
			for (int i = last - 1; i >= current_pos; i--) {
				pos_array[*b] = i;
				if((transform(ps[t-1]).y > y_vals[i]) && (transform(ps[t]).y > y_vals[i])){
					above_array[*b] = 1;
				} else if((transform(ps[t-1]).y < y_vals[i]) && (transform(ps[t]).y < y_vals[i])){
					above_array[*b] = 0;
				} else{
					fprintf(stderr, "issue with over/under");
					abort();
					//pprint(ps[t-1]);
					//pprint(ps[t]);
					//fprintf(stderr, "%f",  y_vals[i]); 
				}
				//fprintf(stderr, "(%d, %d)\n", pos_array[b], above_array[b]);
				(*b)++;
			}
		}
	}
	//fprintf(stderr, "Determine_arrays - first half\n");

	//length of pos_array = b

	//determine position of starting point. print error if fails
	int test = 0;
	for(int i  = 0; i < orbit_length; i++){
		if(dist(transform(ps[0]), transform(marked_orbit[i]))<1e-2){
			pos[0] = i;
			test = 1;
			//fprintf(stderr, "successful starting point assignment\n");
			break;
		}
	}
	if(!test){
		fprintf(stderr, "error in starting point assignment\n");
		pprint(ps[0]);
		abort();
	}
	// determine position of ending point. print error if fails
	test = 0;
	height[(*b)-1] = -1000;
	for(int i  = 0; i < orbit_length; i++){
		if(dist(transform(ps[(*m)-1]), transform(marked_orbit[i]))<1e-2){
			pos[(*b)-1] = i;
			test = 1;
			//fprintf(stderr, "successful endpoint assignment to %d\n", i);
			break;
		}
	}
	if(!test){
		fprintf(stderr, "error in endpoint assignment\n");
		pprint(ps[*m-1]);
		abort();
	}
	//determing position inbetween

	height[0] = -1000;
	for(int i = 1; i<(*b)-1; i++){
		if((pos_array[i] < pos_array[i+1])){
			pos[i] = pos_array[i];
			if(above_array[i+1]){
				height[i]=1;
			} else {
				height[i]=0;
			}
		}
		if((pos_array[i] > pos_array[i+1])){
			pos[i] = pos_array[i+1];
			if(above_array[i+1]){
				height[i]=1;
			} else {
				height[i]=0;
			}
		}
	}
	//length of pos = b
}

void determine_position(double *x_vals, double *y_vals, point p){
	above = 0;
	double x = transform(p).x;
	double y = transform(p).y;
	for (int l = 0; l < orbit_length-1; l++){
		if((x_vals[l] < x) && (x <= x_vals[l+1])) {
			current_pos = l+1;
			if(y > y_vals[l]){
				above = 1;
			}
		}
	}
	if(x > x_vals[orbit_length - 1]){
		current_pos = orbit_length;
		if(y > y_vals[orbit_length - 1]){
			above = 1;
		}
	}
	if(x <= x_vals[0]) {
		current_pos = 0;
		if(y > y_vals[0]){
			above = 1;
		}
	}
	}

int compare(const void *a, const void *b) 
{
    double diff = stereo_proj(rescale(*(point*)a)).x - stereo_proj(rescale(*(point*)b)).x;
    if (diff < 0) return -1;
    else if (diff > 0) return 1;
    else return 0;
}

bool simplify_arrays(int *arr1, int *arr2, int *size) {
    bool simplified = false;
	if ((arr1)[0] == (arr1)[1]){
			for (int j = 1; j < *size - 1; j++) {
                (arr1)[j] = (arr1)[j + 1];
                (arr2)[j] = (arr2)[j + 1];
            }
            (*size) -= 1;
            simplified = true;
			return simplified;
	}

	if ((arr1)[*size - 1] == (arr1)[*size - 2]){
			arr2[*size - 2] = -1000;
            (*size) -= 1;
            simplified = true;
			return simplified;
	}

    for (int i = 0; i < *size - 2; i++) {
        if ((arr1)[i] == (arr1)[i+1] && (arr2)[i] == (arr2)[i+1]) {
            //fprintf(stderr, "made it");
			// Remove the ith and (i+1)th elements
            for (int j = i; j < *size - 2; j++) {
                (arr1)[j] = (arr1)[j + 2];
                (arr2)[j] = (arr2)[j + 2];
            }
            *size -= 2;
            simplified = true;
            return simplified;;
        }
    }
    return simplified;
}

void fully_simplify_arrays(int *arr1, int *arr2, int *size) {
    bool simplified;
    do {
        simplified = simplify_arrays(arr1, arr2, size);
    } while (simplified && *size > 2);
}

//must apply s to SIMPLIFIED sequences
void s(int *arr1, int *arr2, int *size, int i) { //s_i makes COUNTERCLOCKWISE rotation
	 //j=1 since skipping start and <*size-2 since skipping end
	 //after changing array, is it okay to keep iterating?
	 //moving right:
	for (int j = 1; j < *size - 2; j++) { //*size - 1 increments after pointing, *size-- incremenents before pointing
		//fprintf(stderr, "made it");
        if (arr1[j] == i && arr2[j] == 1 && arr1[j + 1] == i + 1 && arr2[j + 1] == 0) {
            arr2[j] = 0;
            arr2[j + 1] = 1;
			continue;
        }
		if (arr1[j] == i && arr2[j] == 0 && arr1[j + 1] == i + 1 && arr2[j + 1] == 1) {
            // Make room for four new entries
            for (int k = *size + 3; k >= j + 5; k--) {
                arr1[k] = arr1[k - 4];
                arr2[k] = arr2[k - 4];
            }
            // Insert the new sequence
            arr1[j + 1] = i + 1; arr2[j + 1] = 0;
            arr1[j + 2] = i + 1; arr2[j + 2] = 1;
            arr1[j + 3] = i; arr2[j + 3] = 0;
            arr1[j + 4] = i; arr2[j + 4] = 1;
			j+=4;
			(*size) +=4;
			continue;
		}

		//moving left:
		if (arr1[j+1] == i && arr2[j+1] == 1 && arr1[j] == i + 1 && arr2[j] == 0) {
            arr2[j+1] = 0;
            arr2[j] = 1;
			continue;
        }
		if (arr1[j+1] == i && arr2[j+1] == 0 && arr1[j] == i + 1 && arr2[j] == 1) {
            // Make room for four new entries
            for (int k = *size + 3; k >= j + 5; k--) {
                arr1[k] = arr1[k - 4];
                arr2[k] = arr2[k - 4];
            }
            // Insert the new sequence
            arr1[j + 4] = i + 1; arr2[j + 4] = 0;
            arr1[j + 3] = i + 1; arr2[j + 3] = 1;
            arr1[j + 2] = i; arr2[j + 2] = 0;
            arr1[j + 1] = i; arr2[j + 1] = 1;
			j+=4;
			(*size) +=4;
			continue;
		}
		//turning around
		if(arr1[j] == i && arr1[j+1] == i && arr1[j+2] == i - 1){
			int first = arr2[j];		
			int second = arr2[j+1];
			for (int k = *size + 1; k >= j + 4; k--) {
                arr1[k] = arr1[k - 2];
                arr2[k] = arr2[k - 2];
            }
            // Insert the new sequence
            arr1[j] = i ; arr2[j] = 0;
			arr1[j+1] = i+1 ; arr2[j+1] = first;
			arr1[j+2] = i+1 ; arr2[j+2] = second;
            arr1[j+3] = i ; arr2[j + 3] = 0;
			j+=4;
			(*size) +=2;
			continue;
		}
		//turning atound
		if(arr1[j] == i+1 && arr1[j+1] == i+1  && arr1[j+2] == i + 2){
			int first = arr2[j];		
			int second = arr2[j+1];
			for (int k = *size + 1; k >= j + 4; k--) {
                arr1[k] = arr1[k - 2];
                arr2[k] = arr2[k - 2];
            }
            // Insert the new sequence
            arr1[j] = i+1 ; arr2[j] = 1;
			arr1[j+1] = i ; arr2[j+1] = first;
			arr1[j+2] = i ; arr2[j+2] = second;
            arr1[j+3] = i+1 ; arr2[j + 3] = 1;
			j+=4;
			(*size) +=2;
			continue;
		}
    } 

	//specific case for size = 2:
	if(arr1[0]==i && arr1[1]==i+1 && *size == 2){
		arr1[0] = i+1;
		arr1[1] = i;
	} else if(arr1[0]==i+1 && arr1[1]==i && *size == 2){
		arr1[0] = i;
		arr1[1] = i+1;
	}

	//change end of path before first move becuase changing size doesn't affect first move.
	//works the same as changing start of path
	if(arr1[*size-1]==i && arr1[*size-2]==i+1 && arr2[*size-2]== 0 ){
		arr2[*size-2] = -1000;
    	(*size)--;
	 } else if(arr1[*size-1]==i && arr1[*size-2]==i+1 && arr2[*size-2]== 1){//ethan changed this
    	arr1[*size+1] = i+1; arr2[*size+1] = -1000;
		arr1[*size] = i; arr2[*size] = 0;
		arr1[*size-1] = i; arr2[*size-1] = 1;
    	(*size) += 2;
	 }else if(arr1[*size-1]==i && arr1[*size-2]==i-1){ //SURE THIS WORKS???
		/*for (int i = *size; i > 0; i--) {
			arr1[i] = arr1[i - 1];
			arr2[i] = arr2[i - 1];
		}*/
    	arr1[*size] = i+1; arr2[*size] = -1000;
		arr1[*size-1] = i; arr2[*size-1] = 0;
    	(*size) += 1;
	 } 
	//accounting for ith spot = i+1
	 else if (arr1[*size-2]==i && arr1[*size-1]==i+1 && arr2[*size-2]== 1 ){ 	
		arr2[*size-2] = -1000;
    	(*size)--;
	 } else if(arr1[*size-2]==i && arr1[*size-1]==i+1 && arr2[1]== 0 ){
    	arr1[*size+1] = i; arr2[*size+1] = -1000;
		arr1[*size] = i+1; arr2[*size] = 1;
		arr1[*size-1] = i+1; arr2[*size-1] = 0;
    	(*size) += 2;
	 } 
	 else if(arr1[*size-2]==i+2 && arr1[*size-1]==i+1){
    	arr1[*size] = i; arr2[*size] = -1000;
		arr1[*size-1] = i+1; arr2[*size-1] = 1;
    	(*size) += 1;
	 }

	//first moves are different: need to test these more!
	//accounting for ith spot = start
	 if(arr1[0]==i && arr1[1]==i+1 && arr2[1]== 0 ){
		for (int i = 0; i < *size - 1; i++) {
        arr1[i] = arr1[i + 1];
		arr2[i] = arr2[i + 1];
    	}
		arr2[0] = -1000;
    	(*size)--;
	 } else if(arr1[0]==i && arr1[1]==i+1 && arr2[1]== 1 ){
		for (int i = *size + 1; i > 0; i--) {
			arr1[i] = arr1[i - 2];
			arr2[i] = arr2[i - 2];
		}
    	// Add the new entry to the beginning
    	arr1[0] = i+1; arr2[0] = -1000;
		arr1[1] = i; arr2[1] = 0;
		arr1[2] = i; arr2[2] = 1;
    	(*size) += 2;
	 }else if(arr1[0]==i && arr1[1]==i-1){
		for (int i = *size; i > 0; i--) {
			arr1[i] = arr1[i - 1];
			arr2[i] = arr2[i - 1];
		}
    	// Add the new entry to the beginning
    	arr1[0] = i+1; arr2[0] = -1000;
		arr1[1] = i; arr2[1] = 0;
    	(*size) += 1;
	 }

	 //accounting for ith spot = i+1
	 else if(arr1[1]==i && arr1[0]==i+1 && arr2[1]== 1 ){
		for (int i = 0; i < *size - 1; i++) {
        	arr1[i] = arr1[i + 1];
			arr2[i] = arr2[i + 1];
    	}
		arr2[0] = -1000;
    	(*size)--;
	 } else if(arr1[1]==i && arr1[0]==i+1 && arr2[1]== 0 ){
		for (int i = *size + 1; i > 0; i--) {
			arr1[i] = arr1[i - 2];
			arr2[i] = arr2[i - 2];
		}
    	// Add the new entry to the beginning
    	arr1[0] = i; arr2[0] = -1000;
		arr1[1] = i+1; arr2[1] = 1;
		arr1[2] = i+1; arr2[2] = 0;
    	(*size) += 2;
	 } 
	 else if(arr1[1]==i+2 && arr1[0]==i+1){
		for (int i = *size + 1; i > 0; i--) {
			arr1[i] = arr1[i - 1];
			arr2[i] = arr2[i - 1];
		}
    	// Add the new entry to the beginning
    	arr1[0] = i; arr2[0] = -1000;
		arr1[1] = i+1; arr2[1] = 1;
    	(*size) += 1;
	 }

	 //fully simplify result

	fully_simplify_arrays(arr1, arr2, size);
}

void S(int *arr1, int *arr2, int *size, int i){
	//flip
	for(int i = 1; i < *size-1; i++){
		if(arr2[i]==0){
			arr2[i]=1;

		} else if(arr2[i]==1){
			arr2[i]=0;
		}
	}
	
	//apply s
	s(arr1, arr2, size, i);

	//flip
	for(int i = 1; i < *size-1; i++){
		if(arr2[i]==0){
			arr2[i]=1;
		} else if(arr2[i]==1){
			//fprintf(stderr, "flipping");
			arr2[i]=0;
		}
	}
}

bool simplifying_move(int *arr1, int *arr2, int *size, int *twistsy, int *len) {
    bool simplified = false;
	//fprintf(stderr, "inside\n");
	//fprintf(stderr, "size is %d\n", *size);
	//array should never have size 0 or 1
	if(*size < 2){
		fprintf(stderr, "ERROR: ARRAY TOO SMALL");
	}
	//return false if array is fully simplified
	if(*size == 2){
		return simplified;
	}

	//simplify otherwise
	if(arr1[0] < arr1[1] && arr2[1] == 0){
		//apply s_{arr1[0]}, append s_{arr1[0]} to twist array
		twistsy[*len] = arr1[0];
		//fprintf(stderr, "in deep %d \n", twists[*len]);
		s(arr1, arr2, size, arr1[0]);
		(*len)++;
		simplified = true;
		return simplified;

	}
	if(arr1[0] < arr1[1] && arr2[1] == 1){
		//apply + append S_{arr1[0]}
		twistsy[*len] = -arr1[0]-1;
		//fprintf(stderr, "in deep %d \n", twists[*len]);
		S(arr1, arr2, size, arr1[0]); //STILL MUST DEFINE S
		(*len)++;
		simplified = true;
		return simplified;
	}
	if(arr1[0] > arr1[1] && arr2[1] == 1){
		//apply + append s_{arr1[1]}
		twistsy[*len] = arr1[1];
		//fprintf(stderr, "in deep %d \n", twists[*len]);
		s(arr1, arr2, size, arr1[1]);
		(*len)++;
		simplified = true;
		return simplified;
	}
	if(arr1[0] > arr1[1] && arr2[1] == 0){
		//apply + append S_{arr1[1]}
		twistsy[*len] = -arr1[1]-1;
		//fprintf(stderr, "in deep %d \n", twists[*len]);
		S(arr1, arr2, size, arr1[1]);
		(*len)++;
		simplified = true;
		return simplified;
	}
    fprintf(stderr, "ERROR: simplification failed");
	return simplified;
}

void star_algorithm(int *arr1, int *arr2, int *size, int *twistsy, int *len, int k){//simplified until path starts at k and moves one step in postive direction.
	// set length of twists to zero
	//fprintf(stderr, "inside star!");
	//*len = 0;
	//int indicator = -5000; //1 = positive twist, 0 = negative twist
	//int ind = -5000; // keeps track of location
	for (int j = 0; j < *size; j++){
			fprintf(stderr, "(%d, %d) \n", arr1[j], arr2[j]);
	}
	bool simplified;
    do {
		fprintf(stderr, "size is %d\n", *size);
        simplified = simplifying_move(arr1, arr2, size, twistsy, len);
			//print example path 
		//for (int j = 0; j < *size; j++){
		//	fprintf(stderr, "(%d, %d) \n", arr1[j], arr2[j]);
		//}
		
    } while (simplified && *size > 2);
	//fprintf(stderr, "inside star - simplified");
	// move to integer
	while(arr1[0] != k || arr1[*size - 1] != k+1){
		//for(int i = 0; i < *size ;i++){
		//	fprintf(stderr, "(%d, %d)\n", arr1[i], arr2[i]);
		//}
		if(arr1[0]>k){
			twistsy[*len] = arr1[0]-1;
			s(arr1, arr2, size, arr1[0]-1);
			(*len)++;
			//fprintf(stderr, "(%d, %d)\n", arr1[0], arr1[*size-1]);
			
			continue;
		}
		if(arr1[*size-1]< k+1){
			twistsy[*len] = arr1[*size] - 1;
			s(arr1, arr2, size, arr1[*size-1]);
			(*len)++;
			//fprintf(stderr, "(%d, %d) \n", arr1[0], arr1[*size-1]);
			continue;
		} 
		if(arr1[*size-1]> k+1){
			twistsy[*len] = arr1[*size-1]-1;
			s(arr1, arr2, size, arr1[*size-1]-1);
			(*len)++;
			//fprintf(stderr, "(%d, %d) \n", arr1[0], arr1[*size-1]);
			continue;
		}
		if(arr1[0]< k){
			twistsy[*len] = arr1[0];
			s(arr1, arr2, size, arr1[0]);
			(*len)++;
			//fprintf(stderr, "(%d, %d) \n", arr1[0], arr1[*size-1]);
			continue;
		}
	}

}

void delete_indices(int *arr1, int *arr2, int *size, int *arr1_temp, int *arr2_temp, int k, int *len) {
	int write_index = 0;
	//important that path starts from right side of contracted path
    for (int i = 0; i < *size; i++) {
        if (arr1[i] >= k) { //strict??
            arr1_temp[write_index] = arr1[i];
            arr2_temp[write_index] = arr2[i];
            write_index++;
        }
    }
    *len = write_index;
}



void process_array(int *arr, int *size, int k) {
	//fprintf(stderr, "inside process");
	int location = k; //location of mega-strand
	for(int i = 0; i < *size; i++) {
		//fprintf(stderr, "location is %d\n", location);
		//fprintf(stderr, "move is %d\n", arr[i]);
		if(arr[i]>=0){
			if(arr[i] == location){
				for(int j = *size - 1 + k; j > i + k ; j--){
					arr[j] = arr[j - k];
				}
				for(int j = i+1; j <= i + k; j++){
					arr[j] = location-(j-i);
				}
				//for(int j = i; j <= i + k; j++){
				//	fprintf(stderr, "appended %d\n", arr[j]);
				//}
				(*size) += k;
				i+=k;
				location++;
				continue;
			} else if(arr[i] == location-1){	
				for(int j = *size - 1 + k; j >= i + k ; j--){
					arr[j] = arr[j - k];
				}
				for(int j = 0; j < k; j++){
					arr[i + j] = location-1-k + j;
				}
				//for(int j = i; j <= i + k; j++){
				//	fprintf(stderr, "appended %d\n", arr[j]);
				//}
				(*size) += k;
				i+=k;
				location--;
				continue;
			}
		} else if(arr[i]<0){
			if(-arr[i]-1 == location){
				for(int j = *size - 1 + k; j > i + k ; j--){
					arr[j] = arr[j - k];
				}
				for(int j = i+1; j <= i + k; j++){
					arr[j] = -location-1+(j-i); //is this right?
				}
				//for(int j = i; j <= i + k; j++){
				//	fprintf(stderr, "appended %d\n", arr[j]);
				//}
				i+=k;
				(*size) += k;
				location++;
				continue;
			} else if(-arr[i]-1 == location-1){	
				for(int j = *size - 1 + k; j >= i + k ; j--){
					arr[j] = arr[j - k];
				}
				for(int j = i; j < i + k; j++){
					arr[j] = -location +k - (j-i); //is this right?
				}
				//for(int j = i; j <= i + k; j++){
				//	fprintf(stderr, "appended %d\n", arr[j]);
				//}
				(*size) += k;
				i+=k;
				location--;
				continue;
			}
		}
		//make space
	}
}

void simplify_final_array(int *arr, int *size){
    bool simplified;
    do {
        simplified = bool_simplify_final_array(arr, size);
    } while (simplified);
}

bool bool_simplify_final_array(int *arr, int *size){
	bool simplified = false;
	for (int i = 0; i < *size-1 ; i++) {
		// If consecutive terms add up to -1, skip them
		if (arr[i] + arr[i + 1] == -1) {
			//delete array elements
			for(int j = i; j < *size - 2; j++){
				arr[j] = arr[j + 2];
			}
			*size = (*size) -2;
			simplified = true;
			return simplified;
		}
	}
	return simplified;
}


void reverse_array(int *arr, int size){
	for (int i = 0; i < size / 2; i++) {
        int temp = arr[i];
        arr[i] = arr[size - 1 - i];
        arr[size - 1 - i] = temp;
	}
}
