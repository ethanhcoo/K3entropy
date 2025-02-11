#include "k3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


static double centerx = CENTERX, centery = CENTERY, radius = RADIUS;
static double length=LENGTH, sep = SEP;
static int drawtop = 1, drawbot = 1, verbose=VERBOSE, spin=0, unravel = 0;

//Very large array to dynamically store paths
point ps[PSMAX];

//very important for determine_arrays
int pos_array[PSMAX]; //HOW BIG SHOULD I MAKE THESE??
int above_array[PSMAX];

int orbit_length = 22;
const char *word[PSMAX];
int above = -100000;
int current_pos = -10000;



int main(int argc, char *argv[]) {
	int i,n;
	double a,b;

	//read in values from mclass.run
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


	//Creates arrays to record twists and keep track of its length.
	char *twists_word[10000]; //twists as characters
	int twists[10000]; //twist with integer coding. i = s_i. -i = S_{i-1}
	int twists_length = 0; //length of twist array

	//Create (orbit_length-1) pos/height arrays to record paths 
	int **array_pos = malloc((orbit_length-1) * sizeof(int*));
	int **array_height = malloc((orbit_length -1) * sizeof(int*));
	int **array_pos_testing = malloc((orbit_length-1) * sizeof(int*));
	int **array_height_testing = malloc((orbit_length -1) * sizeof(int*));
	int path_lengths[orbit_length-1];
	int path_lengths_testing[orbit_length-1];
	
	//Allocate memory for each array. second copy is for testing.
	for (int i = 0; i < orbit_length-1; i++) {
        array_pos[i] = malloc(100000 * sizeof(int)); //100,000 might not be enough - make larger if causing problems!!!!
		array_height[i] = malloc(100000 * sizeof(int));
		array_pos_testing[i] = malloc(100000 * sizeof(int)); 
		array_height_testing[i] = malloc(100000 * sizeof(int));
    }

	//set path_lengths = 1 (important!)
	for (int i = 0; i < orbit_length-1; i++) { 
        path_lengths[i] = 1;
		path_lengths_testing[i] = 1;
    }

	//open postscript file, declare window
	ps_open(argc,argv);
	ps_window(centerx,centery,radius);
	
	//draw stable manifold: this will no longer work since I changed f from sym(x(y(z))) to x(y(z)), and thus f no longer has a fixed point at (0,0,1).
	// draw_manifold(); 

	//draw background
	point starter = lift(.1,.1);
	//draw_background(starter,20000);

	// plot and store finite orbit
	point marked_orbit[orbit_length];
	point marked_point = lift(0.96018, -2.00341);

	//for(int i = 0; i < 100; i++){
		//pprint(marked_point);
		///marked_point = f(marked_point);
	//}
	//abort();
	//lift(1.41050, 0.57300);
	//lift(1.52900, 0.17050);
	//lift(1.00000, 0.50000);
	//lift(2.00100, 0.96100);
	plot_orbit(marked_orbit, marked_point, orbit_length);
	for(int i = 0; i < orbit_length; i++){
		pprint(marked_orbit[i]);
	}
	fprintf(stderr, "done \n");
	//sort finite_orbit
	qsort(marked_orbit, orbit_length, sizeof(point), compare);

	for(int i = 0; i < orbit_length; i++){
		pprint(marked_orbit[i]);
	}
		
	//store x,y coordinates of finite orbit after normalization and stereo projection through 0,0,1
	double x_vals[orbit_length];
	double y_vals[orbit_length];

	for(int i = 0; i < orbit_length; i++){
		x_vals[i] = transform(marked_orbit[i]).x;
		y_vals[i] = transform(marked_orbit[i]).y;
		fprintf(stderr, "( %f, %f)\n", x_vals[i], y_vals[i]);
	}
	
	fprintf(stderr, "Finite orbit stored\n");
	
	//draws straight-line path between q-sorted marked points. Takes its image inder f^2 and record path data
	int example = 0;
	int m = 1;
	
	for(int i = 0; i < orbit_length-1; i++){ //orbit_length-1
		fprintf(stderr, "starting %d\n ", i);
		m = connect_path(marked_orbit[i], marked_orbit[i+1], 1); //1 meaning iterate f once. 
		//for(int j = 0; j<m; j++){
			//fprintf(stderr, "%f\n", transform(ps[j]).x);
		//}
		//abort();
		//for()
		//graph example path
		
		
		draw_path(m-1);
		
		
		 //minus 1 since curt has different conventions
		
		//determine array and modifying path lengths

		determine_arrays(x_vals, y_vals, array_pos[i], array_height[i], &m, marked_orbit, &path_lengths[i]);
		determine_arrays(x_vals, y_vals, array_pos_testing[i], array_height_testing[i], &m, marked_orbit, &path_lengths_testing[i]);  
		
	}
	
	for(int k = 0; k < 20; k++){
		fprintf(stderr, "(%d, %d)\n",pos_array[k], above_array[k]);
	}
	fprintf(stderr, "MADE IT\n");

	for(int k = 0; k < orbit_length-1; k++){
		for(int j = 0; j < path_lengths_testing[k]; j++){
			fprintf(stderr, "(%d, %d)\n", array_pos_testing[k][j], array_height_testing[k][j]);
		}
	}	

	//simplify arrays
	for(int i = 0; i<orbit_length-1; i++){ 
		fully_simplify_arrays(array_pos[i], array_height[i], &path_lengths[i]);
		fully_simplify_arrays(array_pos_testing[i], array_height_testing[i], &path_lengths_testing[i]);
	}
	for(int k = 0; k < orbit_length-1; k++){
		for(int j = 0; j < path_lengths_testing[k]; j++){
			fprintf(stderr, "(%d, %d)\n", array_pos_testing[k][j], array_height_testing[k][j]);
		}
	}	
	//abort();


	//Creates temporary arrays to store path after contracting first through kth vertices
	int arr1_temp[100000]; //size 1000 is safe
	int arr2_temp[100000];
	int size_pos_temp = 0;

	//Creates temporary array to store new twists
	int new_twists[100000];
	int size_temp = 0;


	//TESTING CORNER

	//test 1 passes. Answer: s0.s1.s2.s2.S1.
	/*
	int a_array_pos[60][60] = {
    {1, 2, 3, 3, 2, 1, 0},
    {0, 1, 2},
    {2, 1, 1, 2, 3}
    // The rest will be initialized to 0
	};
	int a_array_height[60][60] = {
		{-1000, 1, 0, 1, 1, 1, -1000},
		{-1000, 0, -1000},
		{-1000, 0, 1, 1, -1000}	
		// The rest will be initialized to 0
	};
	path_lengths[0] = 7;
	path_lengths[1] = 3;
	path_lengths[2] = 5;
	path_lengths_testing[0] = 7;
	path_lengths_testing[1] = 3;
	path_lengths_testing[2] = 5;
	

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 60; j++) {
			array_pos[i][j] = a_array_pos[i][j];
			array_pos_testing[i][j] = a_array_pos[i][j];
			array_height[i][j] = a_array_height[i][j];
			array_height_testing[i][j] = a_array_height[i][j];
		}
	}*/
	//test 2 passes!
	//outputes s_2.s_1.s_0.s_0.s_1.s_2.s_1.s_0.s_0.s_1.s_2.s_1.s_0.s_1.S_2.S_2.s_1
	//solution is S_0s_2s_1S_1s_0S_2s_1
	// check they are the same in SB4 (!!) using:
	/*w = 'S_0s_2s_1S_2s_0S_2s_1'
	q = 's_2.s_1.s_0.s_0.s_1.s_2.s_1.s_0.s_0.s_1.s_2.s_1.s_0.s_1.S_2.S_2.s_1'
	W = S.mapping_class(w)
	Q = S.mapping_class(q)
	if(W == Q): print("yay")*/
	/*
	int a_array_pos[60][60] = {
		{1, 2, 3, 3, 2, 1, 1, 2},
		{2, 1, 1, 2, 3},
		{3, 2, 1, 1, 2, 3, 3, 2, 1, 0}
		// The rest will be initialized to 0
	};
	int a_array_height[60][60] = {
		{-1000, 0, 1, 0, 0, 0, 1, -1000},
		{-1000, 1, 0, 0, -1000},
		{-1000, 0, 0, 1, 0, 1, 0, 0, 0, -1000}	
		// The rest will be initialized to 0
	};
	path_lengths[0] = 8;
	path_lengths[1] = 5;
	path_lengths[2] = 10;
	path_lengths_testing[0] = 8;
	path_lengths_testing[1] = 5;
	path_lengths_testing[2] = 10;
	

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 60; j++) {
			array_pos[i][j] = a_array_pos[i][j];
			array_pos_testing[i][j] = a_array_pos[i][j];
			array_height[i][j] = a_array_height[i][j];
			array_height_testing[i][j] = a_array_height[i][j];
		}
	}
	*/

	/*
	int array[100] = {5, 4, 3, 2, -3, -2};
	int length = 6;
	for(int n = 0; n < length; n++){
		fprintf(stderr, "%d\n", array[n]);
	} 
	simplify_final_array(array, &length);
	for(int n = 0; n < length; n++){
		fprintf(stderr, "new %d\n", array[n]);
	} 
	abort();*/
	/*
	int testing =23 ; 
	//process_array(array, &testing, 2);
	int arr1[100] = {10,9,8,7,6,5,4,3,2,1,1,2,3,4,5,6,7,8,8,7,6,5,4};
	int arr2[100] = {-1000,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,1,-1000};
	for(int n = 0; n < testing; n++){
		fprintf(stderr, "%d, %d\n", arr1[n], arr2[n]);
	} 
	s(arr1, arr2, &testing, 4);
	for(int n = 0; n < testing; n++){
		fprintf(stderr, "%d, %d\n", arr1[n], arr2[n]);
	} 
	
	abort(); */

	//change orbit_length to 4

	//abort();

	//MAIN MACHINERY: iterates through paths, applies twists, stores result in twist array
	fprintf(stderr, "to main machinery:\n");

	for(int i = 0; i < orbit_length-1; i++){ //replace 2 with orbit_length-1
		size_temp = 0;

		// print input data:
		for(int o = 0; o < orbit_length - 1; o++){
			for(int j = 0; j < path_lengths[o]; j++){
				fprintf(stderr, "(%d, %d)\n", array_pos[o][j], array_height[o][j]);
			}
		}
		//applies all previous twists to ith path
		fprintf(stderr, "ITERATION NUMBER %d\n", i);

		for(int j = 0; j < twists_length; j++){
			if(twists[j] >= 0){
				s(array_pos[i], array_height[i], &path_lengths[i], twists[j]);
			}
			if(twists[j] < 0){
				S(array_pos[i], array_height[i], &path_lengths[i], -twists[j]-1); //S not yet defined
			}
		}

		//contracts vertex 0, 1, .., i, stores new data in temporary position/height arrays
		delete_indices(array_pos[i], array_height[i], &path_lengths[i], arr1_temp, arr2_temp, i, &size_pos_temp);
		fully_simplify_arrays(arr1_temp, arr2_temp, &size_pos_temp);

		//stores twist data of contracted path in temporary twist array
		star_algorithm(arr1_temp, arr2_temp, &size_pos_temp, new_twists, &size_temp, i); //orbit_length - 2 > number of contracted strands
				
		//translates twist data on contracted path to twist data on full path
		process_array(new_twists, &size_temp, i);

		//stores new translated twist data in final twist array
		for(int j = 0; j < size_temp; j++){
			twists[twists_length + j] = new_twists[j];
		}
		twists_length += size_temp;

		for(int k = 0; k<=i; k++){
			for(int j = 0; j < size_temp; j++){
				if(new_twists[j] >= 0){
					s(array_pos[k], array_height[k], &path_lengths[k],new_twists[j]);
				}
				if(new_twists[j] < 0){
					S(array_pos[k], array_height[k], &path_lengths[k], -new_twists[j]-1); //S not yet defined}
				}
			}
		}


		size_temp = 0;

		//stores twist data of full path in temporary twist array (just dehn twists around previously contracted path)
		if(path_lengths[i] > 2){ //i.e. not yet reduced
			reverse_array(array_height[i], path_lengths[i]);
			reverse_array(array_pos[i], path_lengths[i]);
		}

		size_temp = 0;
		star_algorithm(array_pos[i], array_height[i], &path_lengths[i], new_twists, &size_temp, i);
		if(size_temp > 0){
			fprintf(stderr, "indside!");
		}
		
		
		//stores new translated twist data in final twist array
		
		for(int j = 0; j < size_temp; j++){
			twists[twists_length + j] = new_twists[j];
		}
		twists_length += size_temp;

		if(size_temp > 0){
			if(twists[twists_length-1] < 0){
				twists[twists_length] = -i-1;
				S(array_pos[i], array_height[i], &path_lengths[i], -i-1);
			} else{
				twists[twists_length] = i;
				s(array_height[i], array_height[i],&path_lengths[i], i);
			} //manually adding last twist.
			twists_length++;
		}
		
		//simplify_final_array(twists, &twists_length);	
	}

	fprintf(stderr, "main machinery complete\n");
	
	simplify_final_array(twists, &twists_length);
	fprintf(stderr, "array is simplified!\n");

	//check machinery
	for(int k = 0; k < orbit_length-1; k++){
		for(int j = 0; j < path_lengths_testing[k]; j++){
			fprintf(stderr, "(%d, %d)\n", array_pos_testing[k][j], array_height_testing[k][j]);
		}	
		for(int j = 0; j < twists_length; j++){
			if(twists[j] >= 0){
				s(array_pos_testing[k], array_height_testing[k], &path_lengths_testing[k], twists[j]);
				fprintf(stderr, "s_%d \n", twists[j]);
			}
			if(twists[j] < 0){
				S(array_pos_testing[k], array_height_testing[k], &path_lengths_testing[k], -twists[j]-1); //S not yet defined
				fprintf(stderr, "S_%d \n", -twists[j]-1);
			}
			for(int j = 0; j < path_lengths_testing[k]; j++){
				fprintf(stderr, "(%d, %d)\n", array_pos_testing[k][j], array_height_testing[k][j]); //ethan added
			}
		}
		fprintf(stderr, "end\n");
		//for(int j = 0; j < path_lengths_testing[k]; j++){
		//	fprintf(stderr, "(%d, %d)\n", array_pos_testing[k][j], array_height_testing[k][j]);
		//}
	}

	int k = 0;

	simplify_final_array(twists, &twists_length);
	fprintf(stderr, "array is simplified!\n");

	fprintf(stderr, "size is %d\n", twists_length);

	//translate final twist array into s_i's. Reverses order! 
	//this either returns mapping class of f or inverse. both have the same entropy
	for(int i = twists_length - 1; i >= 0; i--){
		twists_word[i] = malloc(5); // Allocate memory for "S_i" (including null terminator)
		if(twists[i]>=0){
			sprintf(twists_word[i], "s_%d", twists[i]);
		}
		if(twists[i]<0){
			sprintf(twists_word[i], "S_%d", -twists[i]-1);
		}
	}
	

	//print word!
	fprintf(stderr, "final twists length is %d\n", twists_length);
	//fprintf(stderr, "%s", twists[0]);
	//for (int j = 0; j < twists_length; j++) {
    //    fprintf(stderr, "%s.", twists_word[j]);
    //}
   	fprintf(stderr, "NEW\n");
	for (int j = twists_length-1; j >=0; j--) {
        fprintf(stderr, "%s.", twists_word[j]);
    }


	//close postscript file
	ps_close();

	exit(0);
}

/*TO DO: 
-testing for k = 4, step by step. ENDPOINT is wrong.
-ORDER OF COMPOSITION is off

*/

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
	above_array[0] = -1000; //path_origin //above_array is used as an intermediate.

	//fprintf(stderr, "Determine_arrays - initialized\n");

	for(int t = 1; t < *m; t++){ 
		determine_position(x_vals, y_vals, ps[t]); // if position is between ith and i+1st points, cur_position returns i+1
		//fprintf(stderr, "first position determined\n");
		//fprintf(stderr, "value of b is %d\n", *b);

		//fprintf(stderr, "pos is %d \n", current_pos);
		//fprintf(stderr, "current_pos is %f\n", transform(ps[t]).x);

		if(*b==1){
			pos_array[*b] = current_pos;
			//fprintf(stderr, "%f\n", transform(ps[t]).x);
			above_array[*b] = 0; //this does not matter
			//fprintf(stderr, "current_pos is (%d,%d, %d) \n", pos_array[*b], above_array[*b], *b);
			//pprint(transform(ps[t-1]));
			//pprint(transform(ps[t]));
			//pprint(transform(ps[t+1]));
			(*b)++;
			//fprintf(stderr, "case 1\n");
			continue;
		}
		//fprintf(stderr, "made it");
		if(current_pos - pos_array[(*b)-1] == -1){ //moving left
			fprintf(stderr, "case 2\n");
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
			//fprintf(stderr, "current_pos is (%d,%d) \n", pos_array[*b], above_array[*b]);
			//fprintf(stderr, "current_pos is (%d,%d, %d) \n", pos_array[*b], above_array[*b], *b);
			//pprint(transform(ps[t-1]));
			//pprint(transform(ps[t]));
			//pprint(transform(ps[t+1]));
			(*b)++;
			continue;
		}
		//fprintf(stderr, "made it");
		if(current_pos - pos_array[(*b)-1] == 1){ //moving right
			fprintf(stderr, "case 3\n");
			pos_array[*b] = current_pos;
			//fprintf(stderr, "%d", current_pos);
			if(transform(ps[t]).y > y_vals[pos_array[(*b)-1]]){
				above_array[*b] = 1;
			} else {
				above_array[*b] = 0;
			}
			//fprintf(stderr, "current_pos is (%d,%d, %d) \n", pos_array[*b], above_array[*b], t);
			//pprint(transform(ps[t-1]));
			//pprint(transform(ps[t]));
			//pprint(transform(ps[t+1]));
			(*b)++;
			continue;
		}
		int last = pos_array[(*b)-1];
		if (pos_array[(*b)-1] < current_pos) { //moving right
			//fprintf(stderr, "%d", current_pos);
			fprintf(stderr, "interesting case\n");
			for (int i = last; i <= current_pos-1; i++) { //ethan changed from 'int i = last + 1; i <= current_pos' on jan 16!!
				pos_array[*b] = i+1; //ethan changed from i to i+1
				//fprintf(stderr, "%d", i);
				fprintf(stderr, "%f\n", y_vals[i]);
				if((transform(ps[t-1]).y > y_vals[i]) && (transform(ps[t]).y > y_vals[i])){
					above_array[*b] = 1;
				} else if((transform(ps[t-1]).y < y_vals[i]) && (transform(ps[t]).y < y_vals[i])){
					above_array[*b] = 0;
					fprintf(stderr, "want this twice\n");
				} else{
					fprintf(stderr, "here issue with over/under \n");
					//fprintf(stderr, "(%d,%d) \n", pos_array[*b-1], above_array[*b-1]);
					//pprint(transform(ps[t-1]));
					//pprint(transform(ps[t]));
					//pprint(transform(ps[t+1]));
					//fprintf(stderr, "(%d,%d) \n", pos_array[*b], above_array[*b]);
					//fprintf(stderr, "(%d,%d) \n", pos_array[*b], above_array[*b]);
					abort();
					//pprint(ps[t-1]);
					//pprint(ps[t]);
					//fprintf(stderr, "%f",  y_vals[i]); 
				}
				//fprintf(stderr, "current_pos is (%d,%d) \n", pos_array[*b], above_array[*b]);
				fprintf(stderr, "current_pos is (%d,%d, %d) \n", pos_array[*b], above_array[*b], t);
				pprint(transform(ps[t-1]));
				pprint(transform(ps[t]));
				pprint(transform(ps[t+1]));
				(*b)++;
			}
			continue; //ethan added jan 13th

    	} else if (pos_array[(*b)-1] > current_pos) { //moving left
			for (int i = last - 1; i >= current_pos; i--) {
				pos_array[*b] = i;
				if((transform(ps[t-1]).y > y_vals[i]) && (transform(ps[t]).y > y_vals[i])){
					above_array[*b] = 1;
				} else if((transform(ps[t-1]).y < y_vals[i]) && (transform(ps[t]).y < y_vals[i])){
					above_array[*b] = 0;
				} else{
					fprintf(stderr, "issue with over/under");
					fprintf(stderr, "current_pos is (%d,%d, %d) \n", pos_array[*b], above_array[*b], *b);
					abort();
					//pprint(ps[t-1]);
					//pprint(ps[t]);
					//fprintf(stderr, "%f",  y_vals[i]); 
				}
				//fprintf(stderr, "(%d, %d)\n", pos_array[b], above_array[b]);
				//fprintf(stderr, "current_pos is (%d,%d) \n", pos_array[*b], above_array[*b]);
				//fprintf(stderr, "current_pos is (%d,%d, %d) \n", pos_array[*b], above_array[*b], *b);
				//pprint(transform(ps[t-1]));
				//pprint(transform(ps[t]));
				//pprint(transform(ps[t+1]));
				(*b)++;
			}
			continue; //ethan added jan 13th
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
	height[(*b)-1] = -1000; //why b -1?
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

// arr1 is the first position, arr2 can be either 0 OR 1, only indicates if it's above or below a point
bool simplify_arrays(int *arr1, int *arr2, int *size) {
	//fprintf(stderr, "VALUE OF ARR 1 AND ARR2\n");
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
