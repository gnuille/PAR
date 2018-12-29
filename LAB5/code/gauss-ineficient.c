#include "heat.h"
#include <omp.h>
/*
 * Function to copy one matrix into another
 */

void copy_mat (double *u, double *v, unsigned sizex, unsigned sizey)
{
	for (int i=1; i<=sizex-2; i++)
		for (int j=1; j<=sizey-2; j++) 
			v[ i*sizey+j ] = u[ i*sizey+j ];
}

/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
	double  sum=0.0;
	#pragma omp parallel reduction(+:sum)
{
	int howmany=4;
	//for (int blockid = 0; blockid < howmany; ++blockid) {
	int blockid = omp_get_thread_num();
	int i_start = lowerb(blockid, howmany, sizex);
	int i_end = upperb(blockid, howmany, sizex);
	for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
		for (int j=1; j<= sizey-2; j++) {
			utmp[i*sizey+j]= 0.25 * ( u[ i*sizey     + (j-1) ]+  // left
					u[ i*sizey     + (j+1) ]+  // right
					u[ (i-1)*sizey + j     ]+  // top
					u[ (i+1)*sizey + j     ]); // bottom
			double diff;	
			diff = utmp[i*sizey+j] - u[i*sizey + j];
			sum += diff * diff; 
		}
	}
	//}
}
	return sum;
}

/*
 * Blocked Gauss-Seidel solver: one iteration step
 */
double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
	double unew, diff, sum=0.0;
	int howmany=4;
	#pragma omp parallel for ordered(2) reduction(+:sum) private(diff,unew) 
		for (int i=1; i<= sizex-2; i++) {
			for (int j=1; j<= sizey-2; j++) {
				#pragma omp ordered depend (sink:i,j-1) depend (sink:i-1,j)
				unew= 0.25 * ( u[ i*sizey	+ (j-1) ]+  // left
						u[ i*sizey	+ (j+1) ]+  // right
						u[ (i-1)*sizey	+ j     ]+  // top
						u[ (i+1)*sizey	+ j     ]); // bottom
				diff = unew - u[i*sizey+ j];
				sum += diff * diff; 
				u[i*sizey+j]=unew;
				#pragma omp ordered depend (source)
			}
		}
	return sum;
}
