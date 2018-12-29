#include "heat.h"
#include <omp.h>
#include <math.h>
/*
 * Function to copy one matrix into another
 */

void copy_mat (double *u, double *v, unsigned sizex, unsigned sizey)
{
	#pragma omp parallel
	{
		int howmany = omp_get_num_threads(); 
		int blockid = omp_get_thread_num();
		int i_start=lowerb(blockid,howmany,sizex);
		int i_end = upperb(blockid,howmany,sizex);
		for (int i=max(1,i_start); i<=min(sizex-2,i_end); i++)
			for (int j=1; j<=sizey-2; j++) 
				v[ i*sizey+j ] = u[ i*sizey+j ];
	}
	/*for (int i=1; i<=sizex-2;i++)
		for(int j=1; j<=sizey-2; j++)
			v [i*sizey+j] = u[i*sizey+j];
	*/
}


/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
	double  sum=0.0;
#pragma omp parallel reduction(+:sum)
	{
		int howmany=omp_get_num_threads();
		//int howmany=4;
		//if(howmany>omp_get_num_threads())howmany=omp_get_num_threads();
			
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

struct abc{
	double sum;
	int done;
};

double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
	double sum=0.0;
	struct abc thread_vector[omp_get_max_threads()];
	for(int i=0;i<omp_get_max_threads();++i) {
		thread_vector[i].sum = 0.0;
		thread_vector[i].done = 0;
	}	
	#pragma omp parallel
	{
	
	int P = omp_get_num_threads();
	int root = sqrt(P);
	int howmanyi=root;
	int howmanyj=P/root; 
	int id = omp_get_thread_num();
	int id_i = id/howmanyj;
	int id_j = id%howmanyj;
	int i_start = lowerb(id_i, howmanyi, sizex);
	int i_end   = upperb(id_i, howmanyi, sizex);
	int j_start = lowerb(id_j, howmanyj, sizey);
	int j_end   = upperb(id_j, howmanyj, sizey);
	//printf("id : %d, id_i: %d, id_j: %d -> %d-%d, %d-%d\n",id,id_i,id_j,i_start,i_end,j_start,j_end);
//	printf("id: %d depends of %d,%d and %d,%d \n",id,i_start-1,j_end,i_end,j_start-1);
	//#pragma omp task depend(out:u[i_end*sizey + j_end]) depend(in:u[sizey*(i_start-1)+j_end]) depend(in:u[sizey*i_end+j_start-1])
	//{
	//int tmp=0;
	while((id>howmanyj && thread_vector[id-howmanyj].done==0) || (id_j>0 && thread_vector[id-1].done==0)){
		__asm__("nop");
	}
//	printf("\n %d",id);
	double temp = 0.0;
	for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
		for (int j=max(1, j_start); j<= min(sizey-2, j_end); j++) {
	        	double unew= 0.25 * ( u[ i*sizey	+ (j-1) ]+  // left
				u[ i*sizey	+ (j+1) ]+  // right
				u[ (i-1)*sizey	+ j     ]+  // top
				u[ (i+1)*sizey	+ j     ]); // bottom
			double diff = unew - u[i*sizey+ j];
			temp += diff*diff;
			u[i*sizey+j]=unew;
		} 
	}//for
	#pragma omp critical
	{
	thread_vector[id].done=1;
	sum += temp;
	}
	//}//task
	}//paralel
	//for(int i=0; i<omp_get_max_threads(); ++i) sum+=thread_vector[i].sum;
	return sum;
}
