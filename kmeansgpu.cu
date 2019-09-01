//************************** kmeansgpu.cu ***************************
//*******************Developed by José M. Cecilia*******************
//************************* October 2018************************


#include "cuda.h"
#include "kmeansgpu.h"
#include <curand.h>
#include "device_functions.h"
#include <curand_kernel.h>
#include <math.h>

#define TAMBLOCK 32

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
	if (code != cudaSuccess) {
		const char * error = cudaGetErrorString(code);
		fprintf(stderr, "GPUassert: %s %s %d\n", error, file, line);
		if (abort) exit(code);
	}
}

__global__ void setup_kernel(curandState *state, unsigned long seed) {
	
    int index = threadIdx.x;
	curand_init(seed, index, 0, &state[index]);
    __syncthreads();
}

__global__ void random_init_centroidCUDA(float * cluster_centro_id, float * dataSetMatrix, int clusters, float rows, int columns, curandState *D_state) {
	
    int tx = threadIdx.x;
    int pos=tx*columns;

	int random = ceil(curand_uniform(&D_state[tx])*rows);

    for (int i=0; i<columns; i++){
    	cluster_centro_id[pos+i] = dataSetMatrix[random+i];
        //printf ("El random es %f para el thread %d\n", cluster_centro_id[pos+i], tx);
    }   
}


__device__ float calc_distances(int dim, float *p1, float  *p2) {

    float distance_sq_sum = 0;

    for (int i = 0; i < dim; ++i)
      distance_sq_sum += sqr(p1[i] - p2[i]);

    return distance_sq_sum;

}


__global__ void calc_all_distancesCUDA(int dim, int k, float *d_X, float *centroid, float *dist) {
	
	int tx = threadIdx.x;
	int pos = tx * k;
	
	for(int i=0; i<k; ++i){
		dist[pos+i] = calc_distances(dim, &d_X[pos*dim], &centroid[i*dim]);
	}

}

__global__ void calc_all_distancesCUDA2(int k, float n, int dim, float *d_X, float *centroid, float * dist) {

	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y * k + threadIdx.y * k;
	int index = y + x;
	if(x < k && y < n){
		dist[index] = calc_distances(dim, &d_X[y*dim], &centroid[x*dim]);
	}
}

__global__ void choose_all_clusters_from_distancesCUDA(float n, int k, float *dist, int *cluster_assignment_index){

	int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x<n){
		int best_index = -1;
		float closest_distance = INFINITY;
	
		for (int i = 0; i < k; i++){
			float cur_distance = dist[x*k+i];
			if(cur_distance < closest_distance){
				best_index = i;
				closest_distance = cur_distance;
			}
		}
	
		cluster_assignment_index[x] = best_index;
	}
}

__global__ void copy_assignment_arrayCUDA(int *src, int *tgt) {

	int x = blockIdx.x * blockDim.x + threadIdx.x;

	tgt[x] = src[x];
}

__global__ void calc_short_distanceCUDA(float n, int k, float *dist, int *cluster_assignment_index, float *d_short_dist){

	int x = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(x<n){
		int active_cluster = cluster_assignment_index[x];
		if (active_cluster != -1){
			d_short_dist[x] = dist[x*k+active_cluster];
		}
	}
}

__global__ void suma_arrayCUDA(float n,float * d_short_dist){

	int dim;
	int num = (int) n;
	if(num%2 != 0){
		dim = (num+1)/2;
	} else{
		dim = num/2;
	}
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int index2 = dim + index;
	
	for(int i=0; i<ceil(log2f(num)); i++){
		if(index<dim){
			d_short_dist[index]+=d_short_dist[index2];
			if(dim%2 != 0){
				dim = (dim+1)/2;
			} else{
				dim = dim/2;
			}
			index2 = dim + index;
			syncthreads();	
		}
	}	
}

void calc_total_distance(float n, int k, float *dist, int *cluster_assignment_index, float *d_short_dist){
	
	dim3 block(ceil(n/TAMBLOCK));
	dim3 thread(TAMBLOCK);	
	calc_short_distanceCUDA<<<block, thread>>> (n, k, dist, cluster_assignment_index, d_short_dist); 	
	
	block.x = ceil(n/2/TAMBLOCK);
	suma_arrayCUDA<<<block, thread>>> (n,d_short_dist);

//	return d_short_dist[0];
}

__global__ void init_cluster_centroid(int dim, int * cluster_member_count, float * new_cluster_centroid){

	int x = threadIdx.x;
	
	cluster_member_count[x] = 0;
	
    	for (int i=0; i < dim; i++){
		new_cluster_centroid[x*dim+i] = 0;
	}
}

__global__ void sum_all_points_of_cluster(float n, int dim, float *X, int *cluster_assignment_index, int *cluster_member_count, float *new_cluster_centroid){

	int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < n){
		int active_cluster = cluster_assignment_index[x];
		atomicAdd(&cluster_member_count[active_cluster], 1);
		for(int i = 0; i < dim; i++){
			atomicAdd(&new_cluster_centroid[active_cluster*dim + i], X[x*dim + i]);  
		}
	}
}

__global__ void media_points_of_cluster(int dim, int *cluster_member_count, float *new_cluster_centroid){

	int x = threadIdx.x;
	
	if(cluster_member_count[x] == 0){
		cluster_member_count[x] = 0.00005;
	}
	for(int i = 0; i < dim; i++){
		new_cluster_centroid[x*dim + i] /= cluster_member_count[x];
	}
}

void calc_cluster_centroidsCUDA(int dim, float n, int k, float *X, int *cluster_assignment_index, float *new_cluster_centroid){

	int * cluster_member_count;
	int memsize = sizeof(int) * k;
	cudaMalloc(&cluster_member_count, memsize);
	
    	dim3 block (k/k);
    	dim3 thread (k);
	init_cluster_centroid <<<block, thread>>> (dim, cluster_member_count, new_cluster_centroid);
	
	block.x = ceil(n/TAMBLOCK);
	thread.x = TAMBLOCK;	
	sum_all_points_of_cluster <<<block, thread>>> (n, dim, X,cluster_assignment_index, cluster_member_count, new_cluster_centroid);
	
	block.x = (k/k);
	thread.x = k;
	media_points_of_cluster <<<block, thread>>> (dim, cluster_member_count, new_cluster_centroid);
}

void mostrar_puntos_clusters(int dim, float n, int k, float *X, int *cluster_assignment_index, float *new_cluster_centroid){
	
	int * cluster_member_count;
	int memsize = sizeof(int) * k;
	cudaMalloc(&cluster_member_count, memsize);
	
    	dim3 block (k/k);
    	dim3 thread (k);
	init_cluster_centroid <<<block, thread>>> (dim, cluster_member_count, new_cluster_centroid);
	
	block.x = ceil(n/TAMBLOCK);
	thread.x = TAMBLOCK;	
	sum_all_points_of_cluster <<<block, thread>>> (n, dim, X,cluster_assignment_index, cluster_member_count, new_cluster_centroid);

	
	int * member = (int *) malloc(memsize);
	cudaMemcpy (member, cluster_member_count, memsize, cudaMemcpyDeviceToHost);


	for(int i = 0; i < k; i++){
		printf("%d - %d\n", i+1, member[i]);
	}
    	cudaMemcpy (cluster_member_count, member, memsize, cudaMemcpyHostToDevice);
		
}

extern "C" int kmeansCUDA(int  dim, float *H_X, float n, int k, float *H_cluster_centroid, int iterations, int *H_cluster_assignment_final) {


	float *d_cluster_centroid, *d_X, *d_dist, *d_short_dist;
	int *d_cluster_assignment_final, *d_cluster_assignment_cur, *d_cluster_assignment_prev;
	int memsize;

	memsize = sizeof(float) * k * dim;    	
	cudaMalloc (&d_cluster_centroid, memsize);
    	
	memsize = sizeof(float) * n * dim;
	cudaMalloc (&d_X, memsize);
    	cudaMemcpy (d_X, H_X, memsize, cudaMemcpyHostToDevice);
	
	memsize = sizeof(float) * n * k;
	cudaMalloc (&d_dist, memsize);

	memsize = sizeof(int) * n;
	cudaMalloc(&d_cluster_assignment_final, memsize);
    	cudaMemcpy (d_cluster_assignment_final, H_cluster_assignment_final, memsize, cudaMemcpyHostToDevice);
	cudaMalloc(&d_cluster_assignment_cur, memsize);
	cudaMalloc(&d_cluster_assignment_prev, memsize);

	memsize = sizeof(float) * n;
	cudaMalloc(&d_short_dist, memsize);

	curandState *devStates;
	cudaMalloc (&devStates, k * sizeof(curandState));
	time_t t;
	time(&t);

    	dim3 block (k/k);
    	dim3 thread (k);
	setup_kernel<<<block, thread>>> (devStates, (unsigned long) t );
    	random_init_centroidCUDA<<<block, thread>>> (d_cluster_centroid, d_X, k, n, dim, devStates);
    
    	//dim3 thread (n);
    	//calc_all_distancesCUDA<<<block, thread >>> (dim, k, d_X, d_cluster_centroid, d_dist);
	block.x = ceil(k/TAMBLOCK);
	block.y = ceil(n/TAMBLOCK);	
	thread.x=TAMBLOCK;
	thread.y=TAMBLOCK;
	calc_all_distancesCUDA2<<<block, thread>>> (k, n, dim, d_X, d_cluster_centroid, d_dist);

	block.x = ceil(n/TAMBLOCK);
	block.y = 1;
	thread.x = TAMBLOCK;
	thread.y = 1;	
	choose_all_clusters_from_distancesCUDA <<<block, thread>>> (n, k, d_dist, d_cluster_assignment_cur);
	
	copy_assignment_arrayCUDA<<<block, thread>>> (d_cluster_assignment_cur, d_cluster_assignment_prev);

	calc_total_distance(n, k, d_dist, d_cluster_assignment_cur, d_short_dist);

	memsize = sizeof(float) * n;
	float * H_short_dist = (float *) malloc(memsize);
	cudaMemcpy (H_short_dist, d_short_dist, memsize, cudaMemcpyDeviceToHost);
     	float prev_totD = H_short_dist[0];
    	cudaMemcpy (d_short_dist, H_short_dist, memsize, cudaMemcpyHostToDevice);
	
	int numVariations = 0;

	for(int batch=0; (batch < iterations); ++batch){
		calc_cluster_centroidsCUDA(dim, n, k, d_X, d_cluster_assignment_cur, d_cluster_centroid);
		calc_total_distance(n, k, d_dist, d_cluster_assignment_cur, d_short_dist);
			
		memsize = sizeof(float) * n;
		cudaMemcpy (H_short_dist, d_short_dist, memsize, cudaMemcpyDeviceToHost);
     		float totD = H_short_dist[0];
    		cudaMemcpy (d_short_dist, H_short_dist, memsize, cudaMemcpyHostToDevice);
		
		if(totD >= prev_totD){
			block.x = ceil(n/TAMBLOCK);
			block.y = 1;
			thread.x = TAMBLOCK;
			thread.y = 1;
			copy_assignment_arrayCUDA <<<block, thread>>>(d_cluster_assignment_prev, d_cluster_assignment_cur);

			time(&t);
    			block.x = (k/k);
    			thread.x = k;
			setup_kernel<<<block, thread>>> (devStates, (unsigned long) t );
    			random_init_centroidCUDA<<<block, thread>>> (d_cluster_centroid, d_X, k, n, dim, devStates);

		}
		else{
			block.x = ceil(n/TAMBLOCK);
			block.y = 1;
			thread.x = TAMBLOCK;
			thread.y = 1;
			copy_assignment_arrayCUDA <<<block, thread>>>(d_cluster_assignment_cur, d_cluster_assignment_prev);
			
			block.x = ceil(k/TAMBLOCK);
			block.y = ceil(n/TAMBLOCK);	
			thread.x=TAMBLOCK;
			thread.y=TAMBLOCK;
			calc_all_distancesCUDA2<<<block, thread>>> (k, n, dim, d_X, d_cluster_centroid, d_dist);
	
			block.x = ceil(n/TAMBLOCK);
			block.y = 1;
			thread.x = TAMBLOCK;
			thread.y = 1;	
			choose_all_clusters_from_distancesCUDA <<<block, thread>>> (n, k, d_dist, d_cluster_assignment_cur);
			
			prev_totD = totD;
		}
	}

	block.x = ceil(n/TAMBLOCK);
	block.y = 1;
	thread.x = TAMBLOCK;
	thread.y = 1;
	copy_assignment_arrayCUDA <<<block, thread>>>(d_cluster_assignment_cur, d_cluster_assignment_final);
	
	
	//mostrar_puntos_clusters(dim, n, k, d_X, d_cluster_assignment_final, d_cluster_centroid);

	//printf("Numero de puntos %f.\n",n);
	memsize = sizeof(int) * n;	
	cudaMemcpy (H_cluster_assignment_final, d_cluster_assignment_final, memsize, cudaMemcpyDeviceToHost);
	
	cudaFree(d_cluster_centroid);
	cudaFree(d_X);
	cudaFree(d_dist);
	cudaFree(d_cluster_assignment_cur);
	cudaFree(d_cluster_assignment_prev);
	cudaFree(d_cluster_assignment_final);
	cudaFree(d_short_dist);
	cudaFree(devStates);

	

    	cudaDeviceReset();
    	return 0;
}























