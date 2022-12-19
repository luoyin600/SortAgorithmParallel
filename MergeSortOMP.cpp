#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <omp.h>

#define SIZE 1000000

#define SMALL 32

extern double get_time (void);
void merge (int a[], int size, int temp[]);
void insertion_sort (int a[], int size);
void mergesort_serial (int a[], int size, int temp[]);
void mergesort_parallel_omp (int a[], int size, int temp[], int threads);
void run_omp (int a[], int size, int temp[], int threads);

// OpenMP merge sort with given number of threads
void mergesort_parallel_omp (int a[], int size, int temp[], int threads)
{
  if (threads == 1)
    {
      mergesort_serial (a, size, temp);
    }
  else if (threads > 1)
    {
#pragma omp parallel sections
      {
//                      printf("Thread %d begins recursive section\n", omp_get_thread_num());
#pragma omp section
	{			//printf("Thread %d begins recursive call\n", omp_get_thread_num());
	  mergesort_parallel_omp (a, size / 2, temp, threads / 2);
	}
#pragma omp section
	{			//printf("Thread %d begins recursive call\n", omp_get_thread_num());
	  mergesort_parallel_omp (a + size / 2, size - size / 2,
				  temp + size / 2, threads - threads / 2);
	}
      }
      // Thread allocation is implementation dependent
      // Some threads can execute multiple sections while others are idle 
      // Merge the two sorted sub-arrays through temp
      merge (a, size, temp);
    }
  else
    {
      printf ("Error: %d threads\n", threads);
      return;
    }
}

void mergesort_serial (int a[], int size, int temp[]){
    // Switch to insertion sort for small arrays
    if (size <= SMALL){
        insertion_sort (a, size);
        return;
    }
    mergesort_serial (a, size / 2, temp);
    mergesort_serial (a + size / 2, size - size / 2, temp);
    // Merge the two sorted subarrays into a temp array
    merge (a, size, temp);
}

void merge (int a[], int size, int temp[]){
    int i1 = 0;
    int i2 = size / 2;
    int tempi = 0;
    while (i1 < size / 2 && i2 < size){
        if (a[i1] < a[i2]){
            temp[tempi] = a[i1];
            i1++;
        }
        else{
            temp[tempi] = a[i2];
            i2++;
        }
        tempi++;
    }
    while (i1 < size / 2){
        temp[tempi] = a[i1];
        i1++;
        tempi++;
    }
    while (i2 < size){
        temp[tempi] = a[i2];
        i2++;
        tempi++;
    }
    // Copy sorted temp array into main array, a
    memcpy (a, temp, size * sizeof (int));
}

void insertion_sort (int a[], int size){
    int i;
    for (i = 0; i < size; i++){
        int j, v = a[i];
        for (j = i - 1; j >= 0; j--){
            if (a[j] <= v)
            break;
            a[j + 1] = a[j];
        }
        a[j + 1] = v;
    }
}

// Driver
void run_omp (int a[], int size, int temp[], int threads)
{
    // Enable nested parallelism, if available
    omp_set_nested (1);
    // Parallel mergesort
    mergesort_parallel_omp (a, size, temp, threads);
}

int main()
{
    // Get arguments
    int size = 1000000;	// Array size 
    int threads = 4;	// Requested number of threads
    // Check nested parallelism availability
    omp_set_nested (1);
    if (omp_get_nested () != 1){
        puts ("Warning: Nested parallelism desired but unavailable");
    }

    int threadNum = omp_get_max_threads();
    //printf("threadNum is %d\n", threadNum);

    omp_set_num_threads(threadNum);

    // Array allocation
    int *a = (int*)malloc (sizeof (int) * size);
    int *temp = (int*)malloc (sizeof (int) * size);
  
    printf("Creating Random List of %d elements\n", SIZE);
    for (int i = 0; i < size ; ++i) {
        a[i] =(int) rand() % size;
    }
    printf("Created\n");

    // Sort
    double startTime = omp_get_wtime ();
    run_omp (a, size, temp, threads);
    double finishTimer = omp_get_wtime ();
  
    printf("Total time for %d Clusters : %2.2f sec \n",threadNum, finishTimer-startTime);

    printf("Checking.. \n");
    bool error = false;
    for(int i=0;i<size-1;i++) {
        if (a[i] > a[i+1]){
            error = true;
            printf("error in i=%d \n", i);
        }
    }
    if(error)
        printf("Error..Not sorted correctly\n");
    else
        printf("Correct!\n");  
    return 0;
}