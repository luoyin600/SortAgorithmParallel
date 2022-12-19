#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <omp.h>

#define SIZE 1000000

void parallel_counting_sort(int* arr) {
	int* count;
	int max_value = 0;
	int index = 0;

#pragma omp parallel for reduction(max: max_value)
	for (int i = 0; i < SIZE; i++) {
		max_value = arr[i] > max_value ? arr[i] : max_value;
	}

	count = new int[max_value + 1];
	std::fill(count, count + max_value + 1, 0);

#pragma omp parallel for
	for (int i = 0; i < SIZE; i++) {
#pragma omp atomic
		count[arr[i]]++;
	}

	for (int i = 0; i <= max_value; i++) {
		for (int j = 1; j <= count[i]; j++) {
			arr[index++] = i;
		}
	}
}

bool valid(int* s_arr, int* p_arr) {
	for (int i = 0; i < SIZE; i++) {
		if (s_arr[i] != p_arr[i])
			return false;
	}
	return true;
}

int main() {
	int threadNum = omp_get_max_threads();
	omp_set_num_threads(threadNum);
    size_t size = SIZE;

    int* a = new int[size];
    printf("Creating Random List of %ld elements\n", size);
    for (int i = 0; i < size ; ++i) {
        a[i] =(int) rand() % size;
    }
    printf("Created\n");
    
    double startTime = omp_get_wtime ();
    parallel_counting_sort(a);
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