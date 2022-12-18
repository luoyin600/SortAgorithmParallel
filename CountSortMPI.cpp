#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 1000000

#define RANGE_MIN 0
#define RANGE_MAX SIZE

int key(int item) {
    return item;
}

void array_min_max(const int *array, long long size, int *min, int *max) {
    *min = array[0];
    *max = array[0];

    for (long long i = 0; i < size; i++)
        if (array[i] < *min)
            *min = array[i];
        else if (array[i] > *max)
            *max = array[i];
}

void find_min_max(int *array, long long size, int *min, int *max,
                  int num_proc, int rank, const int gMax, const int gMin)
{
    //int local_min = RANGE_MAX;
    //int local_max = RANGE_MIN;

    int local_min = gMax;
    int local_max = gMin;

    /* Divide the total size evenly among every process. */
    const long long local_size = size / num_proc;

    /*
     * Each process (except the first one) will work on a precise sub-portion of
     * the array, finding a local minimum and local maximum value.
     */
    if (rank > 0)
        array_min_max(&array[(rank - 1) * local_size], local_size, &local_min,
                      &local_max);
    /*
     * Process with rank 0 has to find local_min and local_max for its portion
     * of local_size elements AND every element of the array that was left out
     * from the local_size division.
     */
    else {
        /*
         * This is the number of elements that, when the size is not perfectly
         * divisible by the number of processes, no process would cover.
         * It is, at most, equal to num_proc - 1.
         */
        int num_leftout = size - local_size * num_proc;
        array_min_max(&array[(num_proc - 1) * local_size],
                      local_size + num_leftout, &local_min, &local_max);
    }

    /*
     * Find global min and global max among the local ones and share the result
     * with all processes.
     */
    MPI_Allreduce(&local_min, min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}



void counting_sort(int *array, long long size, int num_proc, int rank, const int gMax, const int gMin) {
    int min = 0;
    int max = 0;

    /* Divide the total size evenly among every process. */
    const long long local_size = size / num_proc;

    find_min_max(array, size, &min, &max, num_proc, rank, gMax, gMin);

    /* Size of the count[] array. */
    const int count_size = max - min + 1;

    /*
     * Each process will operate on its local version of the count[] array.
     * Initialized with all of its items at 0.
     */
    //int *local_count = (int *)safe_alloc(count_size * sizeof(int));

    int *local_count = (int *)malloc(count_size * sizeof(int));
    for (int i = 0; i < count_size; i++)
        local_count[i] = 0;

    /*
     * Each process will only consider the first `local_size` items it finds in
     * the input array starting from an offset that depends on the rank and is,
     * therefore, unique.
     */
    long long local_offset = rank * local_size;
    for (long long i = local_offset; i < local_offset + local_size; i++)
        local_count[key(array[i]) - min] += 1;

    /* ============================== RANK = 0 ============================== */
    if (rank == 0) {
        long long k = 0;

        /*
         * Global (and official) version of the count[] array.
         * Initialized with all of its items at 0.
         */
        int *count = (int *)malloc(count_size * sizeof(int));
        for (int i = 0; i < count_size; i++)
            count[i] = 0;

        /*
         * The size might not always be perfectly divisible by the number of
         * processes; in such cases, a portion of the array might be left out.
         * This is the position in the array where the left out elements start.
         */
        long long index_leftout = local_size * num_proc;
        /* Process 0 has to also consider the left out elements (if any). */
        for (long long i = index_leftout; i < size; i++)
            local_count[key(array[i]) - min] += 1;

        /*
         * Receive all the other local_count[] (except for the one already
         * belonging to process 0) and sum their content to that of count[].
         * 'i' represents the process rank.
         */
        for (int i = 0; i < num_proc; i++) {
            if (i > 0)
                MPI_Recv(local_count, count_size, MPI_INT, i, 2, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            for (int j = 0; j < count_size; j++)
                count[j] += local_count[j];
        }

        /* Final section of the algorithm, not parallelizable. */
        for (int i = min; i < max + 1; i++)
            for (long long j = 0; j < count[i - min]; j++)
                array[k++] = i;

        free(count);
    }

    /* ============================== RANK > 0 ============================== */
    else {
        /* Send local_count[] array back to process 0. */
        MPI_Send(local_count, count_size, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }

    free(local_count);

    /*
     * By the end of the algorithm, the array is only sorted in the process
     * with rank 0; with a call to MPI_Bcast, the sorted copy is sent to all the
     * other processes.
     */
    MPI_Bcast(array, size, MPI_INT, 0, MPI_COMM_WORLD);
}



int main(int argc, char **argv) {

    const long long size = SIZE;
    int *array = (int *)malloc(size * sizeof(int));

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int num_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    double start_timer, finish_timer;
    if (rank == 0) {
        printf("Creating Random List of %d elements\n", size);
        int j = 0;
        for (j = 0; j < size; ++j) {
            array[j] =(int) rand() % size;
        }
        printf("Created\n");
        start_timer = MPI_Wtime();
    }

    counting_sort(array, size, num_proc, rank, size, 0);

    if (rank == 0) {
        finish_timer = MPI_Wtime();
	    printf("Total time for %d Clusters : %2.2f sec \n",num_proc, finish_timer-start_timer);

        bool error = false;
        int i=0;
        for(i=0;i<size-1;i++) { 
            if (array[i] > array[i+1]){
		        error = true;
                printf("error in i=%d \n", i);
            }
        }
        if(error)
            printf("Error..Not sorted correctly\n");
        else
            printf("Correct!\n");  
    }

    MPI_Finalize();
    free(array);

    return EXIT_SUCCESS;
}