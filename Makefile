all:
	mpic++ QuickSortMPI.cpp -o QuickSortMPI
	#mpic++ BubbleSortMPI.cpp -o BubbleSortMPI
	mpic++ MergeSortMPI.cpp -o MergeSortMPI
	#mpic++ SelectionSortMPI.cpp -o SelectionSortMPI
	mpic++ CountSortMPI.cpp -o CountSortMPI
	g++ QuickSortOMP.cpp -o QuickSortOMP -fopenmp
	g++ MergeSortOMP.cpp -o MergeSortOMP -fopenmp
	g++ CountSortOMP.cpp -o CountSortOMP -fopenmp
clean:
	rm QuickSortMPI BubbleSortMPI MergeSortMPI CountSortMPI
