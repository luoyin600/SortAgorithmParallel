all:
	mpic++ QuickSortMPI.cpp -o QuickSortMPI
	#mpic++ BubbleSortMPI.cpp -o BubbleSortMPI
	mpic++ MergeSortMPI.cpp -o MergeSortMPI
	#mpic++ SelectionSortMPI.cpp -o SelectionSortMPI
	mpic++ CountSortMPI.cpp -o CountSortMPI
clean:
	rm QuickSortMPI BubbleSortMPI MergeSortMPI CountSortMPI
