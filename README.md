# spiral
This code produces a square plane of random numbers, that is filled starting from the corners. In this way, using a larger grid number results in a larger plane that however has exactly the same random numbers at the corners than the smaller one (if the same initial seed is used).

The code depends on GSL to generate the random seeds, and on MPI since it is parallel.

To compile:

  mpicc -o pspiral.x pspiral.c -lgsl -lgslcblas

should be sufficient.
