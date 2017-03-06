// for temporary changes and notes search [*]

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>

typedef long long int map_index;
typedef int map_point[2];             // NOTE: this "limits" the plane to a side length of 2^31 points

#define max(x,y) ((x) > (y) ? (x) : (y))

#define x_ 0
#define y_ 1
#define BL 0  // bottom-left
#define TR 1  // top-right

#define dprintf(LEVEL, ...) if((LEVEL) <= verbose_level) fprintf(stderr, __VA_ARGS__);

#define DEBUG_TASK 0
#define DEBUG_LEVEL 0
int verbose_level = -1;
int Ntasks, Task;
int write_spiral_ordinals = 0;

map_index get_map(map_point);
uint random_seed(int *);
int cmp_map_points(const void *, const void *);
int generate_seeds_subregion(map_point, map_point, uint, uint **, uint, map_index *);
int generate_seeds(int, int, map_point [4][2], uint seed, uint **, map_index *);
void transpose_subregion(int, uint, map_point *, map_point *);
int get_plane_subregions(uint, map_point, map_point, map_point [4][2]);
int get_my_region(int, map_point *, map_point *);


map_index get_map(map_point P)
// returns the ordinal number of a given point along the spiral
{
  map_index l, d;
  uint mx, my, c;
  mx = abs(P[x_]);
  my = abs(P[y_]);
  l = 2 * max(mx, my);

  c = (P[y_] > P[x_]) + (P[x_] > 0)*(P[x_] == P[y_]);
  
  d = (c) ? l * 3 + P[x_] + P[y_] : l - P[x_] - P[y_];

  return (l-1)*(l-1) + d;
}


uint random_seed(int *nitems)
// written by R. G. Brown, Duke Univ.
// http://www.sourceware.org/ml/gsl-discuss/2004-q1/msg00071.html
// Generate a random seed from the system
// Called if the user does not provide a seed
{

 unsigned int seed;
 size_t nread;
 struct timeval tv;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
 } else {
   nread = fread(&seed, sizeof(seed), 1, devrandom);
   fclose(devrandom);
 }

 if(nitems != NULL)
   *nitems = nread;
 return(seed);
}


map_index *map;

int cmp_map_points(const void *A, const void *B)
// comparison function to be called by qsort
// performs indexed-sorting
// this version works if map (that is the indexed
//  array) is in its scope
{
  map_index a = map[*(uint*)A];
  map_index b = map[*(uint*)B];

  if(a < b)
    return -1;
  if(a > b)
    return 1;
  return 0;
}



int generate_seeds_subregion(map_point bottom_left, map_point top_right, uint seed, uint **seed_plane, uint offset, map_index *value)
// this function manages the generation of seeds for a single subregion
// - bottom_left, top_right are the corners of the subregion
// - seed is the initial seed for the random number generator
// - seed_plane is a pointer to the pointer that will host the seeds
// - offset is the offset to be accounted while filling the memory pointed by seed_plane
// - xstride is the xlength of the whole memory matrix pointed by seed_plane, it's used to fill it correctly
// - value is a returning value in case of error. Basically it handles the case in which some allocation fails, and
//   it will return the amount of bytes that were requested.
// we assume that the seed plane is stored in memory startign from bottom_left  

{
  map_index N, rnd_counter, Ncnt;
  map_point p;
  uint xlength, ylength, i, j;
  uint *idxs;

  // calculate the region's dimensions
  xlength = top_right[x_] - bottom_left[x_];
  ylength = top_right[y_] - bottom_left[y_];
  N = xlength * ylength;

  // allocate temporary memory for the random seeds
  if((map = calloc( sizeof(map_index), (size_t)N)) == NULL)
    {
      if(value != NULL)
	*value = N * sizeof(map_index);
      return 1;
    }

  // allocate temporary memory for the indexes that will be used to
  // sort the points along the spiral
  if((idxs = calloc( sizeof(uint), (size_t)N)) == NULL)
    {
      if(value != NULL)
	*value = N * sizeof(uint);      
      return 2;
    }

    
  // initialize indexes
  // for the following loop check whether:
  // - leave unroll to the compiler
  // - or explicitly use simd instruction
  for(i = 0; i < N; i++)
    idxs[i] = i;

  dprintf(1, "\t\t-> getting the map..");
  // get the map value from spiral for each point in the region
  for(j = 0; j < ylength; j++)
    {
      p[y_] = bottom_left[y_]+j;
      
      for(i = 0; i < xlength; i++)
	{
	  p[x_] = bottom_left[x_]+i;
	  map[j*xlength + i] = get_map(p);
	}
    }
  dprintf(1, "\n"); fflush(stdout);

  dprintf(1, "\t\t-> sorting %lu spiral points..", (size_t)N);
  // index-sort the map
  qsort(idxs, (size_t)N, sizeof(uint), cmp_map_points);
  dprintf(1, "\n"); fflush(stdout);
  
  // generate random numbers following the spiral traversal order

  if(!write_spiral_ordinals)
    {      
      // set the random number generator
      gsl_rng *R;
      
      R = gsl_rng_alloc(gsl_rng_mt19937);
      gsl_rng_set(R, seed);
      
      // actually generates number
      rnd_counter = 0;
      Ncnt = 0;
      dprintf(1, "\t\t-> generating random seeds..");
      for(i = 0; i < N; i++)
	{
	  if(rnd_counter < map[idxs[i]]-1)
	    for(; rnd_counter < map[idxs[i]]-1; rnd_counter++)
	      gsl_rng_get(R);
	  
	  map[idxs[i]] = (map_index)gsl_rng_get(R);
	  Ncnt++;
	  rnd_counter++;
	}
      dprintf(1, "(%llu randoms for %llu points)\n", rnd_counter, Ncnt); fflush(stdout);
      
      gsl_rng_free(R);  // releases generator
    }
  
  free(idxs); // releases indexes

  // allocate actual memory for the random seeds, if needed
  // The difference with map is that we do not need to store 8-bytes
  // integers, while we needed it to store in the map array the
  // spiral ordered points
  if(*seed_plane == NULL)
    {
      if((*seed_plane = calloc( sizeof(uint), (size_t)N)) == NULL)
	{
	  if(value != NULL)
	    *value = N * sizeof(uint);
	  return 3;
	}
      offset = 0;
    }

  dprintf(1, "\t\t-> copying seeds..");
  //copy the seeds from map into seed_plane, as uint instead of as long long
  uint k;
  for(k = 0, j = 0; j < N; j++, k++)
    {
      if((k > 0) && (j % xlength == 0))
	k += offset;
      (*seed_plane)[k] = (uint)map[j];
    }
  dprintf(1, "\n\t\tdone\n\n"); fflush(stdout);
  
  free(map);
  return 0;
}


int generate_seeds(int N, int Nregions, map_point Regions[4][2], uint seed, uint **seed_plane, map_index *value)
// this function manages the generation of seed for all the subregions
// - N is the side of the whole plane
// - Nregions is the number of subregions
// - Regions is the array of (bottom_left, top_right) coordinates for the subregions
// - seed is the initial seed for the random number generator
// - seed_plane is a pointer to the pointer that will host the seeds
// - value is a returning value in case of error. Basically it handles the case in which some allocation fails, and
//   it will return the amount of bytes that were requested.
  
{
  map_point transposed_subregion[2];
  int i;
  map_index Npoints;
  uint xstride, offset, start, *start_pointer;

  // calculate the total length of the whole region along x
  xstride = Regions[Nregions-1][TR][x_] - Regions[0][BL][x_];
  
  // calculate how many points will be generated
  Npoints = xstride * (Regions[Nregions-1][TR][y_] - Regions[0][BL][y_]);

  // allocate memory for all the seeds
  if((*seed_plane = calloc( sizeof(uint), (size_t)Npoints)) == NULL)
    {
      if(value != NULL)
	*value = Npoints * sizeof(uint);
      return 3;
    }

  // generate seeds in each subregion
  for(i = 0; i < Nregions; i++)
    {
      transposed_subregion[0][x_] = Regions[i][0][x_];
      transposed_subregion[0][y_] = Regions[i][0][y_];
      transposed_subregion[1][x_] = Regions[i][1][x_];
      transposed_subregion[1][y_] = Regions[i][1][y_];

      // transpose subregion in spiral's coordinates
      transpose_subregion(0, N, &transposed_subregion[0], &transposed_subregion[1]);

      if(i & 1) // index is odd -> subregion is on the left part
	offset = Regions[0][TR][x_] - Regions[0][BL][x_]; // offset the insertion on the left by the right part
      else
	offset = Regions[1][TR][x_] - Regions[1][BL][x_];

      start = 0;
      if(i & 2) // index is >= 2 -> subregion is on the higher part
	start = xstride * (Regions[0][TR][y_] - Regions[0][BL][y_]); // offset the starting point by the vertical extent of lower part

      if(i & 1) // index is odd -> subregion is on the right part
	start += Regions[0][TR][x_] - Regions[0][BL][x_]; // further offset the starting point by the horizontal extent of the right part

      start_pointer = &(*seed_plane)[start];
      generate_seeds_subregion(transposed_subregion[0], transposed_subregion[1], seed, &start_pointer, offset, value);
    }

  return 0;
}

void transpose_subregion(int mode, uint N, map_point *bottom_left, map_point *top_right)
// by construction single regions are supposed to lie within quadrants
{
  uint N_2 = N/2;

  if(mode == 0)
    // from plane to spiral
    {     
      if((*bottom_left)[x_] >= N_2)
	{
	  (*bottom_left)[x_] = (*bottom_left)[x_] - N;
	  (*top_right)[x_] = (*top_right)[x_] - N;
	}
      
      if((*bottom_left)[y_] >= N_2)
	{
	  (*bottom_left)[y_] = (*bottom_left)[y_] - N;
	  (*top_right)[y_] = (*top_right)[y_] - N;
	}
    }
  else
    // from spiral back to plane
    {
      if((*bottom_left)[x_] < 0)
	{
	  (*bottom_left)[x_] = N + (*bottom_left)[x_];
	  (*top_right)[x_] = N + (*top_right)[x_];
	}
      
      if((*bottom_left)[y_] < 0)
	{
	  (*bottom_left)[y_] = N + (*bottom_left)[y_];
	  (*top_right)[y_] = N + (*top_right)[y_];
	}
    }
  
  return;
}


int get_plane_subregions(uint N, map_point bottom_left, map_point top_right, map_point Regions[4][2])
{
#define NCOORDS 2
  int c[NCOORDS] = {0, 1};
  uint N_2 = N/2;
  int i, j, k, R, m;
  int Nregions = 1;
  
  for(i = 0; i < NCOORDS; i++)
    {
      Regions[0][BL][i] = bottom_left[i];
      if(bottom_left[i] > N)
	return -1;
    }
  for(i = 0; i < NCOORDS; i++)
    {
      Regions[0][TR][i] = top_right[i];
      if(top_right[i] > N)
	return -1;
    }

  for(k = 0; k < NCOORDS; k++)
    {
      m = 1;
      for(R = 0; R < Nregions; R++)
	if(Regions[R][BL][c[k]] < N_2 &&
	   Regions[R][TR][c[k]] > N_2)
	  {
	    j = Nregions + R;

	    for(i = 0; i < NCOORDS; i++)
	      Regions[j][BL][c[i]] = Regions[R][BL][c[i]];
	    Regions[j][BL][c[k]] = N_2;
	    
	    for(i = 0; i < NCOORDS; i++)
	      Regions[j][TR][c[i]] = Regions[R][TR][c[i]];
	    Regions[R][TR][c[k]] = N_2;


	    m = 2;
	  }

      Nregions *= m;
    }
  

  return Nregions;  
}


int get_my_region(int N, map_point *bottom_left, map_point *top_right)
{
  // give preference to splitting x direction, so that to have slabs

#define shape_factor_limit 100
#define shape_factor 10

  int shape_limit;
  int i, j, mysplit;
  int ysplits = 1;
  int xsplits = 1;
  int *ntasks_for_ysplit;

  if((shape_limit = N / shape_factor) < shape_factor_limit)
    shape_limit = shape_factor_limit;
  
  xsplits = Ntasks;
  while((N / xsplits < shape_limit) && (xsplits > 1))
    xsplits--;

  
  /* while(Ntasks / ysplits > N / shape_factor) */
  /*   ysplits++; */

  ysplits = Ntasks / xsplits + (Ntasks % xsplits > 0);
  
  if(Task == DEBUG_TASK)
    dprintf(DEBUG_LEVEL, "Task %2d: xsplits - ysplits found to be %d - %d\n", Task, xsplits, ysplits);
  
  ntasks_for_ysplit = malloc(sizeof(int)*ysplits);

  for(j = i = 0; i < ysplits; i++)
    j += (ntasks_for_ysplit[i] = xsplits); //Ntasks / ysplits);
  
  if(j < Ntasks)
    ntasks_for_ysplit[ysplits-1] += (Ntasks-j);
  if(j > Ntasks)
    ntasks_for_ysplit[ysplits-1] -= (j - Ntasks);

  if(Task == DEBUG_TASK)
    for(j = i = 0; i < ysplits; i++)
      dprintf(DEBUG_LEVEL, "\tTask %2d: tasks in ysplit[%d] found to be %d\n", Task, i, ntasks_for_ysplit[i]);

  mysplit = Task / xsplits; //(Ntasks / ysplits);
  if(mysplit == ysplits)
    mysplit = ysplits-1;

  int me = Task - (Ntasks / ysplits)*mysplit;
  int step = N / ntasks_for_ysplit[mysplit];
  
  (*bottom_left)[x_] = step * me;
  (*top_right)[x_] = (*bottom_left)[x_] + step;

  if( me == ntasks_for_ysplit[mysplit]-1)
    (*top_right)[x_] = N;

  step = N / ysplits;
  (*bottom_left)[y_] = step * mysplit;
  (*top_right)[y_] = (*bottom_left)[y_] + step;

  if(mysplit == ysplits-1)
    (*top_right)[y_] = N;
  
  return 0;
}


int main(int argc, char **argv)
{
  int Nr, i;
  uint N, *SEEDS, seed;
  map_point bottom_left, top_right;
  map_point subregions[4][2];
  map_index error;

  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &Task);

  //  Ntasks = 1;
  //Task = 0;
  
  N = atoi(*(argv+1));
  if(argc > 2)
    verbose_level = atoi(*(argv+2));
  if(argc > 3)
    {
      int myseed = atoi(*(argv+3));

      if(myseed < 0)
	{
	  write_spiral_ordinals = 1;
	  myseed = -myseed;
	}
      seed = (uint)myseed;      
    }
  else
    {
      if(Task == 0)
	{
	  seed = random_seed(NULL);      
	  printf("using seed: %u\n", seed);
	}
      MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

  get_my_region(N, &bottom_left, &top_right);

  for(i = 0; i < Ntasks; i++)
    {
      if(Task == i)
	{
	  dprintf(1, "[TASK %2d ]\n"
		  "\tmy region is: (%d , %d) (%d , %d) => %u points\n",
		  Task,
		  bottom_left[x_], bottom_left[y_],
		  top_right[x_], top_right[y_],
		  (top_right[x_] - bottom_left[x_] +1)*(top_right[y_] - bottom_left[y_] +1));
	  fflush(stdout);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }

  Nr = get_plane_subregions(N, bottom_left, top_right, subregions);

  dprintf(0, "%d subregion%c: \n\n", Nr, (Nr>1)?'s':' ');
  
  for(i = 0; i < Nr; i++)
    {
      dprintf(0, "\t %d: bottom left : %d, %d\n", i, subregions[i][BL][0], subregions[i][BL][1]);
      dprintf(0, "\t %d: top right : %d, %d\n", i, subregions[i][TR][0], subregions[i][TR][1]);
      dprintf(1, "\t which becomes:\n");

      // transpose to check that transposition is okbottom_left[x_], bottom_left[y_]bottom_left[x_], bottom_left[y_]
      transpose_subregion(0, N, &subregions[i][BL], &subregions[i][TR]);
      dprintf(1, "\t\t bottom left : %d, %d\n", subregions[i][BL][0], subregions[i][BL][1]);
      dprintf(1,"\t\t top right : %d, %d\n", subregions[i][TR][0], subregions[i][TR][1]);

      // transpose back, because transposition will be made by the seeds generating routins
      transpose_subregion(1, N, &subregions[i][BL], &subregions[i][TR]);

      dprintf(0, "\n");
    }


  i = generate_seeds(N, Nr, subregions, seed, &SEEDS, &error);


  if(i > 0)
    printf("*** Task %d :: some problem occured!\n", Task);

  FILE *file;
  char filename[100];

  int k, j;
  int myidx_j, mystride;
  int idx_j;
  
  sprintf(filename, "output_%03d", Task);
  
  for(k = 0; k < Ntasks; k++)
    {
      if(Task == k)
	{
	  file = fopen(filename, "w");
	  mystride = top_right[x_] - bottom_left[x_];
	  
	  for(j = bottom_left[y_]; j < top_right[y_]; j++)
	    {
	      myidx_j = (j - bottom_left[y_]) * mystride;
	      idx_j = j * N;
	      
		for(i = bottom_left[x_]; i < top_right[x_]; i++)
		  fprintf(file,"%u %u %u %u\n", idx_j+i, i, j, SEEDS[myidx_j + (i-bottom_left[x_])]);
	    }
	  fclose(file);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }

  free(SEEDS);

  MPI_Barrier(MPI_COMM_WORLD);
  if(Task == 0)
    {
      char command[200], outname[100], outname2[100];
      int ret;

      if(!write_spiral_ordinals)
	{
	  sprintf(outname, "output.%d_points.%u_seed.%d_tasks", N, seed, Ntasks);
	  sprintf(command, "cat output_* | sort -k1 -n > %s", outname);
	}
      else
	{
	  sprintf(outname, "output.%d_points.%u_seed.%d_tasks.spiral_order", N, seed, Ntasks);
	  sprintf(command, "cat output_* | sort -k4 -n > %s", outname);
	}
      ret = system(command);
      ret = system("rm -f output_*");
    }
  
  MPI_Finalize();
  
  return 0;
}



/* **********************************

How to plot data using gnuplot

- to plot the first N data accordingly to their x,y coordinates
plot $filename u 2:3 every ::::N w p pt 4 ps 0.5 notitle 

- to plot the first N data accordingly to their x,y coordinates normalized to 1
plot $filename u ($2/Np):($3/Np) every ::::N w p pt 4 ps 0.5 notitle 

where Np is the grid number

- to compare datafiles with different Np (say Np1, Np2, Np1 > Np2)
plot $filename1 u ($2/Np1):($3/Np1) every ::::N w p pt 4 ps 0.5 notitle, $filename2 u ($2/Np2):($3/Np2) every ::::N w p pt 4 ps 0.5 notitle


************************************ */
