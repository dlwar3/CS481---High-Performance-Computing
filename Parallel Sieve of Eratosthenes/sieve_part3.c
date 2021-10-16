/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   long long     count;          /* Local prime count */
   double        elapsed_time;   /* Parallel execution time */
   long long     first;          /* Index of first multiple */
   long long     global_count;   /* Global prime count */
   long long     high_value;     /* Highest value on this proc */
   long long     i;
   int           id;             /* Process ID number */
   long long     index;          /* Index of current prime */
   long long     low_value;      /* Lowest value on this proc */
   char          *marked;        /* Portion of 2,...,'n' */
   char          *privatePrimes; /* private copy of primes 3,...,sqrt(n) */
   long long     n;              /* Sieving from 2, ..., 'n' */
   int           p;              /* Number of processes */
   long long     proc0_size;     /* Size of proc 0's subarray */
   long long     prime;          /* Current prime */
   long long     size;           /* Elements in 'marked' */

   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

   low_value = 2 + id*(n-1)/p;
   high_value = 1 + (id+1)*(n-1)/p;
   size = high_value - low_value + 1;

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (n-1)/p;

   if ((2 + proc0_size) < (long long) sqrt((double) n)) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      exit (1);
   }

   /* Allocate this process's share of the array. */

   marked = (char *) malloc (size);

   /* Allocate private copy of primes 3,...,sqrt(n). */

   long long privateSize = (long long) sqrt((double) n);
   privatePrimes = (char *) malloc (privateSize);

   if (marked == NULL || privatePrimes == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit (1);
   }

   for (i = 0; i < size; i++) marked[i] = 0;
   for (i = 0; i < privateSize; i++) privatePrimes[i] = 0;

   prime = 2;
   do {
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime)) first = 0;
            else first = prime - (low_value % prime);
        }
        
        for (i = first; i < size; i += prime) marked[i] = 1;

        long long privateFirst = prime * prime - 2;
        for (i = privateFirst; i<privateSize; i+= prime) privatePrimes[i] = 1;
        
        while (privatePrimes[++index]);
        prime = index + 2;
   } while (prime * prime <= n);
   count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;
   if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
      0, MPI_COMM_WORLD);

   /* Stop the timer */

   elapsed_time += MPI_Wtime();


   /* Print the results */

   if (!id) {
      printf ("There are %d primes less than or equal to %lld\n",
         global_count, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   MPI_Finalize ();
   return 0;
}