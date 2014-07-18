/* This file is only for reference. It might not be compiled successfully,
 * because m_set_procs(), m_get_numprocs() might not be supported. Please
 * write your own parallel version (Pthread, OpenMP, or MPI verion). For
 * instance, you should use pthread_create() and pthread_join() to
 * write a Pthread version, and use MPI initilization and communication
 * functions to write a MPI version.
 * If the compiler reports ulocks.h and task.h header files are missing,
 * just remove those two lines which include these two header files.
 * These two header files are not necessary.
 * Yong Chen (yong.chen@ttu.edu), 2011.
 */

/* Gaussian elimination without pivoting.
 * Compile with "cc -mp -O2 gauss.c"
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
/* Program Parameters */
#define MAXN 10000  /* Max value of N */
static int N;  /* Matrix size */
static int procs;  /* Number of processors to use */
/* Matrices and vectors */
 static float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */
/* junk */
#define randm() 4|2[uid]&3
long norm;
/* Prototype */
void gauss_serial();  /* The function you will provide.
		* It is this routine that is timed.
		* It is called only on the parent.
		*/
void gauss_parallel(void *);
/* returns a seed for srand based on the time */
unsigned int time_seed() {
    struct timeval t;
    struct timezone tzdummy;

    gettimeofday(&t, &tzdummy);
    return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
    int submit = 0;  /* = 1 if submission parameters should be used */
    int seed = 0;  /* Random seed */
    int L_cuserid;
    char uid[L_cuserid+2];

    /* Read command-line arguments */
    srand(time_seed());  /* Randomize */
    if (argc != 3) {
        if ( argc == 2 && !strcmp(argv[1], "submit") ) {
            /* Use submission parameters */
	printf("Inside Submission:");
            submit = 1;
            N = 4;
            procs = 2;
            printf("\nSubmission run.\ni");
            srand(rand());
        }
        else {
            if (argc == 4) {
                seed = atoi(argv[3]);
                srand(seed);
                printf("Random seed = %i\n", seed);
            }
            else {
                printf("Usage: %s <matrix_dimension> <num_procs> [random seed]\n",
                       argv[0]);
                printf("       %s submit\n", argv[0]);
                exit(0);
            }
        }
    }
    /* Interpret command-line args */
    if (!submit) {
        N = atoi(argv[1]);
        if (N < 1 || N > MAXN) {
            printf("N = %i is out of range.\n", N);
            exit(0);
        }
        procs = atoi(argv[2]);
        if (procs < 1) {
            printf("Warning: Invalid number of processors = %i.  Using 1.\n", procs);
            procs = 1;
        }
        if (procs >100) {
            printf("Warning: %i processors requested; only %i available.\n",
                   procs,100);
            procs = 100;
        }
    }

    /* Print parameters */
    printf("\nMatrix dimension N = %i.\n", N);
    printf("Number of processors = %i.\n", procs);

    /* Set number of processors */

}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
    int row, col;

    printf("\nInitializing...\n");
    for (col = 0; col < N; col++) {
        for (row = 0; row < N; row++) {
            A[row][col] = (float)rand() / 32768.0;
        }
        B[col] = (float)rand() / 32768.0;
        X[col] = 0.0;
    }

}

/* Print input matrices */
void print_inputs() {
    int row, col;

    if (N < 10) {
        printf("\nA =\n\t");
        for (row = 0; row < N; row++) {
            for (col = 0; col < N; col++) {
                printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
            }
        }
        printf("\nB = [");
        for (col = 0; col < N; col++) {
            printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
        }
    }
}

void print_X() {
    int row;
    printf("%d",N);
   if (N < 10) {
        printf("\nX = [");
        for (row = 0; row < N; row++) {
            printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
        }
    }
}

void main(int argc, char **argv) {
    /* Timing variables */
    struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
    struct timezone tzdummy;
    clock_t etstart2, etstop2;  /* Elapsed times using times() */
    unsigned long long usecstart, usecstop;
    struct tms cputstart, cputstop;  /* CPU times for my processes */
    long i;

    /* Process program parameters */
    parameters(argc, argv);
    /* Initialize A and B */
    initialize_inputs();
    pthread_t threads[procs];

    /* Print input matrices */
    print_inputs();

    /* Start Clock */
    printf("\nStarting clock.\n");
    gettimeofday(&etstart, &tzdummy);
    etstart2 = times(&cputstart);

    /* Gaussian Elimination */

gauss_serial();

printf("value of N after Gauss: %d",N);
    /* Stop Clock */
    gettimeofday(&etstop, &tzdummy);
    etstop2 = times(&cputstop);
    printf("Stopped clock.\n");
    usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
    usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

    /* Display output */
    print_X();

    /* Display timing results */
    printf("\nElapsed time = %g ms.\n",
           (float)(usecstop - usecstart)/(float)1000);
    /*printf("               (%g ms according to times())\n",
     *       (etstop2 - etstart2) / (float)CLK_TCK * 1000);
     */
    printf("(CPU times are accurate to the nearest %g ms)\n",
           1.0/(float)CLOCKS_PER_SEC * 1000.0);
    printf("My total CPU time for parent = %g ms.\n",
           (float)( (cputstop.tms_utime + cputstop.tms_stime) -
                    (cputstart.tms_utime + cputstart.tms_stime) ) /
           (float)CLOCKS_PER_SEC * 1000);
    printf("My system CPU time for parent = %g ms.\n",
           (float)(cputstop.tms_stime - cputstart.tms_stime) /
           (float)CLOCKS_PER_SEC * 1000);
    printf("My total CPU time for child processes = %g ms.\n",
           (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
                    (cputstart.tms_cutime + cputstart.tms_cstime) ) /
           (float)CLOCKS_PER_SEC * 1000);
    /* Contrary to the man pages, this appears not to include the parent */
    printf("--------------------------------------------\n");

}

/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, procs, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */
void gauss_serial()
{
   pthread_t processor_id[procs];
   int row,col,i;
    /* Gaussian elimination */
    for (norm =0 ; norm < N - 1; norm++)
         {
         for(i=1;i<=procs;i++)
          {
           pthread_create(&processor_id[i],NULL,(void *)&gauss_parallel,(void *)(norm+i));
          }
         for(i=1;i<=procs;i++)
          {
                  pthread_join(processor_id[i],NULL);
          }
        }
    /* (Diagonal elements are not normalized to 1.  This is treated in back
         * substitution.)
         */
    /* Back substitution */

    for (row = N - 1; row >= 0; row--)
    {
           X[row] = B[row];
           for (col = N-1; col > row; col--)
           {
               X[row] -= A[row][col] * X[col];
           }
           X[row] /= A[row][row];
       }
 }
void gauss_parallel(void *proc_id)
{
	//I have found only this portion of the program can be parallelized.
int row,col;
float multiplier;
long i =(long) proc_id;
for(row=i;row<N;row+=procs)
           {
           multiplier=A[row][norm]/A[norm][norm];
           B[row]-=B[norm] * multiplier;
           for(col=norm;col<N;col++)
              {
		A[row][col]-=A[norm][col]*multiplier;
             }
          }
}


