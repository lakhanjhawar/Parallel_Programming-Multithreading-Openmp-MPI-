#include<mpi.h>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
double **A, **B, **C;
double *a, *b, *c, *temple_a, *temple_b;
int dg, dl, dl2,p, sp;
int my_rank, my_row, my_col;
MPI_Status status;
void initializer()
{
   int i,j;
    srand((unsigned int)time(NULL));
    for(i=0; i<dg; i++)
      for(j=0; j<dg ; j++)
	  {
    	    A[i][j] = rand();
            B[i][j] = rand();
            C[i][j] = 0.0;
	  }
}
void replacer()
{
   int i,j,k,l;
   int p_imin,p_imax,p_jmin,p_jmax;
   for(k=0; k<p; k++)
   {
	  p_jmin = (k % sp) * dl;
  	  p_jmax = (k % sp + 1) * dl-1;
	  p_imin = (k - (k % sp))/sp * dl;
	  p_imax = ((k - (k % sp))/sp +1) *dl -1;
               l = 0;
      for(i=p_imin; i<=p_imax; i++)
      {
      	  for(j=p_jmin; j<=p_jmax; j++)
      	  {
          temple_a[l] = A[i][j];
	      temple_b[l] = B[i][j];
	      l++;
          }
      }
      if(k==0)
      {
         memcpy(a, temple_a, dl2 * sizeof(double));
	     memcpy(b, temple_b, dl2 * sizeof(double));
      } else
      {
      MPI_Send(temple_a, dl2, MPI_DOUBLE, k, 1, MPI_COMM_WORLD);
	  MPI_Send(temple_b, dl2, MPI_DOUBLE, k, 2, MPI_COMM_WORLD);
      }
   }
}
void process_cleaner()
{
   int i,j,i2,j2,k;
   int p_imin,p_imax,p_jmin,p_jmax;


   for (i=0;i<dl;i++)
	 for(j=0;j<dl;j++)
	   C[i][j]=c[i*dl+j];

   for (k=1;k<p;k++)
   {

       MPI_Recv(c, dl2, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);

       p_jmin = (k % sp) *dl;
       p_jmax = (k % sp + 1) *dl-1;
       p_imin =  (k - (k % sp))/sp     *dl;
       p_imax = ((k - (k % sp))/sp +1) *dl -1;

       i2=0;

       for(i=p_imin; i<=p_imax; i++)
       {
           j2=0;
           for(j=p_jmin; j<=p_jmax; j++)
           {
               C[i][j]=c[i2*dl+j2];
               j2++;
           }
           i2++;
       }
   }
}

void MatrixMultiply(int n, double *a, double *b, double *c)
{
    int i, j, k;

    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            for (k=0; k<n; k++)
                c[i*n+j] += a[i*n+k]*b[k*n+j];
}
void cannon_matrix_multiply(int n, double *a, double *b, double *c, MPI_Comm comm)
{
      int i;
      int nlocal;
      int npes, dims[2], periods[2];
      int myrank, my2drank, mycoords[2];
      int uprank, downrank, leftrank, rightrank, coords[2];
      int shiftsource, shiftdest;
      MPI_Status status;
      MPI_Comm comm_2d;
      MPI_Comm_size(comm, &npes);
      MPI_Comm_rank(comm, &myrank);
      dims[0] = dims[1] = sqrt(npes);
      periods[0] = periods[1] = 1;
      MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);
      MPI_Comm_rank(comm_2d, &my2drank);
      MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);
      MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
      MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);
      nlocal = n/dims[0];
      MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
      MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_DOUBLE, shiftdest,
          1, shiftsource, 1, comm_2d, &status);
      MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
      MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_DOUBLE,
          shiftdest, 1, shiftsource, 1, comm_2d, &status);
    for (i=0; i<dims[0]; i++)
    {
        MatrixMultiply(nlocal, a, b, c);
        MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_DOUBLE,
        leftrank, 1, rightrank, 1, comm_2d, &status);
        MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_DOUBLE,
        uprank, 1, downrank, 1, comm_2d, &status);
    }
      MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource, &shiftdest);
      MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_DOUBLE,
          shiftdest, 1, shiftsource, 1, comm_2d, &status);
      MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource, &shiftdest);
      MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_DOUBLE,
          shiftdest, 1, shiftsource, 1, comm_2d, &status);
      MPI_Comm_free(&comm_2d);
}

int main(int argc, char *argv[])
{
   int i,j;
   double begin, last;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   MPI_Request reqs[4];
   sp = sqrt(p);


   if (sp*sp != p)
   {
      if (my_rank == 0)
	  printf("Runs only for perfect square \n");
      MPI_Finalize();
      exit(1);
   }

   if (argc != 2)
   {
      if (my_rank == 0)
          printf("");
      MPI_Finalize();
      exit(1);
   }
   dg  = atoi(argv[1]);
   dl  = dg / sp;
   dl2 = dl * dl;
   my_col =  my_rank % sp ;
   my_row = (my_rank-my_col) / sp ;
   a = (double *)malloc( dl2 * sizeof(double) );
   b = (double *)malloc( dl2 * sizeof(double) );
   c = (double *)malloc( dl2 * sizeof(double) );
   for(i=0; i<dl2 ; i++)
     c[i] = 0.0;
   temple_a = (double *)malloc( dl2 * sizeof(double) );
   temple_b = (double *)malloc( dl2 * sizeof(double) );

  if (my_rank == 0)
       begin = MPI_Wtime();

  if (my_rank == 0)
   {

      A = (double **)malloc( dg * sizeof(double*) );
      B = (double **)malloc( dg * sizeof(double*) );
      C = (double **)malloc( dg * sizeof(double*) );

      for(i=0; i<dg; i++)
      {
         A[i] = (double *)malloc( dg * sizeof(double) );
         B[i] = (double *)malloc( dg * sizeof(double) );
         C[i] = (double *)malloc( dg * sizeof(double) );
      }
      initializer();
      replacer();
   }
   else
   {
       MPI_Irecv(a, dl2, MPI_DOUBLE, 0 , 1, MPI_COMM_WORLD, &reqs[0]);
       MPI_Irecv(b, dl2, MPI_DOUBLE, 0 , 2, MPI_COMM_WORLD, &reqs[1]);
       for(j=0;j<2;j++)
      {
        MPI_Wait(&reqs[j],&status);
      }
   }

cannon_matrix_multiply(dg,a, b,c, MPI_COMM_WORLD);

   if(my_rank == 0)
   {
     process_cleaner();
  }
   else
   {
      MPI_Isend(c,dl2,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&reqs[2]);
      MPI_Wait(&reqs[2],&status);
   }

   MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank==0)
      {
       last = MPI_Wtime();
       printf("Total time taken: %f\n", last - begin);
      }

   MPI_Finalize();

   return 0;
}
