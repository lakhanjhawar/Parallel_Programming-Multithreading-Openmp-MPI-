#include <stdio.h>
#include<omp.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#define MAX 10000
unsigned int time_seed() {
    struct timeval t;
    struct timezone tzdummy;

    gettimeofday(&t, &tzdummy);
    return (unsigned int)(t.tv_usec);
}
void main()
{
	struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
		    struct timezone tzdummy;
		    clock_t etstart2, etstop2;  /* Elapsed times using times() */
		    unsigned long long usecstart, usecstop;
		    struct tms cputstart, cputstop;  /* CPU times for my processes */
   static int dim;
  int i,j,k;
   static int sum=0;
    static int a[MAX][MAX],b[MAX][MAX],c[MAX][MAX];
  printf("Specify the dimension value\n");
  scanf("%d",&dim);
  if (dim>10000)
  {
	  printf("enter the range less than 10000");
  }
  //printf("Enter the elements of first matrix\n");
  for (i=0;i<dim;i++)
  {
    for (j=0;j<dim;j++)
    {
      a[i][j]=rand();
    	//scanf("%d",&a[i][j]);
    }
  }
    //printf("Enter the elements of second matrix\n");

    for (i=0;i<dim;i++)
    {
      for (j=0;j<dim;j++)
      {
       b[i][j]=rand();
    	//  scanf("%d",&b[i][j])
      }
    }
    printf("\nStarting clock.\n");
            gettimeofday(&etstart, &tzdummy);
            etstart2 = times(&cputstart);
/*for (i=0;i<dim;i++)
        {
          for (j=0;j<dim;j++)
          {
        	  printf("%d\t",a[i][j]);
          }
printf("\n");
          }
    for (i=0;i<dim;i++)
            {
              for (j=0;j<dim;j++)
              {
            	  printf("%d\t",b[i][j]);
              }
              printf("\n");
              }*/


#pragma omp parallel shared(a,b,c,dim) private(i,j,k)
{
#pragma omp for schedule(static)
for(i=0;i<dim;i++)
{
for(j=0;j<dim;j++)
{
c[i][j]=0;
}
}
	for(i=0;i<dim;i++)
	{
		#pragma omp for schedule(static)
                    for(j=0;j<dim;j++)
			for(k=0;k<dim;k++)

				c[i][j]+=a[i][k]*b[k][j];

	}
}

   /* printf("Product of entered matrices:-\n");

    for (i=0;i<dim;i++)
    {
      for (j=0;j<dim;j++)
      {
        printf("%d\n", c[i][j]);

    }
  }*/
    /* Stop Clock */
          gettimeofday(&etstop, &tzdummy);
          etstop2 = times(&cputstop);
          printf("Stopped clock.\n");
          usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
          usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

          /* Display output */

          /* Display timing results */
          printf("\nElapsed time = %g ms.\n",(float)(usecstop - usecstart)/(float)1000);

}
