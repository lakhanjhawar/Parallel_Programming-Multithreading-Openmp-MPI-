#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"
#include "math.h"
int M,N,m,p,l,my_rank;
float *A,*B,*x;
double starttime,time1,time2;
MPI_Status status;
#define a(x,y) a[x*M+y]
#define b(x) b[x]
#define A(x,y) A[x*M+y]
#define B(x) B[x]
int main(int argc, char **argv)
{
    int i,j,t,k,my_rank,group_size,i1,i2,v,w,tem,row,col,*shift;
    float temp,*sum,*f,lmax,*a,*b;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&group_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if(argc==1)
    {
    	if(my_rank==0)
    		printf("Please specify the dimensions and make sure that the number is divisible by number of processes\n");
	      MPI_Finalize();
	      exit(1);
   }
    M = atoi(argv[1]);
    if(M%group_size!=0)
    {
    	if(my_rank==0)
    	    	printf("");
    		    MPI_Finalize();
    		    exit(1);
    }
    A=(float *)malloc(sizeof(float)*M*M);
    B=(float *)malloc(sizeof(float)*M);
    x=(float *)malloc(sizeof(float)*M);
    for (col = 0; col < M; col++)
        {
            for (row = 0; row < M; row++)
            {
                A[row*M+col] = (float)rand()/ 32768.0;
            }
            B[col] = (float)rand()/ 32768.0;
            x[col] = 0.0;
        }
    p=group_size;
    starttime=MPI_Wtime();
    MPI_Bcast(&M,1,MPI_INT,0,MPI_COMM_WORLD);
    m=M/p;
    if (M%p!=0)
    	m++;
    f=(float*)malloc(sizeof(float)*(M+1));
    a=(float*)malloc(sizeof(float)*m*M);
    b=(float*)malloc(sizeof(float)*m);
    sum=(float*)malloc(sizeof(float)*m);
    shift=(int*)malloc(sizeof(int)*M);
    if (a==NULL||b==NULL||f==NULL||sum==NULL||x==NULL||shift==NULL)
        printf("Memory not allocated\n");
    for(i=0;i<M;i++)
        shift[i]=i;


    if (my_rank==0)
    {
        for(i=0;i<m;i++)
            for(j=0;j<M;j++)
                a(i,j)=A(i*p,j);

        for(i=0;i<m;i++)
            b(i)=B(i*p);
    }

    if (my_rank==0)
    {
        for(i=0;i<M;i++)
            if ((i%p)!=0)
        {
            i1=i%p;
            i2=i/p+1;

            MPI_Send(&A(i,0),M,MPI_FLOAT,i1,i2,MPI_COMM_WORLD);
            MPI_Send(&B(i),1,MPI_FLOAT,i1,i2,MPI_COMM_WORLD);
        }
    }
    else
    {
        for(i=0;i<m;i++)
        {
            MPI_Recv(&a(i,0),M,MPI_FLOAT,0,i+1,MPI_COMM_WORLD,&status);
            MPI_Recv(&b(i),1,MPI_FLOAT,0,i+1,MPI_COMM_WORLD,&status);
        }
    }

    time1=MPI_Wtime();
    for(i=0;i<m;i++)
        for(j=0;j<p;j++)
    {
        if (my_rank==j)
        {
            v=i*p+j;
            lmax=a(i,v);
            l=v;

            for(k=v+1;k<M;k++)
                if (fabs(a(i,k))>lmax)
            {
                lmax=a(i,k);
                l=k;
            }

            if (l!=v)
            {
                for(t=0;t<m;t++)
                {
                    temp=a(t,v);
                    a(t,v)=a(t,l);
                    a(t,l)=temp;
                }

                tem=shift[v];
                shift[v]=shift[l];
                shift[l]=tem;
            }

            for(k=v+1;k<M;k++)
                a(i,k)=a(i,k)/a(i,v);

            b(i)=b(i)/a(i,v);
            a(i,v)=1;

            for(k=v+1;k<M;k++)
                f[k]=a(i,k);
            f[M]=b(i);


            MPI_Bcast(&f[0],M+1,MPI_FLOAT,my_rank,MPI_COMM_WORLD);

            MPI_Bcast(&l,1,MPI_INT,my_rank,MPI_COMM_WORLD);
        }
        else
        {
            v=i*p+j;
            MPI_Bcast(&f[0],M+1,MPI_FLOAT,j,MPI_COMM_WORLD);
            MPI_Bcast(&l,1,MPI_INT,j,MPI_COMM_WORLD);

            if (l!=v)
            {
                for(t=0;t<m;t++)
                {
                    temp=a(t,v);
                    a(t,v)=a(t,l);
                    a(t,l)=temp;
                }

                tem=shift[v];
                shift[v]=shift[l];
                shift[l]=tem;
            }
        }

        if (my_rank<=j)
            for(k=i+1;k<m;k++)
        {
            for(w=v+1;w<M;w++)
                a(k,w)=a(k,w)-f[w]*a(k,v);
            b(k)=b(k)-f[M]*a(k,v);
        }

        if (my_rank>j)
            for(k=i;k<m;k++)
        {
            for(w=v+1;w<M;w++)
                a(k,w)=a(k,w)-f[w]*a(k,v);
            b(k)=b(k)-f[M]*a(k,v);
        }
    }

    for(i=0;i<m;i++)
        sum[i]=0.0;

    for(i=m-1;i>=0;i--)
        for(j=p-1;j>=0;j--)
            if (my_rank==j)
            {
                x[i*p+j]=(b(i)-sum[i])/a(i,i*p+j);

                MPI_Bcast(&x[i*p+j],1,MPI_FLOAT,my_rank,MPI_COMM_WORLD);

                for(k=0;k<i;k++)
                    sum[k]=sum[k]+a(k,i*p+j)*x[i*p+j];
            }
            else
            {
        MPI_Bcast(&x[i*p+j],1,MPI_FLOAT,j,MPI_COMM_WORLD);

        if (my_rank>j)
            for(k=0;k<i;k++)
                sum[k]=sum[k]+a(k,i*p+j)*x[i*p+j];

        if (my_rank<j)
            for(k=0;k<=i;k++)
                sum[k]=sum[k]+a(k,i*p+j)*x[i*p+j];
    }

    if (my_rank!=0)
        for(i=0;i<m;i++)
            MPI_Send(&x[i*p+my_rank],1,MPI_FLOAT,0,i,MPI_COMM_WORLD);
    else
        for(i=1;i<p;i++)
            for(j=0;j<m;j++)
                MPI_Recv(&x[j*p+i],1,MPI_FLOAT,i,j,MPI_COMM_WORLD,&status);
    time2=MPI_Wtime();
    if (my_rank==0)
    {
        printf("\n");
        printf("Total time = %f seconds\n",time2-time1);
    }

    MPI_Finalize();
        free(a);
        free(b);
        free(x);
        free(f);
   return(0);
}
