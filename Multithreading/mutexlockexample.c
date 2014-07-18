 #include <pthread.h>
#include<stdio.h>  
 void *find_min(void *list_ptr);
 pthread_mutex_t minimum_value_lock;
   int minimum_value, partial_list_size;
 int  main() 
{
    /* declare and initialize data structures and list */
static long  list[100000000],i,n;
 minimum_value = 0;
  printf("enter the number of threads for which you want to run this program \n");
  scanf("%d",&n);
  pthread_t thread1[n];
 pthread_mutex_init(&minimum_value_lock, NULL); 
  /* initialize lists, list_ptr, and partial_list_size */ 
for(i=0;i<100;i++)
{
list[i]=random(100000000);
//printf("the value is %d",list[i]);
}

partial_list_size =100000000/n;
      /* create and join threads here */ 
for(i=0;i<n;i++)
pthread_create( &thread1[i], NULL, find_min, &list);
for(i=0;i<n;i++)
 pthread_join( thread1[i], NULL);
 //pthread_join( thread2, NULL);
printf("The minimum value out of 100 million integers is%d\n ",minimum_value);
printf("NOTE: This program takes a list of size 100million integers as asked in the programming assighnment\n");
printf("The number of threads running on this program is %d\n:",n);
return 0;
  } 
  void *find_min(void *list_ptr) 
{
     int *partial_list_pointer, my_min, i; 
      my_min = 0;
      partial_list_pointer = (int *) list_ptr; 
      for (i = 0; i < partial_list_size; i++) 
          if (partial_list_pointer[i] < my_min) 
         my_min = partial_list_pointer[i];   
   /* lock the mutex associated with minimum_value and update the variable as required */
      pthread_mutex_lock(&minimum_value_lock);
      if (my_min < minimum_value)     
      minimum_value = my_min;
      /* and unlock the mutex */
      pthread_mutex_unlock(&minimum_value_lock);
  pthread_exit(0);  
 }
