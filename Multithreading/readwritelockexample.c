 #include <pthread.h>
#include<stdio.h>
 int minimum_value, partial_list_size;
 pthread_mutex_t minimum_value_lock;
typedef struct
{
       int readers;
       int writer;
       pthread_cond_t readers_proceed;
       pthread_cond_t writer_proceed;
       int pending_writers;
       pthread_mutex_t read_write_lock;

  } mylib_rwlock_t;
 void mylib_rwlock_init (mylib_rwlock_t *l)
{
   l -> readers = l -> writer = l -> pending_writers = 0;
    pthread_mutex_init(&(l -> read_write_lock), NULL);
   pthread_cond_init(&(l -> readers_proceed), NULL);
     pthread_cond_init(&(l -> writer_proceed), NULL);
  }
 void mylib_rwlock_rlock(mylib_rwlock_t *l)
{
     /* if there is a write lock or pending writers, perform condition wait.. else increment count of readers and grant read lock */
     pthread_mutex_lock(&(l -> read_write_lock));
     while ((l -> pending_writers > 0) || (l -> writer > 0))
         pthread_cond_wait(&(l -> readers_proceed), &(l -> read_write_lock));
    l -> readers ++;
    pthread_mutex_unlock(&(l -> read_write_lock));
  }
  void mylib_rwlock_wlock(mylib_rwlock_t *l)
{
  /* if there are readers or writers, increment pending writers count and wait. On being woken, decrement pending writers 34     count and increment writer count */
pthread_mutex_lock(&(l -> read_write_lock));
     while ((l -> writer > 0) || (l -> readers > 0)) {
       l -> pending_writers ++;
        pthread_cond_wait(&(l -> writer_proceed),&(l -> read_write_lock));
     }
     l -> pending_writers --;
     l -> writer ++;
     pthread_mutex_unlock(&(l -> read_write_lock));
  }
  void mylib_rwlock_unlock(mylib_rwlock_t *l)
 {
    /* if there is a write lock then unlock, else if there are
read locks, decrement count of read locks. If the count  is 0 and there is a pending writer, let it through, else if there are pending readers, let them all go through */
  pthread_mutex_lock(&(l -> read_write_lock));
     if (l -> writer > 0)
        l -> writer = 0;
     else if (l -> readers > 0)
        l -> readers --;
     pthread_mutex_unlock(&(l -> read_write_lock));
    if ((l -> readers == 0) && (l -> pending_writers > 0))
        pthread_cond_signal(&(l -> writer_proceed));
     else if (l -> readers > 0)
        pthread_cond_broadcast(&(l -> readers_proceed));
  }
  void *find_min_rw(void *list_ptr)
  {
         int *partial_list_pointer, my_min, i;
         mylib_rwlock_t read_write_lock;
         my_min = 0;
         partial_list_pointer = (int *) list_ptr;
         for (i = 0; i < partial_list_size; i++)
             if (partial_list_pointer[i] < my_min)
                 my_min = partial_list_pointer[i];
         /* lock the mutex associated with minimum_value and  update the variable as required */       mylib_rwlock_rlock(&read_write_lock);
        if (my_min < minimum_value)
  {
      mylib_rwlock_unlock(&read_write_lock);
            mylib_rwlock_wlock(&read_write_lock);
            minimum_value = my_min;
       }
        /* and unlock the mutex */
        mylib_rwlock_unlock(&read_write_lock);
        pthread_exit(0);
    }
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
 pthread_create( &thread1[i], NULL, find_min_rw, &list);
 for(i=0;i<n;i++)
  pthread_join( thread1[i], NULL);
  //pthread_join( thread2, NULL);
 printf("The minimum value out of 100 million integers is%d\n ",minimum_value);
 printf("NOTE: This program takes a list of size 100million integers as asked in the programming assighnment\n");
 printf("The number of threads running on this program is %d\n:",n);
 return 0;
   }
