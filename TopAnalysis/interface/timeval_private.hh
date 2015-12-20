#ifndef _timeval_private_hh_
#define _timeval_private_hh_
#include <stdio.h>
#include <sys/time.h>


struct timeval_private:public timeval 
{
  //time_t      tv_sec;     /* seconds */
  //suseconds_t tv_usec;    /* microseconds */
  timeval_private();
  timeval_private operator-(timeval_private& other);
  double return_msec();
  double return_sec();
  void print_sec();
  
  void print();
};

#endif
