#include "TopLJets2015/TopAnalysis/interface/timeval_private.hh"

timeval_private::timeval_private()
{
  gettimeofday(this, 0);
}

timeval_private timeval_private::operator-(timeval_private& other)
{
  timeval_private ret;
  ret.tv_sec = tv_sec - other.tv_sec;
  ret.tv_usec = tv_usec - other.tv_usec;
  return ret;
}
double timeval_private::return_msec()
{
  return tv_sec*1000 + tv_usec/1E3;
}

double timeval_private::return_sec()
{
  return tv_sec + tv_usec/1E6;
}
void timeval_private::print_sec()
{
  printf("seconds %f\n", return_sec());
}

void timeval_private::print()
{
  printf("sec %lu, usec %lu\n", tv_sec, tv_usec);
}

