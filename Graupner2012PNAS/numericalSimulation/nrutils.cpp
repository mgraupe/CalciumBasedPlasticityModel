// Nrutil.cc
//
// Original author:  Numerical Recipes, Yurii Il'inskii
//
// Modifications:
// 04 Dec 1996: Modified by Ronald Kumon to remove obsolete #pragma
//              construction; also realigned text and added comments;

#ifndef NRUTIL_CC
#define NRUTIL_CC
#define NR_END 1
#define FREE_ARG char*

#include "nrutils.hpp"
#include <iostream>
#include <stdlib.h>
#include <sstream>
//#include <stdio.h>
//#include <stddef.h>

using namespace std;

void nrerror(char *error_text)
     // Numerical recipes standard error handler
{
  cerr << "Numerical Recipes run-time error...\n";
  cerr << error_text<<"\n";
  cerr << "...now exiting to system...\n";
  exit(1);
}

float *vector(long nl, long nh)
// Allocate a float vector with subscript range v[nl..nh]
{
  float *v;

  v=(float *)malloc((size_t) (nh-nl+1+NR_END)*sizeof(float));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

double *dvector(long nl, long nh)
     // Allocates a vector of type double with range [nl..nh] 
{
  double *v;

  v=(double *)malloc((size_t)(nh-nl+1+NR_END)*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl+NR_END;  
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
     // Allocates a matrix of type double with range [nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

 // Allocate pointers to rows
  m=(double**) malloc((size_t)((nrow+NR_END)*sizeof(double)));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;
	
  // Allocate rows and set pointers to them
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  // Return pointer to array of pointers to rows
  return m;
}

void free_vector(float *v, long nl, long nh)
  // Free a float vector allocated with vector()
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
  // Free a vector of type double allocated by dvector()
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
  // Free a matrix of type double allocated by dmatrix()
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

#endif //NRUTIL_CC 
