// nrutil.hpp
//
// Original author:  Yurii Il'inskii
//
// Modifications:
// 04 Dec 1996: Modified by Ronald Kumon to remove obsolete #pragma
//              construction; also realigned text and added comments

#ifndef NRUTIL_H
#define NRUTIL_H

extern void nrerror(char *error_text);
extern float *vector(long nl, long nh);
extern double *dvector(long nl, long nh);
extern double **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_vector(float *v, long nl, long nh);
extern void free_dvector(double *v, long nl, long nh);
extern void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

#endif //NRUTIL_H
