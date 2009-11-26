#ifndef EIGENFUNCTIONS_H
#define EIGENFUNCTIONS_H
#include <complex>
#include "./prob_engine.h"
/* Functions */
using std::complex;
inline void zhetrd3(complex<double> A[3][3], complex<double> Q[3][3],
                    double d[3], double e[2]);
int zheevc3(complex<double> A[3][3], double w[3]);
int zheevq3(complex<double> A[3][3], complex<double> Q[3][3], double w[3]);
int zheevh3(complex<double> A[3][3], complex<double> Q[3][3], double w[3]);

#endif /*EIGENFUNCTIONS_H*/
