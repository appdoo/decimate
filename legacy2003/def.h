#ifndef DEF_H
#define DEF_H

#include <float.h>
#include <math.h>

// Definition of one and zero
#define IsDZero(a) (fabs(a)       < 10.0 * DBL_MIN)
#define IsDOne(a)  (fabs(a - 1.0) < 10.0 * DBL_EPSILON)

#define IsFZero(a) (fabs(a)       < 10.0 * FLT_MIN)
#define IsFOne(a)  (fabs(a - 1.0) < 10.0 * FLT_EPSILON)

#endif