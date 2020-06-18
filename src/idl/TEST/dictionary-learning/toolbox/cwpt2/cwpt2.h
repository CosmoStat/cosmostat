/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#ifndef __CWPT2_H_INCLUDED
#define __CWPT2_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>

typedef double FLOAT;

#define INLINE              /* inline is not supported by Matlab compiler */
/* #define INLINE inline */

#define CWPT2_SCALE_MAX         10

#define MEX_ERROR(a) {sprintf(mexErrMsg,a);goto mex_err_quit;}

#include "atrous.h"
#include "btree.h"
#include "dyadic2.h"
#include "quad.h"

#endif
