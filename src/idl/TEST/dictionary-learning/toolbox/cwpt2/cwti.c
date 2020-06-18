/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

/*=================================================================
 *
 * This is a MEX-file for Matlab, it serves as an interface between
 * the Matlab environment and external C functions.
 *
 * See function-specific documentation in the corresponding .m file.
 *
 * NOTE: Do not include comments starting by double slash in source files
 *       Functions in linked files have to be declared extern
 *
 * USEFUL Functions:
 *      mxIsDouble, mxIsClass, mxGetPr, mxGetScalar
 *      mxGetM, mxGetN, mxGetNumberOfDimensions, mxGetDimensions
 *      mxCreateDoubleMatrix, mxCreateNumericArray
 *      mexErrMsgTxt
 *
 *=================================================================*/

#include "mex.h"

/* INCLUDE HERE THE AUXILIARY FUNCTION PROTOTYPES */

#include "cwpt2.h"

/* IN THE FOLLOWING TRACK 'IN' AND 'OUT' IF CHANGING NUMBER OF INPUTS OR OUTPUTS */
/* Input Arguments */
#define NUM_IN              1
#define IN_1    prhs[0]

/* Output Arguments */
#define NUM_OUT             1
#define OUT_1   plhs[0]

/* This is the mexfunc part. Its prototype is NEVER changed */
void mexFunction(
    int nlhs, mxArray *plhs[],     /* left  hand side i.e. out */
    int nrhs, const mxArray*prhs[] /* right hand side i.e. in */
    )
{
    char mexErrMsg[150];
    int j;
    int sig_size;

    double * input_pointer;
    double * signal;
    int nb_scales;

    int dimensions[2];              /* For creating an output matrix */
    FLOAT * * subbands;

    int result;

    int out_runner;

    for (out_runner=0;out_runner<nlhs;out_runner++)
        plhs[out_runner] = NULL;

    /* Check for proper number of arguments */
    if ((nrhs != NUM_IN) || (nlhs != NUM_OUT))
        MEX_ERROR("Wrong number of input or output arguments.");

    /* Checking input argument types */
    if (!mxIsClass(IN_1, "double"))
        MEX_ERROR("The first argument must be 'real' and of type 'double'.");

    /* checking dimensions of input arguments */
    /* the order is okay, M = the first raster order */
    sig_size = mxGetM(IN_1);
    nb_scales = mxGetN(IN_1) - 1;

    if ((sig_size < 2) || (sig_size > 16384))
        MEX_ERROR("The first argument size is invalid.");

    if ((nb_scales < 1) || (nb_scales > 10))
        MEX_ERROR("Invalid 'nb_scales' value.");

    /* retrieving input arguments */
    input_pointer    = (double *) mxGetPr(IN_1);

    /* Although the input array is created as a matrix, we mask it by a double-referencing argument */

    subbands = (FLOAT * *) calloc(nb_scales+1, sizeof(FLOAT*));
    if (!subbands)
        MEX_ERROR("Memory allocation error.");

    subbands[0] = input_pointer;
    for (j=1; j < nb_scales+1; j++)
        subbands[j] = subbands[j-1] + sig_size;

    /* Create a matrix for the return argument : determining parameters ... */
    dimensions[0] = sig_size;
    dimensions[1] = 1;
    OUT_1 = mxCreateNumericArray(2, dimensions, mxDOUBLE_CLASS, mxREAL);

    /* Assign pointers to the output parameters */
    signal = (double *) mxGetPr(OUT_1);

    /* Do the actual computations in a subroutine */
    result = ReverseDyadicTransform(subbands, sig_size, nb_scales, signal, NULL, NULL, 0, NULL, 0);

    free(subbands);

    if (result)
        MEX_ERROR("Error in CWPT2 module.");

    return;
mex_err_quit:
    for (out_runner=0;out_runner<nlhs;out_runner++)
    {
        if (plhs[out_runner])
            mxDestroyArray(plhs[out_runner]);
        plhs[out_runner] = mxCreateString(mexErrMsg);
    }
/*    mexErrMsgTxt(); */
}
