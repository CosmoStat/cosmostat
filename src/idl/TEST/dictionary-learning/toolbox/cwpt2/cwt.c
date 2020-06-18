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
#define NUM_IN              2
#define IN_1    prhs[0]
#define IN_2    prhs[1]

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
    int x_size, y_size, sig_size;
    int j;
    int nb_scales;

    double * signal;

    int dimensions[2];              /* For creating an output matrix */
    double * output_pointer;
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

    if (!mxIsClass(IN_2, "double"))
        MEX_ERROR("The second argument must be 'real' and of type 'double'.");

    /* checking dimensions of input arguments */
    /* the order is okay, M = the first raster order */
    x_size = mxGetM(IN_1);
    y_size = mxGetN(IN_1);

    if ((x_size != 1) && (y_size !=1))
        MEX_ERROR("The first argument must be a vector.");

    sig_size = x_size * y_size;

    x_size = mxGetM(IN_2);
    y_size = mxGetN(IN_2);

    if ((x_size != 1) || (y_size !=1))
        MEX_ERROR("The second argument must be a scalar.");

    /* retrieving input arguments */
    signal    = (double *) mxGetPr(IN_1);
    nb_scales = (int) mxGetScalar(IN_2);

    if ((nb_scales < 1) || (nb_scales > 10))
        MEX_ERROR("Invalid 'nb_scales' value.");

    /* Create a matrix for the return argument : determining parameters ... */
    dimensions[0] = sig_size;
    dimensions[1] = nb_scales + 1;
    OUT_1 = mxCreateNumericArray(2, dimensions, mxDOUBLE_CLASS, mxREAL);

    /* Assign pointers to the output parameters */
    output_pointer = (double *) mxGetPr(OUT_1);

    /* Although the output array is created as a matrix, we mask it by a double-referencing argument */

    subbands = (FLOAT * *) calloc(nb_scales+1, sizeof(FLOAT*));
    if (!subbands)
        MEX_ERROR("Memory allocation error.");

    subbands[0] = output_pointer;
    for (j=1; j < nb_scales+1; j++)
        subbands[j] = subbands[j-1] + sig_size;

    /* Do the actual computations in a subroutine */
    result = ForwardDyadicTransform(signal, sig_size, nb_scales, subbands, NULL, 0, NULL, 0);

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
