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
#define NUM_IN              2       /* or 6 */
#define IN_1    prhs[0]
#define IN_2    prhs[1]
#define IN_3    prhs[2]
#define IN_4    prhs[3]
#define IN_5    prhs[4]
#define IN_6    prhs[5]

/* Output Arguments */
#define NUM_OUT             1
#define OUT_1   plhs[0]

#define  mymxGetPr mxGetData

/* This is the mexfunc part. Its prototype is NEVER changed */
void mexFunction(
    int nlhs, mxArray *plhs[],     /* left  hand side i.e. out */
    int nrhs, const mxArray*prhs[] /* right hand side i.e. in */
    )
{
    char mexErrMsg[150];
    int x_size, y_size;
    int btree_size;

    double * signal;
    char * btree;
    double * l, * h;
    int loff, hoff;

    int dimensions[3];              /* For creating an output matrix */
    double * output_pointer;

    int result;
    int out_runner;

    for (out_runner=0;out_runner<nlhs;out_runner++)
        plhs[out_runner] = NULL;

    /* Check for proper number of arguments */
    if (((nrhs != NUM_IN) && (nrhs != 6)) || (nlhs != NUM_OUT))
        MEX_ERROR("Wrong number of input or output arguments.");

    /* Checking input argument types */
    if (!mxIsClass(IN_1, "double"))
        MEX_ERROR("The first argument must be 'real' and of type 'double'.");
    if (!mxIsClass(IN_2, "uint8"))
        MEX_ERROR("The first argument must be of type 'uint8'.");
    if (nrhs == 6)
    {
        if (!mxIsClass(IN_3, "double"))
            MEX_ERROR("The argument 3 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_4, "double"))
            MEX_ERROR("The argument 4 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_5, "double"))
            MEX_ERROR("The argument 5 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_6, "double"))
            MEX_ERROR("The argument 6 must be 'real' and of type 'double'.");
    }

    /* checking dimensions of input arguments */
    /* the order is okay, M = the first raster order */
    x_size = mxGetM(IN_1);
    y_size = mxGetN(IN_1);

    if ((x_size < 2) || (y_size < 2) || (x_size > 16384) || (y_size > 16384))
        MEX_ERROR("The first argument size is invalid.");

    if ((mxGetM(IN_2) != 1) && (mxGetN(IN_2) != 1))
        MEX_ERROR("The second argument must be a vector.");

    btree_size = mxGetM(IN_2) * mxGetN(IN_2);

    if (nrhs == 6)
    {
        if ((mxGetM(IN_3) != 1) && (mxGetN(IN_3) != 1))
            MEX_ERROR("The argument 3 must be a vector.");
        if ((mxGetM(IN_4) != 1) || (mxGetN(IN_4) != 1))
            MEX_ERROR("The argument 4 must be a scalar.");
        if ((mxGetM(IN_5) != 1) && (mxGetN(IN_5) != 1))
            MEX_ERROR("The argument 5 must be a vector.");
        if ((mxGetM(IN_6) != 1) || (mxGetN(IN_6) != 1))
            MEX_ERROR("The argument 6 must be a scalar.");
    }

    /* retrieving input arguments */
    signal    = (double *) mymxGetPr(IN_1);
    btree     = (char *) mymxGetPr(IN_2);
    if (nrhs == 6)
    {
        l         = (double *) mymxGetPr(IN_3);
        loff      = (int) mxGetScalar(IN_4);
        h         = (double *) mymxGetPr(IN_5);
        hoff      = (int) mxGetScalar(IN_6);
    }
    else
    {
        l = h = NULL;
        loff = hoff = 0;
    }

    /* checking read arguments */
    if (!btree_is_valid(btree,btree_size))
        MEX_ERROR("Invalid btree.");
    if (nrhs == 6)
    {
        if (mxGetM(IN_3) * mxGetN(IN_3) != (2 * loff + 1))
            MEX_ERROR("Invalid low-pass filter offset.");
        if (mxGetM(IN_5) * mxGetN(IN_5) != (2 * hoff + 1))
            MEX_ERROR("Invalid high-pass filter offset.");
    }

    /* Create a matrix for the return argument : determining parameters ... */
    dimensions[0] = x_size;
    dimensions[1] = y_size;
    dimensions[2] = btree_get_nb_leaves(btree, btree_size, 3);
    OUT_1 = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS, mxREAL);

    /* Assign pointers to the output parameters */
    output_pointer = (double *) mymxGetPr(OUT_1);

    /* Do the actual computations in a subroutine */
    result = forward3(signal, x_size, y_size, btree, btree_size, output_pointer, l, loff, h, hoff);

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
