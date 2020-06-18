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


#define  mymxGetPr mxGetData

/* IN THE FOLLOWING TRACK 'IN' AND 'OUT' IF CHANGING NUMBER OF INPUTS OR OUTPUTS */
/* Input Arguments */
#define NUM_IN              2       /* or 10 */
#define IN_1    prhs[0]
#define IN_2    prhs[1]
#define IN_3    prhs[2]
#define IN_4    prhs[3]
#define IN_5    prhs[4]
#define IN_6    prhs[5]
#define IN_7    prhs[6]
#define IN_8    prhs[7]
#define IN_9    prhs[8]
#define IN_10    prhs[9]

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
    int btree_size;
    int x_size, y_size;
    int first_scale_shifted;

    double * input_pointer;
    const int * input_dimensions;
    double * signal;
    char * btree;
    double * la, * ha;
    int laoff, haoff;
    double * lr, * hr;
    int lroff, hroff;

    int dimensions[2];              /* For creating an output matrix */

    int result;

    int out_runner;

    for (out_runner=0;out_runner<nlhs;out_runner++)
        plhs[out_runner] = NULL;

    /* Check for proper number of arguments */
    if (((nrhs != NUM_IN) && (nrhs != 10)) || (nlhs != NUM_OUT))
        MEX_ERROR("Wrong number of input or output arguments.");

    /* Checking input argument types */
    if (!mxIsClass(IN_1, "double"))
        MEX_ERROR("The first argument must be 'real' and of type 'double'.");
    if (!mxIsClass(IN_2, "uint8"))
        MEX_ERROR("The first argument must be of type 'uint8'.");
    if (nrhs == 10)
    {
        if (!mxIsClass(IN_3, "double"))
            MEX_ERROR("The argument 3 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_4, "double"))
            MEX_ERROR("The argument 4 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_5, "double"))
            MEX_ERROR("The argument 5 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_6, "double"))
            MEX_ERROR("The argument 6 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_7, "double"))
            MEX_ERROR("The argument 7 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_8, "double"))
            MEX_ERROR("The argument 8 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_9, "double"))
            MEX_ERROR("The argument 9 must be 'real' and of type 'double'.");
        if (!mxIsClass(IN_10, "double"))
            MEX_ERROR("The argument 10 must be 'real' and of type 'double'.");
    }

    if (mxGetNumberOfDimensions(IN_1)!=3)
        MEX_ERROR("The first argument has invalid dimension list.");

    /* checking dimensions of input arguments */
    /* the order is okay, M = the first raster order */
    input_dimensions = mxGetDimensions(IN_1);
    x_size = input_dimensions[0];
    y_size = input_dimensions[1];

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
        if ((mxGetM(IN_7) != 1) && (mxGetN(IN_7) != 1))
            MEX_ERROR("The argument 7 must be a vector.");
        if ((mxGetM(IN_8) != 1) || (mxGetN(IN_8) != 1))
            MEX_ERROR("The argument 8 must be a scalar.");
        if ((mxGetM(IN_9) != 1) && (mxGetN(IN_9) != 1))
            MEX_ERROR("The argument 9 must be a vector.");
        if ((mxGetM(IN_10) != 1) || (mxGetN(IN_10) != 1))
            MEX_ERROR("The argument 10 must be a scalar.");
    }

    /* retrieving input arguments */
    input_pointer    = (double *) mymxGetPr(IN_1);
    btree = (char *) mymxGetPr(IN_2);
    if (nrhs == 10)
    {
        la         = (double *) mymxGetPr(IN_3);
        laoff      = (int) mxGetScalar(IN_4);
        ha         = (double *) mymxGetPr(IN_5);
        haoff      = (int) mxGetScalar(IN_6);
        lr         = (double *) mymxGetPr(IN_7);
        lroff      = (int) mxGetScalar(IN_8);
        hr         = (double *) mymxGetPr(IN_9);
        hroff      = (int) mxGetScalar(IN_10);
    }
    else
    {
        la = ha = lr = hr = NULL;
        laoff = haoff = lroff = hroff = 0;
    }

    /* checking read arguments */
    if (!btree_is_valid(btree,btree_size))
        MEX_ERROR("Invalid btree.");
    if (nrhs == 10)
    {
        if (mxGetM(IN_3) * mxGetN(IN_3) != (2 * laoff + 1))
            MEX_ERROR("Invalid low-pass analysis filter offset.");
        if (mxGetM(IN_5) * mxGetN(IN_5) != (2 * haoff + 1))
            MEX_ERROR("Invalid high-pass analysis filter offset.");
        if (mxGetM(IN_7) * mxGetN(IN_7) != (2 * lroff + 1))
            MEX_ERROR("Invalid low-pass reconstruction filter offset.");
        if (mxGetM(IN_9) * mxGetN(IN_9) != (2 * hroff + 1))
            MEX_ERROR("Invalid high-pass reconstruction filter offset.");
    }

    first_scale_shifted = 0;
    if (input_dimensions[2] != btree_get_nb_leaves(btree, btree_size, 3))
    {
        if (input_dimensions[2] == btree_get_nb_leaves(btree, btree_size, 3) - 2)
            first_scale_shifted = 1;
        else
            MEX_ERROR("Invalid size of the first argument: does not correspond to btree.");
    }

    /* Create a matrix for the return argument : determining parameters ... */
    dimensions[0] = x_size;
    dimensions[1] = y_size;
    OUT_1 = mxCreateNumericArray(2, dimensions, mxDOUBLE_CLASS, mxREAL);

    /* Assign pointers to the output parameters */
    signal = (double *) mymxGetPr(OUT_1);

    /* Do the actual computations in a subroutine */
    result = reverse3(input_pointer, x_size, y_size, btree, btree_size, signal,
        la, laoff, ha, haoff, lr, lroff, hr, hroff, first_scale_shifted);

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
