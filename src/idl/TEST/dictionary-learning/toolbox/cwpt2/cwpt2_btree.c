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

/* #define  mymxGetPr mxGetPr */
#define  mymxGetPr mxGetData

/* This is the mexfunc part. Its prototype is NEVER changed */
void mexFunction(
    int nlhs, mxArray *plhs[],     /* left  hand side i.e. out */
    int nrhs, const mxArray*prhs[] /* right hand side i.e. in */
    )
{
    char mexErrMsg[150];
    int btree_type;
    int nb_scales;
    int dimensions[2];              /* For creating an output matrix */
    char * output_pointer;
    int out_runner;

    for (out_runner=0;out_runner<nlhs;out_runner++)
        plhs[out_runner] = NULL;

    /* Check for proper number of arguments */
    if ((nrhs != NUM_IN) || (nlhs != NUM_OUT))
        MEX_ERROR("Wrong number of input or output arguments.")

    /* Checking input argument types */
    if (!mxIsClass(IN_1, "double"))
        MEX_ERROR("The first argument must be 'real' and of type 'double'.")
    if (!mxIsClass(IN_2, "double"))
        MEX_ERROR("The second argument must be 'real' and of type 'double'.")

    /* Checking dimensions of input arguments */
    if ( mxGetM(IN_1) * mxGetN(IN_1) != 1 )
        MEX_ERROR("The first argument must be a scalar.")
    if ( mxGetM(IN_2) * mxGetN(IN_2) != 1 )
        MEX_ERROR("The second argument must be a scalar.")

    /* retrieving input arguments */
    nb_scales     = (int) mxGetScalar(IN_1);
    btree_type    = (int) mxGetScalar(IN_2);

    if ((nb_scales < 1) || (nb_scales > 10))
        MEX_ERROR("Invalid 'nb_scales' parameter.")

    /* Create a matrix for the return argument : determining parameters ... */
    dimensions[0] = btree_get_size(nb_scales);
    dimensions[1] = 1;
    OUT_1 = mxCreateNumericArray(2, dimensions, mxUINT8_CLASS, mxREAL);

    /* Assign pointers to the output parameters */
    output_pointer = (char *) mymxGetPr(OUT_1);

    switch (btree_type)
    {
    case 1:
        btree_fill_full(output_pointer, nb_scales);
        break;
    case 2:
        btree_fill_wavelet(output_pointer, nb_scales);
        break;
    default:
        MEX_ERROR("Unknown btree type. Please, choose 1 (full) or 2 (wavelet).")
    }

    return;
mex_err_quit:
    for (out_runner=0;out_runner<nlhs;out_runner++)
    {
        if (plhs[out_runner])
            mxDestroyArray(plhs[out_runner]);
        plhs[out_runner] = mxCreateString(mexErrMsg);
    }
    /* mexErrMsgTxt(mexErrMsg);*/
}
