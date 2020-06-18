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
    char * btree;
    int btree_size;
    int nb_leaves;

    int param_x_size, param_y_size;
    int children_per_parent;

    int result;
    int dimensions[2];              /* For creating an output matrix */
    int * output_pointer;

    int out_runner;

    for (out_runner=0;out_runner<nlhs;out_runner++)
        plhs[out_runner] = NULL;

    /* Check for proper number of arguments */
    if ((nrhs != NUM_IN) || (nlhs != NUM_OUT))
        MEX_ERROR("Wrong number of input or output arguments.");

    /* Checking input argument types */
    if (!mxIsClass(IN_1, "uint8"))
        MEX_ERROR("The first argument must be of type 'uint8'.");
    if (!mxIsClass(IN_1, "double"))
        MEX_ERROR("The second argument must be of type 'double'.");

    /* Checking dimensions of input arguments */
    param_x_size = mxGetM(IN_1);
    param_y_size = mxGetN(IN_1);
    if ((param_x_size != 1) && (param_y_size != 1))
        MEX_ERROR("The second argument must be a vector.");

    btree_size = param_x_size * param_y_size;

    if ((mxGetM(IN_2) != 1) || (mxGetN(IN_2) != 1))
        MEX_ERROR("The second argument must be a scalar.");

    /* retrieving input arguments */
    btree     = (char *) mxGetPr(IN_1);

    if (!btree_is_valid(btree,btree_size))
        MEX_ERROR("Invalid btree.");

    children_per_parent = (int) mxGetScalar(IN_2);
    if ((children_per_parent != 3) && (children_per_parent !=4))
        MEX_ERROR("Only quad- or tri-trees are supported.");

    /* Create a matrix for the return argument : determining parameters ... */
    nb_leaves = btree_get_nb_leaves(btree, btree_size, children_per_parent);

    dimensions[0] = nb_leaves;
    dimensions[1] = 1;
    OUT_1 = mxCreateNumericArray(2, dimensions, mxINT32_CLASS, mxREAL);

    /* Assign pointers to the output parameters */
    output_pointer = (int *) mxGetPr(OUT_1);

    result = btree_get_leaf_order(btree, btree_size, output_pointer, nb_leaves);

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
