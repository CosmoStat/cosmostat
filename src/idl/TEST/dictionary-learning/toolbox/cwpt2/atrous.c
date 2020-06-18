/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#include "cwpt2.h"

/*--------------------------------------------------------------*/
/* Symmetric extension convolution with 'à trous' technique.    */
/*--------------------------------------------------------------*/

INLINE void atrous (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale)
{
    int i,k;
    int fos= filter_offset*scale;

    /* Left boundary */

    for (i=0; i<n && i<fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset= i+scale*k;
            while (offset<0) offset+= 2*n-2;
            while (offset>=2*n-2) offset-= 2*n-2;
            if (offset>=n) offset= 2*n-2-offset;
            sum+= src[offset]*filter[k+filter_offset];
        }
        dst[i]= sum;
    }

    /* Middle */

    for (; i<n-fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
            sum+= src[i+scale*k]*filter[k+filter_offset];
        dst[i]= sum;
    }

    /* Right boundary */

    for (; i<n; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset= i+scale*k;
            while (offset<0) offset+= 2*n-2;
            while (offset>=2*n-2) offset-= 2*n-2;
            if (offset>=n) offset= 2*n-2-offset;
            sum+= src[offset]*filter[k+filter_offset];
        }
        dst[i]= sum;
    }
}

/*--------------------------------------------------------------*/
/* Same as before, but adds the result to dst[] instead of      */
/* setting dst[] to be the result.                              */
/*--------------------------------------------------------------*/

INLINE void add_atrous (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale)
{
    int i,k;
    int fos= filter_offset*scale;

    /* Left boundary */

    for (i=0; i<n && i<fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset= i+scale*k;
            while (offset<0) offset+= 2*n-2;
            while (offset>=2*n-2) offset-= 2*n-2;
            if (offset>=n) offset= 2*n-2-offset;
            sum+= src[offset]*filter[k+filter_offset];
         }
        dst[i]+= sum;
    }

    /* Middle */

    for (; i<n-fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
            sum+= src[i+scale*k]*filter[k+filter_offset];
        dst[i]+= sum;
    }

    /* Right boundary */

    for (; i<n; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset= i+scale*k;
            while (offset<0) offset+= 2*n-2;
            while (offset>=2*n-2) offset-= 2*n-2;
            if (offset>=n) offset= 2*n-2-offset;
            sum+= src[offset]*filter[k+filter_offset];
        }
        dst[i]+= sum;
    }
}

/*--------------------------------------------------------------*/
/* Two functions same as above, but using a stride parameter    */
/* for vertical filtering                                       */
/*--------------------------------------------------------------*/

INLINE void atrous_stride (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale, int stride)
{
    int i,k;
    int fos= filter_offset*scale;

    /* Left boundary */

    for (i=0; i<n && i<fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset = i+scale*k;
            while (offset < 0)
                offset += 2*n-2;
            while (offset >= 2*n-2)
                offset -= 2*n-2;
            if (offset >= n)
                offset = 2*n-2-offset;

            offset *= stride;       /* TODO: OPTIMIZE? */
            sum+= src[offset]*filter[k+filter_offset];
        }
        dst[i*stride]= sum; /* TODO: OPTIMIZE? */
    }

    /* Middle */

    for (; i<n-fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            /* TODO: OPTIMIZE? */
            sum += src[(i+scale*k)*stride] * filter[k+filter_offset];
        }
        dst[i*stride]= sum; /* TODO: OPTIMIZE? */
    }

    /* Right boundary */

    for (; i<n; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset= i+scale*k;
            while (offset<0) offset+= 2*n-2;
            while (offset>=2*n-2) offset-= 2*n-2;
            if (offset>=n) offset= 2*n-2-offset;

            offset *= stride;       /* TODO: OPTIMIZE? */
            sum+= src[offset]*filter[k+filter_offset];
        }
        dst[i*stride]= sum; /* TODO: OPTIMIZE? */
    }
}

INLINE void add_atrous_stride (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale, int stride)
{
    int i,k;
    int fos= filter_offset*scale;

    /* Left boundary */

    for (i=0; i<n && i<fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset= i+scale*k;
            while (offset<0) offset+= 2*n-2;
            while (offset>=2*n-2) offset-= 2*n-2;
            if (offset>=n) offset= 2*n-2-offset;
            sum+= src[offset*stride]*filter[k+filter_offset];
         }
        dst[i*stride]+= sum;
    }

    /* Middle */

    for (; i<n-fos; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
            sum+= src[(i+scale*k)*stride]*filter[k+filter_offset];
        dst[i*stride]+= sum;
    }

    /* Right boundary */

    for (; i<n; i++)
    {
        FLOAT sum= 0;
        for (k=-filter_offset; k<=filter_offset; k++)
        {
            int offset= i+scale*k;
            while (offset<0) offset+= 2*n-2;
            while (offset>=2*n-2) offset-= 2*n-2;
            if (offset>=n) offset= 2*n-2-offset;
            sum+= src[offset*stride]*filter[k+filter_offset];
        }
        dst[i*stride]+= sum;
    }
}
