/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#ifndef __ATROUS_H_INCLUDED
#define __ATROUS_H_INCLUDED

INLINE void atrous (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale);

INLINE void add_atrous (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale);

INLINE void atrous_stride (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale, int stride);

INLINE void add_atrous_stride (FLOAT *dst, FLOAT *src, int n,
    FLOAT *filter, int filter_offset, int scale, int stride);

#endif

