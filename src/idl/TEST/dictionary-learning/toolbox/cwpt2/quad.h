/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#ifndef __QUAD_H_INCLUDED
#define __QUAD_H_INCLUDED

void quad_forward(FLOAT *signal, int x_size, int y_size, int scale,
    FLOAT * * trans_data,
    FLOAT * buffer, FLOAT *l, int loff, FLOAT *h, int hoff,
    int low_high_inv_x, int low_high_inv_y);

void quad_reverse(
    FLOAT * * trans_data,
    int x_size, int y_size, int scale,
    FLOAT *signal, FLOAT *buffer, FLOAT *l, int loff, FLOAT *h, int hoff,
    int low_high_inv_x, int low_high_inv_y);

void triway_forward(FLOAT *signal, int x_size, int y_size, int scale,
    FLOAT * * trans_data,
    FLOAT * buffer, FLOAT *l, int loff, FLOAT *h, int hoff,
    int low_high_inv_x, int low_high_inv_y);

void triway_reverse(
    FLOAT * * trans_data,
    int x_size, int y_size, int scale,
    FLOAT * signal, FLOAT * * buffer,
    FLOAT *la, int laoff, FLOAT *ha, int haoff,
    FLOAT *lr, int lroff, FLOAT *hr, int hroff,
    int low_high_inv_x, int low_high_inv_y);

#endif
    
