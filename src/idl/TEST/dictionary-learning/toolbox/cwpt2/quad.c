/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#include "cwpt2.h"

/*--------------------------------------------------------------*/
/* 2D Forward Dyadic Wavelet Transform at 1 scale (quad)        */
/*--------------------------------------------------------------*/

/*
   subband
    0  1    LL HL       (for instance, 1 = HL = high_X, low_Y)
    2  3    LH HH

    1) X low-filtering              of the signal   into the buffer
    2) Y low and high filtering     of the buffer   into subbands 0 and 2
    3) X high-filtering             of the signal   into the buffer
    4) Y low and high filtering     of the buffer   into subbands 1 and 3

    All filters are of the same scale.
*/

void quad_forward(FLOAT *signal, int x_size, int y_size, int scale,
    FLOAT * * trans_data,
    FLOAT * buffer, FLOAT *l, int loff, FLOAT *h, int hoff,
    int low_high_inv_x, int low_high_inv_y)
{
    int offset;
    FLOAT * signal_off;

    FLOAT * trans_data_0_off, * trans_data_1_off, * trans_data_2_off, * trans_data_3_off;
    FLOAT * buffer_off;

    FLOAT * l_x, * h_x, * l_y, * h_y;
    int loff_x, hoff_x, loff_y, hoff_y;

    CWPT2_LOW_HIGH_INV

    /* analysis */

    /* step 1: X low-pass signal -> buffer */
    buffer_off = buffer;
    signal_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(buffer_off, signal_off, x_size, l_x, loff_x, scale);
        buffer_off += x_size;
        signal_off += x_size;
    }

    /* step 2: Y low-pass -> subband 0; Y high-pass -> subband 2 */
    trans_data_0_off = trans_data[0];
    trans_data_2_off = trans_data[2];
    buffer_off = buffer;
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(trans_data_0_off,   buffer_off, y_size, l_y, loff_y, scale, x_size);
        atrous_stride(trans_data_2_off, buffer_off, y_size, h_y, hoff_y, scale, x_size);
        trans_data_0_off ++;
        trans_data_2_off ++;
        buffer_off ++;
    }

    /* step 3: X high-pass signal -> buffer */
    buffer_off = buffer;
    signal_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(buffer_off, signal_off, x_size, h_x, hoff_x, scale);
        buffer_off += x_size;
        signal_off += x_size;
    }

    /* step 4: Y low-pass -> subband 1; Y high-pass -> subband 3 */
    trans_data_1_off = trans_data[1];
    trans_data_3_off = trans_data[3];
    buffer_off = buffer;
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(trans_data_1_off, buffer_off, y_size, l_y, loff_y, scale, x_size);
        atrous_stride(trans_data_3_off, buffer_off, y_size, h_y, hoff_y, scale, x_size);
        trans_data_1_off ++;
        trans_data_3_off ++;
        buffer_off ++;
    }

}

/*--------------------------------------------------------------*/
/* 2D Reverse Dyadic Wavelet Transform at 1 scale (quad)        */
/*--------------------------------------------------------------*/

/*
   subband
    0  1    LL HL       (for instance, 1 = HL = high_X, low_Y)
    2  3    LH HH

    1) Y low and high filtering     of subbands 1 and 3     into the buffer
    2) X high-filtering             of the buffer           into the signal
    3) Y low and high filtering     of subbands 0 and 2     into the buffer
    4) X low-filtering              of the buffer           into the signal

    All filters are of the same scale.
*/

void quad_reverse(
    FLOAT * * trans_data,
    int x_size, int y_size, int scale,
    FLOAT *signal, FLOAT *buffer, FLOAT *l, int loff, FLOAT *h, int hoff,
    int low_high_inv_x, int low_high_inv_y)
{
    int offset;
    FLOAT * signal_off;

    FLOAT * trans_data_0_off, * trans_data_1_off, * trans_data_2_off, * trans_data_3_off;
    FLOAT * buffer_off;

    FLOAT * l_x, * h_x, * l_y, * h_y;
    int loff_x, hoff_x, loff_y, hoff_y;

    CWPT2_LOW_HIGH_INV

    /* reconstruction */

    /* step 1: Y low-pass subband 1 -> buffer; Y high-pass subband 3 -> buffer*/
    trans_data_1_off = trans_data[1];
    trans_data_3_off = trans_data[3];
    buffer_off = buffer;
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(buffer_off, trans_data_1_off, y_size, l_y, loff_y, scale, x_size);
        add_atrous_stride(buffer_off, trans_data_3_off, y_size, h_y, hoff_y, scale, x_size);
        trans_data_1_off ++;
        trans_data_3_off ++;
        buffer_off ++;
    }

    /* step 2: X high-pass buffer -> signal */
    buffer_off = buffer;
    signal_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(signal_off, buffer_off, x_size, h_x, hoff_x, scale);
        buffer_off += x_size;
        signal_off += x_size;
    }

    /* step 3: Y low-pass subband 0 -> buffer; Y high-pass subband 2 -> buffer*/
    trans_data_0_off = trans_data[0];
    trans_data_2_off = trans_data[2];
    buffer_off = buffer;
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(buffer_off, trans_data_0_off, y_size, l_y, loff_y, scale, x_size);
        add_atrous_stride(buffer_off, trans_data_2_off, y_size, h_y, hoff_y, scale, x_size);
        trans_data_0_off ++;
        trans_data_2_off ++;
        buffer_off ++;
    }

    /* step 4: X low-pass buffer -> signal */
    buffer_off = buffer;
    signal_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        add_atrous(signal_off, buffer_off, x_size, l_x, loff_x, scale);
        buffer_off += x_size;
        signal_off += x_size;
    }
}

/*--------------------------------------------------------------*/
/* 2D Forward Dyadic Wavelet Transform at 1 scale (3-way)       */
/*--------------------------------------------------------------*/

/*
   subband    filters
    0  1       LL H0
    2  1,2     0H

    1) X low-filtering              of the signal   into the buffer
    2) Y low-filtering              of the buffer   into subband 0
    3) X high-filtering             of the signal   into subband 1
    4) Y high-filtering             of the signal   into subband 2

    All filters are of the same scale.
*/

void triway_forward(FLOAT *signal, int x_size, int y_size, int scale,
    FLOAT * * trans_data,
    FLOAT * buffer, FLOAT *l, int loff, FLOAT *h, int hoff,
    int low_high_inv_x, int low_high_inv_y)
{
    int offset;
    FLOAT * signal_off;

    FLOAT * trans_data_0_off, * trans_data_1_off, * trans_data_2_off;
    FLOAT * buffer_off;

    FLOAT * l_x, * h_x, * l_y, * h_y;
    int loff_x, hoff_x, loff_y, hoff_y;

    CWPT2_LOW_HIGH_INV

    /* analysis */

    /* step 1: X low-pass signal -> buffer */
    buffer_off = buffer;
    signal_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(buffer_off, signal_off, x_size, l_x, loff_x, scale);
        buffer_off += x_size;
        signal_off += x_size;
    }

    /* step 2: Y low-pass buffer -> subband 0 */
    trans_data_0_off = trans_data[0];
    buffer_off = buffer;
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(trans_data_0_off,   buffer_off, y_size, l_y, loff_y, scale, x_size);
        trans_data_0_off ++;
        buffer_off ++;
    }

    /* step 3: X high-pass signal -> subband 1 */
    trans_data_1_off = trans_data[1];
    signal_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(trans_data_1_off, signal_off, x_size, h_x, hoff_x, scale);
        trans_data_1_off += x_size;
        signal_off += x_size;
    }

    /* step 4: Y high-pass signal -> subband 2 */
    trans_data_2_off = trans_data[2];
    signal_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        atrous_stride(trans_data_2_off, signal_off, y_size, h_y, hoff_y, scale, x_size);
        trans_data_2_off ++;
        signal_off ++;
    }

}

/*--------------------------------------------------------------*/
/* 2D Reverse Dyadic Wavelet Transform at 1 scale (3-way)       */
/*--------------------------------------------------------------*/

/*
   subband    filters               true subband
    0  1       LL H0                   0  1
    2  1,2     0H                      2  3

    1) X analysis low-filtering              of the subband 2   into the buffer 0   (true subband 2)
    2) Y reconstruction high-filtering       of the buffer  0   into the buffer 1
    3) Y reconstruction low-filtering        of the subband 0   into the buffer 1 (+)
    4) X reconstruction low-filtering        of the buffer  1   into the signal

    5) Y analysis low-filtering              of the subband 1   into the buffer 0   (true subband 1)
    6) Y reconstruction low-filtering        of the buffer  0   into the buffer 1

    7) X analysis high-filtering / 2         of the subband 2   into the buffer 0   (1 half of the true subband 3)
    8) Y analysis high-filtering / 2         of the subband 1   into the buffer 0 (+) (2 half of the true subband 3)
    9) Y reconstruction high-filtering       of the buffer  0   into the buffer 1 (+)
   10) X reconstruction high-filtering       of the buffer  1   into the signal (+)

    All filters are of the same scale.
*/

/*
    NOTE: THE PROVIDED ANALYSIS HIGH-PASS FILTER MUST BE ALREADY DIVIDED BY 2
*/    
void triway_reverse(
    FLOAT * * trans_data,
    int x_size, int y_size, int scale,
    FLOAT * signal, FLOAT * * buffer,
    FLOAT *la, int laoff, FLOAT *ha, int haoff,
    FLOAT *lr, int lroff, FLOAT *hr, int hroff,
    int low_high_inv_x, int low_high_inv_y)
{
    int offset;

    FLOAT * source_off;
    FLOAT * source2_off;
    FLOAT * dest_off;

    FLOAT * la_x, * ha_x, * la_y, * ha_y;
    int laoff_x, haoff_x, laoff_y, haoff_y;
    FLOAT * lr_x, * hr_x, * lr_y, * hr_y;
    int lroff_x, hroff_x, lroff_y, hroff_y;

    CWPT2_LOW_HIGH_INV_A
    CWPT2_LOW_HIGH_INV_R

    /* reconstruction */

/*  1) X analysis low-filtering              of the subband 2   into the buffer 0   (true subband 2)                   */
    source_off = trans_data[2];
    dest_off = buffer[0];
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(dest_off, source_off, x_size, la_x, laoff_x, scale);
        dest_off += x_size;
        source_off += x_size;
    }

/*  2) Y reconstruction high-filtering       of the buffer  0   into the buffer 1                                      */
/*  3) Y reconstruction low-filtering        of the subband 0   into the buffer 1 (+)                                  */
    source_off = buffer[0];
    source2_off = trans_data[0];
    dest_off = buffer[1];
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(dest_off, source_off, y_size, hr_y, hroff_y, scale, x_size);
        add_atrous_stride(dest_off, source2_off, y_size, lr_y, lroff_y, scale, x_size);
        dest_off ++;
        source_off ++;
        source2_off ++;
    }

/*  4) X reconstruction low-filtering        of the buffer  1   into the signal                                        */
    source_off = buffer[1];
    dest_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(dest_off, source_off, x_size, lr_x, lroff_x, scale);
        dest_off += x_size;
        source_off += x_size;
    }

/*  5) Y analysis low-filtering              of the subband 1   into the buffer 0   (true subband 1)                   */
    source_off = trans_data[1];
    dest_off = buffer[0];
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(dest_off, source_off, y_size, la_y, laoff_y, scale, x_size);
        dest_off ++;
        source_off ++;
    }

/*  6) Y reconstruction low-filtering        of the buffer  0   into the buffer 1                                      */
    source_off = buffer[0];
    dest_off = buffer[1];
    for (offset = 0; offset < x_size; offset++)
    {
        atrous_stride(dest_off, source_off, y_size, lr_y, lroff_y, scale, x_size);
        dest_off ++;
        source_off ++;
    }

/*  7) X analysis high-filtering / 2         of the subband 2   into the buffer 0   (1 half of the true subband 3)     */
    source_off = trans_data[2];
    dest_off = buffer[0];
    for (offset = 0; offset < y_size; offset++)
    {
        atrous(dest_off, source_off, x_size, ha_x, haoff_x, scale);
        dest_off += x_size;
        source_off += x_size;
    }

/*  8) Y analysis high-filtering / 2         of the subband 1   into the buffer 0 (+) (2 half of the true subband 3)   */
    source_off = trans_data[1];
    dest_off = buffer[0];
    for (offset = 0; offset < x_size; offset++)
    {
        add_atrous_stride(dest_off, source_off, y_size, ha_y, haoff_y, scale, x_size);
        dest_off ++;
        source_off ++;
    }

/*  9) Y reconstruction high-filtering       of the buffer  0   into the buffer 1 (+)                                  */
    source_off = buffer[0];
    dest_off = buffer[1];
    for (offset = 0; offset < x_size; offset++)
    {
        add_atrous_stride(dest_off, source_off, y_size, hr_y, hroff_y, scale, x_size);
        dest_off ++;
        source_off ++;
    }

/* 10) X reconstruction high-filtering       of the buffer  1   into the signal (+)                                    */
    source_off = buffer[1];
    dest_off = signal;
    for (offset = 0; offset < y_size; offset++)
    {
        add_atrous(dest_off, source_off, x_size, hr_x, hroff_x, scale);
        dest_off += x_size;
        source_off += x_size;
    }

}
