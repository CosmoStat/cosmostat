/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#ifndef __DYADIC2_H_INCLUDED
#define __DYADIC2_H_INCLUDED

/* swaps b and c. a must be of the same type */
#define CWPT2_SWAP(a,b,c)   a=b;b=c;c=a

#define MAX_BUFFER_HEAP_SIZE            (4*CWPT2_SCALE_MAX)

/* The maximum number of necessary buffers:
   for forward transform: 0, 4, 8, 11, 14, 17, ...
   for reverse transform: 0, 4, 8, 12, 16, 20, ...
   (for nb_scales = 1, 2, 3, 4, 5, 6, ...)

   The heap implementation below performs malloc only if necessary (all allocated
   zones are being used). One should "release" a memory block if it is not used.
   Effective memory freeing is performed when destroying the heap.

   Note: relies on the value of 'sig_size' for the size of the memory zones.
*/

#define DECLARE_BUFFER_HEAP int current_buffer_heap_size;       \
    int buffer_heap_runner;                                     \
    FLOAT * * buffer_heap;                                      \
    int buffer_heap_full;                                       \
    int * buffer_heap_states;

#define MALLOC_BUFFER_HEAP current_buffer_heap_size = 0;                            \
    buffer_heap = NULL;                                                             \
    buffer_heap_states = NULL;                                                      \
    buffer_heap = (FLOAT * *)malloc(sizeof(FLOAT *) * MAX_BUFFER_HEAP_SIZE);        \
    if (!buffer_heap)                                                               \
        goto quit;                                                                  \
    buffer_heap_states = (int *)malloc(sizeof(int) * MAX_BUFFER_HEAP_SIZE);         \
    if (!buffer_heap_states)                                                        \
        goto quit;                                                                  \
    for (buffer_heap_runner=0;buffer_heap_runner<MAX_BUFFER_HEAP_SIZE;buffer_heap_runner++)   \
    {                                                                                         \
        buffer_heap[buffer_heap_runner] = NULL;                                               \
        buffer_heap_states[buffer_heap_runner] = 1;                                           \
    }

#define BUFFER_HEAP_GET(a) buffer_heap_full = 1;                                                    \
        for (buffer_heap_runner=0; buffer_heap_runner<current_buffer_heap_size; buffer_heap_runner++)    \
            if (buffer_heap_states[buffer_heap_runner] == 0)                                        \
            {                                                                                       \
                buffer_heap_states[buffer_heap_runner] = 1;                                         \
                a = buffer_heap_runner;                                                             \
                buffer_heap_full = 0;                                                               \
                break;                                                                              \
            }                                                                                       \
        if (buffer_heap_full)                                                                       \
        {                                                                                           \
            if (current_buffer_heap_size == MAX_BUFFER_HEAP_SIZE)                                   \
                goto quit;                                                                          \
            a = current_buffer_heap_size;                                                           \
            buffer_heap[current_buffer_heap_size] = (FLOAT *) malloc(sig_size * sizeof(FLOAT));     \
            if (!buffer_heap[buffer_heap_runner])                                                   \
                goto quit;                                                                          \
            current_buffer_heap_size ++;                                                            \
        }

#define BUFFER_HEAP_RELEASE(a) buffer_heap_states[(a)] = 0;

#define FREE_BUFFER_HEAP if (buffer_heap)                                                           \
    {                                                                                               \
        for (buffer_heap_runner=0;buffer_heap_runner<MAX_BUFFER_HEAP_SIZE;buffer_heap_runner++)     \
            if (buffer_heap[buffer_heap_runner])                                                    \
                free(buffer_heap[buffer_heap_runner]);                                              \
        free(buffer_heap);                                                                          \
    }                                                                                               \
    if (buffer_heap_states)                                                                         \
        free(buffer_heap_states);

#define CWPT2_LOW_HIGH_INV_PARAM(l_x,loff_x,h_x,hoff_x,l_y,loff_y,h_y,hoff_y,l,loff,h,hoff)   if (low_high_inv_x)    \
    {                                               \
        l_x    = h;                                 \
        loff_x = hoff;                              \
        h_x    = l;                                 \
        hoff_x = loff;                              \
    }                                               \
    else                                            \
    {                                               \
        l_x    = l;                                 \
        loff_x = loff;                              \
        h_x    = h;                                 \
        hoff_x = hoff;                              \
    }                                               \
    if (low_high_inv_y)                             \
    {                                               \
        l_y    = h;                                 \
        loff_y = hoff;                              \
        h_y    = l;                                 \
        hoff_y = loff;                              \
    }                                               \
    else                                            \
    {                                               \
        l_y    = l;                                 \
        loff_y = loff;                              \
        h_y    = h;                                 \
        hoff_y = hoff;                              \
    }

#define CWPT2_LOW_HIGH_INV   CWPT2_LOW_HIGH_INV_PARAM(l_x,loff_x,h_x,hoff_x,l_y,loff_y,h_y,hoff_y,l,loff,h,hoff)
#define CWPT2_LOW_HIGH_INV_A CWPT2_LOW_HIGH_INV_PARAM(la_x,laoff_x,ha_x,haoff_x,la_y,laoff_y,ha_y,haoff_y,la,laoff,ha,haoff)
#define CWPT2_LOW_HIGH_INV_R CWPT2_LOW_HIGH_INV_PARAM(lr_x,lroff_x,hr_x,hroff_x,lr_y,lroff_y,hr_y,hroff_y,lr,lroff,hr,hroff)

int dyadic2forward(FLOAT * signal, int x_size, int y_size, char * btree, int btree_size, FLOAT * trans_data, FLOAT *l, int loff, FLOAT *h, int hoff);
int dyadic2reverse(FLOAT * trans_data, int x_size, int y_size, char * btree, int btree_size, FLOAT * signal, FLOAT *l, int loff, FLOAT *h, int hoff, int first_scale_shifted);

int forward3(FLOAT *signal, int x_size, int y_size, char * btree, int btree_size,
    FLOAT * trans_data, FLOAT *l, int loff, FLOAT *h, int hoff);
int reverse3(FLOAT * trans_data, int x_size, int y_size, char * btree, int btree_size,
    FLOAT *signal, FLOAT *la, int laoff, FLOAT *ha, int haoff,
    FLOAT *lr, int lroff, FLOAT *hr, int hroff, int first_scale_shifted);

int ForwardDyadicTransform(FLOAT *signal, int sigsize,
    int scales, FLOAT * * subbands, FLOAT *l, int loff, FLOAT *h, int hoff);
int ReverseDyadicTransform(FLOAT **subbands, int sigsize, int scales,
    FLOAT *signal, FLOAT *buffer, FLOAT *l, int loff, FLOAT *h, int hoff);

#endif
