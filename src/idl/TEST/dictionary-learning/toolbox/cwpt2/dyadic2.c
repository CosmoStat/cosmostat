/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#include "cwpt2.h"

FILE * debug_file;
int debug_mode = 0;

/*--------------------------------------------------------------*/
/* Renormalized filters for dyadic filtering                    */
/*--------------------------------------------------------------*/

/* Analysis low-pass filter                                     */

static FLOAT bi9[]={
 0.03782845550699547033,
-0.02384946501937996663,
-0.1106244044184234859,
 0.3774028556126538536,
 0.8526986790094033264,
 0.3774028556126538536,
-0.1106244044184233194,
-0.02384946501937999092,
 0.03782845550699546339
};

/* Analysis high-pass filter                                    */

static FLOAT bi7h[]={
 0.06453888262893847649,
-0.04068941760955852721,
-0.4180922732222120963,
 0.7884856164056645023,
-0.4180922732222124294,
-0.04068941760955835374,
 0.06453888262893842098,
};

/* Reconstruction low-pass filter                               */

static FLOAT bi7[]={
-0.03226944131446923825,
-0.02034470880477926361,
 0.2090461366111060482,
 0.3942428082028322511,
 0.2090461366111062147,
-0.02034470880477917687,
-0.03226944131446921049,
};

/* Reconstruction high-pass filter                              */

static FLOAT bi9h[]={
 0.01891422775349773516,
 0.01192473250968998331,
-0.05531220220921174296,
-0.1887014278063269268,
 0.4263493395047016632,
-0.1887014278063269268,
-0.0553122022092116597,
 0.01192473250968999546,
 0.01891422775349773169,
};

/* Offsets (index of the centre coefficient)                    */

static int off9=4, off9h=4, off7=3, off7h=3;


/*
    Forward Continuous Wavelet Packet Transform routine.
    Supports arbitrary btree structures.
*/

int dyadic2forward(FLOAT *signal, int x_size, int y_size, char * btree, int btree_size,
    FLOAT * trans_data, FLOAT *l, int loff, FLOAT *h, int hoff)
{
    int result;
    int sig_size;
    int notused;

    DECLARE_BUFFER_HEAP

    FLOAT * buffer = NULL;

    int scale;              /* 1, 2, 4, 8, ... */
    FLOAT * dst[4];
    FLOAT * src;
    int low_high_inv_x, low_high_inv_y;

    int * btree_stack = NULL;
    int node_to_process;
    int src_buf_no, dst_buf_no;
    int child_no;
    int next_final;

if (debug_mode) debug_file = fopen("debug.txt","at");
if (debug_mode) fprintf(debug_file,"FORWARD===\n");

    if ((l==NULL) || (h==NULL))
    {
        l = bi9;
        loff = off9;
        h = bi7h;
        hoff = off7h;
    }

    result = 1;         /* default: malloc failure */
    MALLOC_BUFFER_HEAP

    sig_size = x_size * y_size;
    buffer = (FLOAT *) malloc(sig_size * sizeof(FLOAT));
    if (buffer == NULL)
        goto quit;

    btree_stack = btree_create_stack(btree_size);
    if (!btree_stack)
        goto quit;

    result = 2;         /* default : buffer failure */

    next_final = 0;
    if (btree_push(btree_stack, 0, 1, -1, 0, 0, -1, -1, -1))                  /* to process: the tree root */
        goto quit;

    while(!btree_is_stack_empty(btree_stack))
    {
        if (btree_pop(btree_stack, &node_to_process, &scale, &src_buf_no, &low_high_inv_x, &low_high_inv_y, &notused, &notused, &notused))
            goto quit;

        if (debug_mode) fprintf(debug_file,"popping node %d, scale %d, inv %d %d\n",node_to_process, scale, low_high_inv_x,low_high_inv_y);

        /* choosing source ... */
        if (src_buf_no == -1)
            src = signal;
        else
            src = buffer_heap[src_buf_no];

        /* choosing destinations ... */
        /* children: 1 2
                     3 4 where top left corner corresponds to low frequencies */

        /* Terminal nodes will be numbered from left to right */
        for (child_no = 1; child_no <= 4; child_no ++)
            if ( btree [ node_to_process * 4 + child_no ] == 0 )       /* the child node is terminal */
            {
                if (debug_mode) fprintf(debug_file,"found terminal node: %d\n",node_to_process * 4 + child_no);
                dst[child_no - 1] = trans_data + next_final * sig_size;
                next_final ++;
            }

        /* Non-terminal nodes will be put into stack from right to left so that the left ones are processed first */
        for (child_no = 4; child_no >= 1; child_no --)
            if ( btree [ node_to_process * 4 + child_no ] != 0 )       /* the child node is non-terminal */
            {
                BUFFER_HEAP_GET(dst_buf_no)
                dst[child_no - 1] = buffer_heap[dst_buf_no];
                /* for children 1,3 (low x-filtering) there is no x filter inversion
                   for children 2,4 (high x-filtering) there is an x filter inversion
                   for children 1,2 (low y-filtering) there is no y filter inversion
                   for children 3,4 (high y_filtering) there is an y filter inversion */
                if (debug_mode) fprintf(debug_file,"pushing node %d, scale %d, inv %d %d\n",node_to_process * 4 + child_no, scale *2, (child_no + 1) % 2,(child_no - 1) / 2);
                if (btree_push(btree_stack, node_to_process * 4 + child_no, scale * 2, dst_buf_no,
                    (child_no + 1) % 2,
                    (child_no - 1) / 2,
                    -1, -1, -1
                    ))
                    goto quit;
            }

        /* processing ... */

        if (debug_mode) fprintf(debug_file,"processing at scale %d, inv %d %d\n", scale, low_high_inv_x, low_high_inv_y);

        quad_forward(src, x_size, y_size, scale, dst, buffer, l, loff, h, hoff,
            low_high_inv_x, low_high_inv_y);

        /* releasing source ... */
        if (src_buf_no != -1)
        {
            BUFFER_HEAP_RELEASE(src_buf_no)
        }
    }

    result = 0;         /* succeeded */

quit:
    btree_free_stack(btree_stack);

    if (buffer)
        free(buffer);

    FREE_BUFFER_HEAP

if (debug_mode) fclose(debug_file);

    return result;
}

/*
    Reverse Continuous Wavelet Packet Transform routine.
    Supports arbitrary btree structures.
*/

int dyadic2reverse(FLOAT * trans_data, int x_size, int y_size, char * btree, int btree_size,
    FLOAT *signal, FLOAT *l, int loff, FLOAT *h, int hoff, int first_scale_shifted)
{
    int result;
    int sig_size;

    DECLARE_BUFFER_HEAP

    FLOAT * buffer = NULL;

    int scale;              /* 1, 2, 4, 8, ... */
    FLOAT * src[4];
    FLOAT * dst;
    int low_high_inv_x, low_high_inv_y;

    int * btree_stack = NULL;
    int parent, child1, child2, child3, child4;
    int first_undef_child;
    int dst_buf_no;
    int leaf_index;

if (debug_mode) debug_file = fopen("debug.txt","at");
if (debug_mode) fprintf(debug_file,"REVERSE==========\n");

    if ((l==NULL) || (h==NULL))
    {
        l = bi7;
        loff = off7;
        h = bi9h;
        hoff = off9h;
    }

    result = 1;         /* default: malloc failure */
    MALLOC_BUFFER_HEAP

    sig_size = x_size * y_size;
    buffer = (FLOAT *) malloc(sig_size * sizeof(FLOAT));
    if (buffer == NULL)
        goto quit;

    btree_stack = btree_create_stack(btree_size);
    if (!btree_stack)
        goto quit;

    result = 2;         /* default : buffer failure */

    leaf_index = 0;

    /* analyzing children of the root :
    if further division is needed, mark them as undefined,
    otherwise (terminal node), assign the next leaf index (+1 because 0 is reserved for 'undefined') */
    if (first_scale_shifted)
    {
        /* child1: always transformed (it was root at simple resolution) */
        /* other children : no data; we should intercept the last reconstruction phase */
        child1 = child2 = child3 = child4 = 0;
    }
    else
    {
        child1 = (btree[1]) ? 0 : (++leaf_index);
        child2 = (btree[2]) ? 0 : (++leaf_index);
        child3 = (btree[3]) ? 0 : (++leaf_index);
        child4 = (btree[4]) ? 0 : (++leaf_index);
    }

    if (btree_push(btree_stack, 0, child1, child2, child3, child4, 1, 0, 0))
        goto quit;


    while(!btree_is_stack_empty(btree_stack))
    {
        if (btree_pop(btree_stack, &parent, &child1, &child2, &child3, &child4, &scale, &low_high_inv_x, &low_high_inv_y))
            goto quit;

        if (debug_mode) fprintf(debug_file,"Popped parent: %d, children : %d %d %d %d ; scale %d ; inv %d %d\n",
            parent, child1, child2, child3, child4, scale, low_high_inv_x, low_high_inv_y);

        /* children: 1 2
                     3 4 where top left corner corresponds to low frequencies */

        /* are all children of this node already defined (final leaves or already computed buffers) ? */
        if (child1 && child2 && child3 && child4)
        {
            /* in this case, we can reconstruct this quad-tree */
            /* 1) choosing sources : negative values = (buffers + 1), positive values = (final leaves + 1)*/
            src[0] = (child1 > 0) ? (trans_data + (child1 - 1) * sig_size) : (buffer_heap[ -child1 - 1]);
            src[1] = (child2 > 0) ? (trans_data + (child2 - 1) * sig_size) : (buffer_heap[ -child2 - 1]);
            src[2] = (child3 > 0) ? (trans_data + (child3 - 1) * sig_size) : (buffer_heap[ -child3 - 1]);
            src[3] = (child4 > 0) ? (trans_data + (child4 - 1) * sig_size) : (buffer_heap[ -child4 - 1]);

            /* 2) choosing destination */
            if (parent == 0 || ((parent == 1) && first_scale_shifted))      /* 1 = root if scales are shifted */
                dst = signal;
            else
            {
                BUFFER_HEAP_GET(dst_buf_no)
                if (debug_mode) fprintf(debug_file,"Got buffer %d\n", dst_buf_no);
                dst = buffer_heap[dst_buf_no];
            }

            /* 3) processing ... */

            if (debug_mode) fprintf(debug_file,"Processing at scale %d inv %d %d\n", scale, low_high_inv_x, low_high_inv_y);

            quad_reverse(src, x_size, y_size, scale, dst, buffer, l, loff, h, hoff,
                low_high_inv_x, low_high_inv_y);

            /* 4) releasing "buffer" sources if there were any ... */
            if (child1 < 0) { BUFFER_HEAP_RELEASE(-child1 - 1) if (debug_mode) fprintf(debug_file,"Released buffer %d\n", -child1 - 1); }
            if (child2 < 0) { BUFFER_HEAP_RELEASE(-child2 - 1) if (debug_mode) fprintf(debug_file,"Released buffer %d\n", -child2 - 1); }
            if (child3 < 0) { BUFFER_HEAP_RELEASE(-child3 - 1) if (debug_mode) fprintf(debug_file,"Released buffer %d\n", -child3 - 1); }
            if (child4 < 0) { BUFFER_HEAP_RELEASE(-child4 - 1) if (debug_mode) fprintf(debug_file,"Released buffer %d\n", -child4 - 1); }

            if (debug_mode) fprintf(debug_file, "BUFSTATE: %d %d \n",buffer_heap_states[0], buffer_heap_states[1]);

            if ((parent == 1) && first_scale_shifted)
                goto ok_quit;

            /* 5) modifying the info in its parent node */
            if (parent == 0)
            {
                if (!btree_is_stack_empty(btree_stack))     /* should NEVER happen : we have finished ! */
                    goto quit;
            }
            else
            {
                if (btree_pop(btree_stack, &parent, &child1, &child2, &child3, &child4, &scale, &low_high_inv_x, &low_high_inv_y))
                    goto quit;

                /* by the definition of the processing order,
                   the first undefined child of its parent (normally, there should be) is our newly created buffer */
                if (child1 == 0) child1 = -dst_buf_no - 1;
                else if (child2 == 0) child2 = -dst_buf_no - 1;
                else if (child3 == 0) child3 = -dst_buf_no - 1;
                else if (child4 == 0) child4 = -dst_buf_no - 1;
                else goto quit;         /* should NEVER happen */

                if (btree_push(btree_stack, parent, child1, child2, child3, child4, scale, low_high_inv_x, low_high_inv_y))
                    goto quit;
            }
        }
        else        /* there are some undefined children */
        {
            /* pushing back this node, we are not going to process it for the moment ...
               right now, we will push into the stack data about its first undefined child ONLY,
               other children are processed later (tree processing: depth search) */
            if (btree_push(btree_stack, parent, child1, child2, child3, child4, scale, low_high_inv_x, low_high_inv_y))
                goto quit;

            if (child1 == 0) first_undef_child = 1;
            else if (child2 == 0) first_undef_child = 2;
            else if (child3 == 0) first_undef_child = 3;
            else if (child4 == 0) first_undef_child = 4;
            else goto quit;         /* should NEVER happen */

            parent = parent * 4 + first_undef_child;

            /* analyzing children of the first undefined child of the previously popped node :
               if further division is needed, mark them as undefined,
               otherwise (terminal node), assign the next leaf index (+1 because 0 is reserved for 'undefined') */
            child1 = (btree[parent * 4 + 1]) ? 0 : (++leaf_index);
            child2 = (btree[parent * 4 + 2]) ? 0 : (++leaf_index);
            child3 = (btree[parent * 4 + 3]) ? 0 : (++leaf_index);
            child4 = (btree[parent * 4 + 4]) ? 0 : (++leaf_index);

            /* for children 1,3 (low x-filtering) there is no x filter inversion
            for children 2,4 (high x-filtering) there is an x filter inversion
            for children 1,2 (low y-filtering) there is no y filter inversion
            for children 3,4 (high y_filtering) there is an y filter inversion */

            if (btree_push(btree_stack, parent, child1, child2, child3, child4, scale * 2,
                (first_undef_child + 1) % 2,
                (first_undef_child - 1) / 2
                ))
                goto quit;
        }
    }

ok_quit:
    result = 0;         /* succeeded */

quit:
    btree_free_stack(btree_stack);

    if (buffer)
        free(buffer);

    FREE_BUFFER_HEAP

if (debug_mode) fclose(debug_file);

    return result;
}

/*
   The following functions provide the 1D Continuous WAVELET Transform support
   at a given scale. They are not used any more, but they are maintained
   for compatibility reasons with the original file: dyadic_transform.c

   Example: scales = 3:
        subbands[0] = signal * high(1)
        subbands[1] = signal * high(2)
        subbands[2] = signal * high(4)
        subbands[3] = signal * low(4)

*/

/*--------------------------------------------------------------*/
/* Forward Dyadic Wavelet Transform                             */
/* Returns the array of arrays of coefficients                  */
/*--------------------------------------------------------------*/

int ForwardDyadicTransform(FLOAT *signal, int sigsize,
    int scales, FLOAT * * subbands, FLOAT *l, int loff, FLOAT *h, int hoff)
{
    int j,scale;
    FLOAT *in=signal;

    if (l==NULL)
    {
        l = bi9;
        loff = off9;
    }
    if (h==NULL)
    {
        h = bi7h;
        hoff = off7h;
    }

    /* low-pass filterings */

    for (j=0, scale=1; j<scales-1; j++, scale*=2)
    {
        atrous(subbands[j], in, sigsize, l, loff, scale);
        in=subbands[j];
    }

    atrous(subbands[scales], in, sigsize, l, loff, scale);

    /* high-pass filterings */

    for (; j>0; j--, scale/=2)
        atrous(subbands[j], subbands[j-1], sigsize, h, hoff, scale);

    atrous(subbands[0], signal, sigsize, h, hoff, scale);

    return 0;
}

/*--------------------------------------------------------------*/
/* Reverse Dyadic Wavelet Transform                             */
/* Returns the array of samples                                 */
/*--------------------------------------------------------------*/

int ReverseDyadicTransform(FLOAT **subbands, int sigsize, int scales,
    FLOAT *signal, FLOAT *buffer, FLOAT *l, int loff, FLOAT *h, int hoff)
{
    int j,scale;
    FLOAT *current, *next, *sw;
    FLOAT *ownbuffer;

    if (l==NULL)
    {
        l = bi7;
        loff = off7;
    }

    if (h==NULL)
    {
        h = bi9h;
        hoff = off9h;
    }

    /* allocate missing array. abort and return 1 on failure */

    if (buffer)
        ownbuffer = buffer;
    else
        ownbuffer = calloc(sigsize, sizeof(FLOAT));

    if(!ownbuffer)
        return 1;           /* malloc error */

    /* buffer and index preparation */

    current = signal;
    next = ownbuffer;

    for (j=0, scale=1; j<scales; j++, scale*=2)
    {
        CWPT2_SWAP(sw, current, next);
    }

    /* reconstruction */

    j--;
    scale /= 2;

    atrous(next, subbands[scales], sigsize, l, loff, scale);
    add_atrous(next, subbands[scales-1], sigsize, h, hoff, scale);

    for (j--, scale/=2; j>=0; j--, scale/=2)
    {
        CWPT2_SWAP(sw, current, next);
        atrous(next, current, sigsize, l, loff, scale);
        add_atrous(next, subbands[j], sigsize, h, hoff, scale);
    }

quit:
    if (!buffer)
        free(ownbuffer);
    return 0;
}


/*
    The following functions implement the
    Continuous 3-way Wavelet Packet Transforms. They are based on triway_forward
    and triway_reverse functions supporting multi-scale and packets (low<->high filter inversion for x and/or y).

    However, they HAVE NOT BEEN TESTED in case of packets.
*/

int forward3(FLOAT *signal, int x_size, int y_size, char * btree, int btree_size,
    FLOAT * trans_data, FLOAT *l, int loff, FLOAT *h, int hoff)
{
    int result;
    int sig_size;
    int notused;

    DECLARE_BUFFER_HEAP

    FLOAT * buffer = NULL;

    int scale;              /* 1, 2, 4, 8, ... */
    FLOAT * dst[3];
    FLOAT * src;
    int low_high_inv_x, low_high_inv_y;

    int * btree_stack = NULL;
    int node_to_process;
    int src_buf_no, dst_buf_no;
    int child_no;
    int next_final;

if (debug_mode) debug_file = fopen("debug.txt","at");
if (debug_mode) fprintf(debug_file,"FORWARD===\n");

    if ((l==NULL) || (h==NULL))
    {
        l = bi9;
        loff = off9;
        h = bi7h;
        hoff = off7h;
    }

    result = 1;         /* default: malloc failure */
    MALLOC_BUFFER_HEAP

    sig_size = x_size * y_size;
    buffer = (FLOAT *) malloc(sig_size * sizeof(FLOAT));
    if (buffer == NULL)
        goto quit;

    btree_stack = btree_create_stack(btree_size);
    if (!btree_stack)
        goto quit;

    result = 2;         /* default : buffer failure */

    next_final = 0;
    if (btree_push(btree_stack, 0, 1, -1, 0, 0, -1, -1, -1))                  /* to process: the tree root */
        goto quit;

    while(!btree_is_stack_empty(btree_stack))
    {
        if (btree_pop(btree_stack, &node_to_process, &scale, &src_buf_no, &low_high_inv_x, &low_high_inv_y, &notused, &notused, &notused))
            goto quit;

        if (debug_mode) fprintf(debug_file,"popping node %d, scale %d, inv %d %d\n",node_to_process, scale, low_high_inv_x,low_high_inv_y);

        /* choosing source ... */
        if (src_buf_no == -1)
            src = signal;
        else
            src = buffer_heap[src_buf_no];

        /* choosing destinations ... */
        /* children: 1 2
                     3 X where top left corner corresponds to low frequencies */
        /* Still, the btree structure is the same as before (aligned at 4) */

        if ( btree [ node_to_process * 4 + 4] == 1 )       /* the fourth child is non-terminal : makes NO sense in 3-way transform */
            goto quit;

        /* Terminal nodes will be numbered from left to right */
        for (child_no = 1; child_no <= 3; child_no ++)
            if ( btree [ node_to_process * 4 + child_no ] == 0 )       /* the child node is terminal */
            {
                if (debug_mode) fprintf(debug_file,"found terminal node: %d\n",node_to_process * 4 + child_no);
                dst[child_no - 1] = trans_data + next_final * sig_size;
                next_final ++;
            }

        /* Non-terminal nodes will be put into stack from right to left so that the left ones are processed first */
        for (child_no = 3; child_no >= 1; child_no --)
            if ( btree [ node_to_process * 4 + child_no ] != 0 )       /* the child node is non-terminal */
            {
                BUFFER_HEAP_GET(dst_buf_no)
                dst[child_no - 1] = buffer_heap[dst_buf_no];
                /* for child 1 (low x- and y-filtering) there is no x nor y filter inversion
                   for child 2 (high x-filtering) there is an x filter inversion
                   for child 3 (high y_filtering) there is an y filter inversion */
                if (debug_mode) fprintf(debug_file,"pushing node %d, scale %d, inv %d %d\n",node_to_process * 4 + child_no, scale *2, (child_no + 1) % 2,(child_no - 1) / 2);
                if (btree_push(btree_stack, node_to_process * 4 + child_no, scale * 2, dst_buf_no,
                    (child_no + 1) % 2,
                    (child_no - 1) / 2,
                    -1, -1, -1
                    ))
                    goto quit;
            }

        /* processing ... */

        if (debug_mode) fprintf(debug_file,"processing at scale %d, inv %d %d\n", scale, low_high_inv_x, low_high_inv_y);

        triway_forward(src, x_size, y_size, scale, dst, buffer, l, loff, h, hoff,
            low_high_inv_x, low_high_inv_y);

        /* releasing source ... */
        if (src_buf_no != -1)
        {
            BUFFER_HEAP_RELEASE(src_buf_no)
        }
    }

    result = 0;         /* succeeded */

quit:
    btree_free_stack(btree_stack);

    if (buffer)
        free(buffer);

    FREE_BUFFER_HEAP

if (debug_mode) fclose(debug_file);

    return result;
}


int reverse3(FLOAT * trans_data, int x_size, int y_size, char * btree, int btree_size,
    FLOAT *signal,
    FLOAT *la, int laoff, FLOAT *ha_orig, int haoff,
    FLOAT *lr, int lroff, FLOAT *hr, int hroff,
    int first_scale_shifted)
{
    int result;
    int sig_size;
    int notused;

    DECLARE_BUFFER_HEAP

    FLOAT * buffer[2];
    FLOAT * ha = NULL;

    int scale;              /* 1, 2, 4, 8, ... */
    FLOAT * src[3];
    FLOAT * dst;
    int low_high_inv_x, low_high_inv_y;

    int * btree_stack = NULL;
    int parent, child1, child2, child3;
    int first_undef_child;
    int dst_buf_no;
    int leaf_index;
    int runner;

if (debug_mode) debug_file = fopen("debug.txt","at");
if (debug_mode) fprintf(debug_file,"REVERSE==========\n");

    if ((la==NULL) || (ha_orig==NULL) || (lr==NULL) || (hr==NULL))
    {
        la = lr = bi7;
        laoff = lroff = off7;
        ha_orig  = hr = bi9h;
        haoff = hroff = off9h;
    }

    result = 1;         /* default: malloc failure */
    MALLOC_BUFFER_HEAP

    /* we should pre-divide the high-pass analysis filter by 2 */
    ha = (FLOAT *) malloc ((2 * haoff + 1) * sizeof(FLOAT));
    if (ha == NULL)
        goto quit;
    for (runner = 0; runner < (2 * haoff + 1); runner ++)
        ha[runner] = ha_orig[runner] / 2;

    sig_size = x_size * y_size;
    buffer[0] = (FLOAT *) malloc(sig_size * sizeof(FLOAT));
    buffer[1] = (FLOAT *) malloc(sig_size * sizeof(FLOAT));
    if ((buffer[0] == NULL) || (buffer[1] == NULL))
        goto quit;

    btree_stack = btree_create_stack(btree_size);
    if (!btree_stack)
        goto quit;

    result = 2;         /* default : buffer failure */

    leaf_index = 0;

    /* analyzing children of the root :
    if further division is needed, mark them as undefined,
    otherwise (terminal node), assign the next leaf index (+1 because 0 is reserved for 'undefined') */
    if (first_scale_shifted)
    {
        /* child1: always transformed (it was root at simple resolution) */
        /* other children : no data; we should intercept the last reconstruction phase */
        child1 = child2 = child3 = 0;
    }
    else
    {
        child1 = (btree[1]) ? 0 : (++leaf_index);
        child2 = (btree[2]) ? 0 : (++leaf_index);
        child3 = (btree[3]) ? 0 : (++leaf_index);
        if (btree[4])
            goto quit;
    }

    if (btree_push(btree_stack, 0, child1, child2, child3, -1, 1, 0, 0))
        goto quit;


    while(!btree_is_stack_empty(btree_stack))
    {
        if (btree_pop(btree_stack, &parent, &child1, &child2, &child3, &notused, &scale, &low_high_inv_x, &low_high_inv_y))
            goto quit;

        if (debug_mode) fprintf(debug_file,"Popped parent: %d, children : %d %d %d ; scale %d ; inv %d %d\n",
            parent, child1, child2, child3, scale, low_high_inv_x, low_high_inv_y);

        /* children: 1 2
                     3 4 where top left corner corresponds to low frequencies */

        /* are all children of this node already defined (final leaves or already computed buffers) ? */
        if (child1 && child2 && child3)
        {
            /* in this case, we can reconstruct this quad-tree */
            /* 1) choosing sources : negative values = (buffers + 1), positive values = (final leaves + 1)*/
            src[0] = (child1 > 0) ? (trans_data + (child1 - 1) * sig_size) : (buffer_heap[ -child1 - 1]);
            src[1] = (child2 > 0) ? (trans_data + (child2 - 1) * sig_size) : (buffer_heap[ -child2 - 1]);
            src[2] = (child3 > 0) ? (trans_data + (child3 - 1) * sig_size) : (buffer_heap[ -child3 - 1]);

            /* 2) choosing destination */
            if (parent == 0 || ((parent == 1) && first_scale_shifted))      /* 1 = root if scales are shifted */
                dst = signal;
            else
            {
                BUFFER_HEAP_GET(dst_buf_no)
                if (debug_mode) fprintf(debug_file,"Got buffer %d\n", dst_buf_no);
                dst = buffer_heap[dst_buf_no];
            }

            /* 3) processing ... */

            if (debug_mode) fprintf(debug_file,"Processing at scale %d inv %d %d\n", scale, low_high_inv_x, low_high_inv_y);

            triway_reverse(src, x_size, y_size, scale, dst, buffer,
                la, laoff, ha, haoff, lr, lroff, hr, hroff,
                low_high_inv_x, low_high_inv_y);

            /* 4) releasing "buffer" sources if there were any ... */
            if (child1 < 0) { BUFFER_HEAP_RELEASE(-child1 - 1) if (debug_mode) fprintf(debug_file,"Released buffer %d\n", -child1 - 1); }
            if (child2 < 0) { BUFFER_HEAP_RELEASE(-child2 - 1) if (debug_mode) fprintf(debug_file,"Released buffer %d\n", -child2 - 1); }
            if (child3 < 0) { BUFFER_HEAP_RELEASE(-child3 - 1) if (debug_mode) fprintf(debug_file,"Released buffer %d\n", -child3 - 1); }

            if (debug_mode) fprintf(debug_file, "BUFSTATE: %d %d \n",buffer_heap_states[0], buffer_heap_states[1]);

            if ((parent == 1) && first_scale_shifted)
                goto ok_quit;

            /* 5) modifying the info in its parent node */
            if (parent == 0)
            {
                if (!btree_is_stack_empty(btree_stack))     /* should NEVER happen : we have finished ! */
                    goto quit;
            }
            else
            {
                if (btree_pop(btree_stack, &parent, &child1, &child2, &child3, &notused, &scale, &low_high_inv_x, &low_high_inv_y))
                    goto quit;

                /* by the definition of the processing order,
                   the first undefined child of its parent (normally, there should be) is our newly created buffer */
                if (child1 == 0) child1 = -dst_buf_no - 1;
                else if (child2 == 0) child2 = -dst_buf_no - 1;
                else if (child3 == 0) child3 = -dst_buf_no - 1;
                else goto quit;         /* should NEVER happen */

                if (btree_push(btree_stack, parent, child1, child2, child3, -1, scale, low_high_inv_x, low_high_inv_y))
                    goto quit;
            }
        }
        else        /* there are some undefined children */
        {
            /* pushing back this node, we are not going to process it for the moment ...
               right now, we will push into the stack data about its first undefined child ONLY,
               other children are processed later (tree processing: depth search) */
            if (btree_push(btree_stack, parent, child1, child2, child3, -1, scale, low_high_inv_x, low_high_inv_y))
                goto quit;

            if (child1 == 0) first_undef_child = 1;
            else if (child2 == 0) first_undef_child = 2;
            else if (child3 == 0) first_undef_child = 3;
            else goto quit;         /* should NEVER happen */

            parent = parent * 4 + first_undef_child;

            /* analyzing children of the first undefined child of the previously popped node :
               if further division is needed, mark them as undefined,
               otherwise (terminal node), assign the next leaf index (+1 because 0 is reserved for 'undefined') */
            child1 = (btree[parent * 4 + 1]) ? 0 : (++leaf_index);
            child2 = (btree[parent * 4 + 2]) ? 0 : (++leaf_index);
            child3 = (btree[parent * 4 + 3]) ? 0 : (++leaf_index);
            if (btree[parent * 4 + 4])
                goto quit;

            /* for children 1,3 (low x-filtering) there is no x filter inversion
            for child 2 (high x-filtering) there is an x filter inversion
            for children 1,2 (low y-filtering) there is no y filter inversion
            for child 3 (high y_filtering) there is an y filter inversion */

            if (btree_push(btree_stack, parent, child1, child2, child3, -1, scale * 2,
                (first_undef_child + 1) % 2,
                (first_undef_child - 1) / 2
                ))
                goto quit;
        }
    }

ok_quit:
    result = 0;         /* succeeded */

quit:
    btree_free_stack(btree_stack);

    if (ha)
        free(ha);

    if (buffer[0])
        free(buffer[0]);
    if (buffer[1])
        free(buffer[1]);

    FREE_BUFFER_HEAP

if (debug_mode) fclose(debug_file);

    return result;
}
