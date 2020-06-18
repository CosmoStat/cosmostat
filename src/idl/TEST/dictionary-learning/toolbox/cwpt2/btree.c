/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#include "cwpt2.h"

/* gets total number of elements for a btree structure at given maximum number of scales */
int btree_get_size(int nb_scales)
{
    int j;
    int btree_size;
    int btree_j_size;

    if (nb_scales < 0 || nb_scales > CWPT2_SCALE_MAX)
        return 0;

    btree_size = 1;
    btree_j_size = 1;

    for (j=1; j <= nb_scales; j++)
    {
        btree_j_size *= 4;
        btree_size += btree_j_size;
    }

    return btree_size;
}

int btree_get_nb_scales(int btree_size)
{
    int j;

    for (j=1; j<=CWPT2_SCALE_MAX; j++)
        if (btree_size == btree_get_size(j))
            return j;

    return 0;
}

int btree_get_next_size(int btree_size)
{
    return btree_get_size(btree_get_nb_scales(btree_size) + 1);
}

/* gets the total number of leaves (terminal nodes) = 1 + (children_per_parent - 1) * (number of divisions) */
int btree_get_nb_leaves(char * btree, int btree_size, int children_per_parent)
{
    int i;
    int nb_leaves;

    nb_leaves = 1;
    for (i=0; i<btree_size; i++)
        if (btree[i])
            nb_leaves+=(children_per_parent-1);

    return nb_leaves;
}

/* filling an (already created) btree array with values corresponding to full decomposition */
void btree_fill_full(char * btree, int nb_scales)
{
    int j;
    int j_start_index, j_curr_index;

    for (j=0, j_start_index = 0; j<nb_scales; j++)
    {
        j_curr_index = j_start_index;
        j_start_index *= 4;
        j_start_index ++;
        for (; j_curr_index < j_start_index; j_curr_index ++)
            btree[j_curr_index] = 1;
    }

    /* deepest scale : filled by zeros */
    j_curr_index = j_start_index;
    j_start_index *= 4;
    j_start_index ++;
    for (; j_curr_index < j_start_index; j_curr_index ++)
        btree[j_curr_index] = 0;
}

/* filling an (already created) btree array with values corresponding to standard wavelet decomposition */
void btree_fill_wavelet(char * btree, int nb_scales)
{
    int j;
    int j_start_index, j_curr_index;
    for (j=0, j_start_index = 0; j<nb_scales; j++)
    {
        btree[j_start_index] = 1;
        j_curr_index = j_start_index + 1;
        j_start_index *= 4;
        j_start_index ++;
        for (; j_curr_index < j_start_index; j_curr_index ++)
            btree[j_curr_index] = 0;
    }

    /* deepest scale : filled by zeros */
    j_curr_index = j_start_index;
    j_start_index *= 4;
    j_start_index ++;
    for (; j_curr_index < j_start_index; j_curr_index ++)
        btree[j_curr_index] = 0;
}

/* checks the validity of a given btree structure */
int btree_is_valid(char * btree, int btree_size)
{
    int i, j;
    char * btree_copy;

    j = btree_get_nb_scales(btree_size);
    if (j == 0)         /* invalid btree_size */
        return 0;

    if (btree[0] != 1)                              /* the root MUST be divided */
        return 0;

    for (i=0; i<btree_size; i++)
        if ((btree[i]!=0) && (btree[i]!=1))         /* bad values */
            return 0;

    /* all nodes at the last scale should be terminal */

    for (i=btree_get_size(j-1); i<btree_size; i++)
        if (btree[i] == 1)
            return 0;

    /* checking that there are no "1" in the dead branches */
    btree_copy = (char *)malloc(sizeof(char) * btree_size);
    if (!btree_copy)
        return 0;           /* malloc error */
    memcpy(btree_copy, btree, btree_size);
    if (btree_fill_living_branches(btree_copy, btree_size, 2))
    {
        free(btree_copy);
        return 0;           /* internal error */
    }

    for (i=0; i<btree_size; i++)
        if (btree_copy[i]==1)         /* 1 in a dead branch */
            return 0;

    free(btree_copy);

    /* at least one node at the scale before the last should be non-terminal */

    for (i=btree_get_size(j-2); i<btree_get_size(j-1); i++)
        if (btree[i] == 1)
            return 1;

    return 0;
}

int btree_get_leaf_order(char * btree, int btree_size, int * leaf_order, int nb_leaves)
{
    int * btree_stack;
    int result;
    int notused;
    int child_no, node_to_process, next_final;

    result = 1;         /* default: malloc failure */

    btree_stack = btree_create_stack(btree_size);
    if (!btree_stack)
        goto quit;

    result = 2;         /* default : buffer failure */

    next_final = 0;
    if (btree_push(btree_stack, 0, -1, -1, -1, -1, -1, -1, -1))                  /* to process: the tree root */
        goto quit;

    while(!btree_is_stack_empty(btree_stack))
    {
        if (btree_pop(btree_stack, &node_to_process, &notused, &notused, &notused, &notused, &notused, &notused, &notused))
            goto quit;

        /* Terminal nodes will be numbered from left to right */
        for (child_no = 1; child_no <= 4; child_no ++)
            if ( btree [ node_to_process * 4 + child_no ] == 0 )       /* the child node is terminal */
            {
                if (next_final == nb_leaves)
                    goto quit;                      /* should NEVER happen */
                leaf_order[next_final++] = node_to_process * 4 + child_no;
            }

        /* Non-terminal nodes will be put into stack from right to left so that the left ones are processed first */
        for (child_no = 4; child_no >= 1; child_no --)
            if ( btree [ node_to_process * 4 + child_no ] != 0 )       /* the child node is non-terminal */
                if (btree_push(btree_stack, node_to_process * 4 + child_no, -1, -1, -1, -1, -1, -1, -1 ))
                    goto quit;

    }

    if (next_final != nb_leaves)        /* should NEVER happen either: expecting more leaves... */
        goto quit;

    result = 0;         /* succeeded */

quit:
    btree_free_stack(btree_stack);

    return result;
}

int btree_fill_living_branches(char * btree, int btree_size, char value)
{
    int * btree_stack;
    int result;
    int notused;
    int child_no, node_to_process;

    result = 1;         /* default: malloc failure */

    btree_stack = btree_create_stack(btree_size);
    if (!btree_stack)
        goto quit;

    result = 2;         /* default : buffer failure */

    if (btree_push(btree_stack, 0, -1, -1, -1, -1, -1, -1, -1))                  /* to process: the tree root */
        goto quit;

    while(!btree_is_stack_empty(btree_stack))
    {
        if (btree_pop(btree_stack, &node_to_process, &notused, &notused, &notused, &notused, &notused, &notused, &notused))
            goto quit;

        btree[node_to_process] = value;

        /* Non-terminal nodes will be put into stack from right to left so that the left ones are processed first */
        for (child_no = 4; child_no >= 1; child_no --)
            if ( btree [ node_to_process * 4 + child_no ] != 0 )       /* the child node is non-terminal */
                if (btree_push(btree_stack, node_to_process * 4 + child_no, -1, -1, -1, -1, -1, -1, -1 ))
                    goto quit;

    }

    result = 0;         /* succeeded */

quit:
    btree_free_stack(btree_stack);

    return result;
}

/* BTREE STACK ********************************************************/

int * btree_create_stack(int btree_stack_size)
{
    int * btree_stack;
    /* 5 = number of fields */
    btree_stack = (int *) malloc (sizeof(int) * (btree_stack_size * BTREE_NB_FIELDS + 3));
    if (btree_stack)
    {
        btree_stack[0] = btree_stack_size;
        btree_stack[1] = 0;                     /* current number of records */
        btree_stack[2] = BTREE_NB_FIELDS;
    }

    return btree_stack;
}

int btree_is_stack_empty(int * btree_stack)
{
    if (btree_stack[1] == 0)
        return 1;
    return 0;
}

int btree_push(int * btree_stack, int f1, int f2, int f3, int f4, int f5, int f6, int f7, int f8)
{
    int * current_record;
    if (btree_stack[1] == btree_stack[0])
        return 1;                   /* error : stack is full */

    btree_stack[1] ++;
    current_record = btree_stack + (BTREE_NB_FIELDS * btree_stack[1] + 3);

    current_record[0] = f1;
    current_record[1] = f2;
    current_record[2] = f3;
    current_record[3] = f4;
    current_record[4] = f5;
    current_record[5] = f6;
    current_record[6] = f7;
    current_record[7] = f8;

    return 0;
}

int btree_pop(int * btree_stack, int * f1, int * f2, int * f3, int * f4, int * f5, int * f6, int * f7, int * f8)
{
    int * current_record;
    if (btree_stack[1] == 0)
        return 1;                   /* error : stack is empty */

    current_record = btree_stack + (BTREE_NB_FIELDS * btree_stack[1] + 3);

    * f1 = current_record[0];
    * f2 = current_record[1];
    * f3 = current_record[2];
    * f4 = current_record[3];
    * f5 = current_record[4];
    * f6 = current_record[5];
    * f7 = current_record[6];
    * f8 = current_record[7];

    btree_stack[1] --;

    return 0;
}

void btree_free_stack(int * btree_stack)
{
    if (btree_stack)
        free(btree_stack);
}
