/*
 2D Continuous Wavelet Packet Transform package
 (c) 2002-2005 Let It Wave, all rights reserved
*/

#ifndef __BTREE_H_INCLUDED
#define __BTREE_H_INCLUDED

#define BTREE_NB_FIELDS         8

int btree_get_size(int nb_scales);
int btree_get_nb_scales(int btree_size);
int btree_get_next_size(int btree_size);

int btree_get_nb_leaves(char * btree, int btree_size, int children_per_parent);
void btree_fill_full(char * btree, int nb_scales);
void btree_fill_wavelet(char * btree, int nb_scales);
int btree_is_valid(char * btree, int btree_size);
int btree_get_leaf_order(char * btree, int btree_size, int * leaf_order, int nb_leaves);
int btree_fill_living_branches(char * btree, int btree_size, char value);

int * btree_create_stack(int btree_stack_size);
int btree_is_stack_empty(int * btree_stack);
int btree_push(int * btree_stack, int f1, int f2, int f3, int f4, int f5, int f6, int f7, int f8);
int btree_pop(int * btree_stack, int * f1, int * f2, int * f3, int * f4, int * f5, int * f6, int * f7, int * f8);
void btree_free_stack(int * btree_stack);

#endif
