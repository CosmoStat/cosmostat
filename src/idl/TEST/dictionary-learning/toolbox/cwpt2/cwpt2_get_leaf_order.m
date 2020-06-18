function leaf_order = cwpt2_get_leaf_order(btree, nb_children)
% function leaf_order = cwpt2_get_leaf_order(btree, nb_children)
%
% Returns the throughout order of the leaves (subbands of the transformed
% signal corresponding to the terminal nodes) defined by a given btree.
% This order corresponds to the order of the results provided by the
% cwpt2 function.
%
% btree (B x 1 or 1 x B uint8 vector)
%     The btree structure to be used by cwpt2 and cwpt2i functions.
%
% nb_children (real double scalar)
%     The number of children per parent. Must be set to 3 or 4.
%
% leaf_order (L x 1 int32 vector)
%     The leaf order. For each leaf, the leaf_order array gives its number
%     corresponding to the throughout enumeration of the full transform tree.
%     Its enumeration is defined as follows:
%
%     ===== FOR QUAD-TREES =====
%     root = 0
%     for a node N, its children (LL, HL, LH, HH) are marked
%                                (4*N+1, 4*N+2, 4*N+3, 4*N+4)
%
%     i.e.
%     scale 1 (LL, HL, LH, HH) = 1, 2, 3, 4
%     scale 2 (LL1 -> (LL,HL,LH,HH), HL1 -> (LL,..,HH), LH1 -> (..), HH1 -> (..)) =
%                      5, 6, 7, 8            9, ..,12,  13, .. 16,   17, .. 20,
%     scale 3 = 21 .. 84
%     ...
%
%     ===== FOR TRI-TREES =====
%     root = 0
%     for a node N, its children (LL, H0, 0H) are marked
%                                (4*N+1, 4*N+2, 4*N+3)
%
%     i.e.
%     scale 1 (LL, H0, 0H) = 1, 2, 3
%     scale 2 (LL1 -> (LL,H0,0H   ), H01 -> (LL,..,0H     ), 0H1 -> (..)) =
%                      5, 6, 7,              9, ..,11,       13, .. 15,
%     scale 3 = 21 .. 63
%     scale 4 = 85 ...
%     ...
%
%     Note: the first filtering direction is the first
%     memory raster direction of data. For instance, LH means for Matlab images
%     low-filtering along columns (first dimension) and high-filtering along
%     rows (second dimension).
%
%     The leaf order is computed as follows:
%     one performs the depth left-to-right search of the transform tree. Subbands are
%     added to the output N-D matrix as soon as they are computed, so leaves
%     of lower scales are always output earlier than leaves of higher scales
%     (which are children of the processed node). "Left-to-right" means
%     the LL, HL, LH, HH order (or LL, H0, 0H for 3-way transform).
%
% See also:
%     cwpt2, cwpt2_btree
%

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
