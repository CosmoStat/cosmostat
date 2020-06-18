function transformed_image = c3wpt(original_image, btree, alf, alf_off, ahf, ahf_off)
% function transformed_image = c3wpt(original_image, btree, alf, alf_off, ahf, ahf_off)
%
% Forward 2D 3-way Continuous Wavelet Packet Transform.
%
% original_image (M x N real double matrix)
%     The image to be transformed. 2 <= M, N <= 16384
%
% btree (B x 1 or 1 x B uint8 vector)
%     The btree structure defining the wavelet packet transform tree.
%     It should be defined as follows:
%
%     (attention, C-like array element enumeration)
%
%     btree[n] determines if the node  n  should be transformed further.
%              0 = terminal node
%              1 = non-terminal node, should be transformed.
%
%     btree[0] = 1 (corresponds to the root of the transform tree.
%                   it is always transformed.)
%
%     The children of a node  n  are nodes 4*n+1, 4*n+2, 4*n+3
%     (their order: LL, H0, 0H; see below). So, scale 1 nodes
%     are 1..3; scale 2: 5..7,9..11,13..15; scale 3: 21..63; scale 4: 85...
%
%     Attention: if the maximum scale of the transform tree is K, one
%     should provide a K-scale btree where the last scale can not contain
%     non-terminal nodes (all values are set to zero).
%
%     Notes: 1) all nodes of dead branches must be set to zero.
%            2) the scale before the last one must contain at least one
%               non-terminal node.
%            3) nodes 4, 8, 12, 16, 17..20, ... (correspond to the HH children in
%               classic transform) must be set to zero.
%
%     Use cwpt2_btree function to generate specific btrees.
%
% alf, ahf (double vectors) [optional]
%     Analysis low-pass and high-pass filters
%
% alf_off, ahf_off (double scalars) [optional]
%     Corresponding offsets of the center coefficients (>=0, C-like notation)
%
% If these arguments are omitted, 7-9 filters are used.
%
% transformed_image (M x N x L real double N-D matrix)
%     The transformed full-resolution subbands. The order of subbands is defined
%     by the depth left-to-right search of the transform tree. Subbands are
%     added to the output N-D matrix as soon as they are computed, so leaves
%     of lower scales are always output earlier than leaves of higher scales
%     (which are children of the processed node). "Left-to-right" means
%     LL, H0, 0H order where the first filtering direction is the first
%     memory raster direction of data. For instance, 0H means for Matlab images
%     no filtering along columns (first dimension) and high-filtering along
%     rows (second dimension).
%
%     Use cwpt2_get_leaf_order function to get the explicit order of the subbands.
%
% See also:
%     c3wpti, cwpt2, cwpt2_btree, cwpt2_get_leaf_order
%

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
