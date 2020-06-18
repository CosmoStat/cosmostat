function nb_scales = btree_nb_scales(btree)
% function nb_scales = btree_nb_scales(btree)
%
% Gets the maximum scale number of the given btree.
%
% See also:
%     btree_parent, btree_child_type, btree_freqzone,
%     btree_length
%

nb_scales = btree_scale(length(btree) - 1);
