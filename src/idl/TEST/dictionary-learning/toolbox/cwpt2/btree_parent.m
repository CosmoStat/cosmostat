function parent = btree_parent(node)
% function parent = btree_parent(node)
%
% Gets the number of the parent node.
%
% See also:
%     btree_scale, btree_child_type, btree_freqzone,
%     btree_length, btree_nb_scales
%

if (node < 0)
    error('Invalid node number.');
end;

if (node == 0)
    error('Root has no parent.');
end;

parent = floor((node - 1) / 4);
