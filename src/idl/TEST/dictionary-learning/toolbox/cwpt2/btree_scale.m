function scale = btree_scale(node)
% function scale = btree_scale(node)
%
% Gets the scale of the node.
%
% See also:
%     btree_parent, btree_child_type, btree_freqzone,
%     btree_length, btree_nb_scales
%

if (node < 0)
    error('Invalid node number.');
end;

scale_start = 1;
for scale = 1:10;
    if (node<scale_start)
        scale = scale - 1;
        return;
    end;
    scale_start = scale_start * 4 + 1;
end;

error('Scale is larger than 10');
