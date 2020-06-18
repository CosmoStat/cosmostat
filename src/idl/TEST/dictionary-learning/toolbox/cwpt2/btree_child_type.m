function [child_no, filter_1, filter_2] = btree_child_type(node)
% function [child_no, filter_1, filter_2] = btree_child_type(node)
%
% Gets the "child number" of the node:
%   1 = LL
%   2 = HL
%   3 = LH
%   4 = HH
%
% and its filtering type along the first and the second memory
% raster directions: 0 = low-pass; 1 = high_pass
%
% See also:
%     btree_parent, btree_scale, btree_freqzone,
%     btree_length, btree_nb_scales
%

if (node < 0)
    error('Invalid node number.');
end;

if (node == 0)
    error('Root is not a child.');
end;

child_no = mod(node, 4);
switch child_no
case 1
    filter_1 = 0;
    filter_2 = 0;
case 2
    filter_1 = 1;
    filter_2 = 0;
case 3
    filter_1 = 0;
    filter_2 = 1;
case 0
    child_no = 4;
    filter_1 = 1;
    filter_2 = 1;
end
