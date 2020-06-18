function f = btree_freqzone(node)
% function f = btree_freqzone(node)
%
% Gets the frequency zone corresponding to the given node.
% f = [f_1_start, f1_end, f2_start, f2_end]
% 0 <= f_i <= 1
% f_i = 1 corresponds to the Nyquist frequency.
%
% See also:
%     btree_freqzone_fft,
%     btree_parent, btree_scale, btree_child_type,
%     btree_length, btree_nb_scales
%

if (node < 0)
    error('Invalid node number.');
end;

if (node == 0)
    f = [0 1 0 1];
    return;
end;

f_parent = btree_freqzone(btree_parent(node));
f_parent_mid_1 = (f_parent(1) + f_parent(2)) / 2;
f_parent_mid_2 = (f_parent(3) + f_parent(4)) / 2;

[child_no, filter_1, filter_2] = btree_child_type(node);

if filter_1         % high-pass in direction 1
    f(1) = f_parent_mid_1;
    f(2) = f_parent(2);
else                % low-pass in direction 1
    f(1) = f_parent(1);
    f(2) = f_parent_mid_1;
end;

if filter_2         % high-pass in direction 2
    f(3) = f_parent_mid_2;
    f(4) = f_parent(4);
else                % low-pass in direction 2
    f(3) = f_parent(3);
    f(4) = f_parent_mid_2;
end;
