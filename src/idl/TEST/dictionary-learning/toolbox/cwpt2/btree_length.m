function len = btree_length(nb_scales)
% function len = btree_length(nb_scales)
%
% Gets the required length for a btree of the given maximum scale.
% Note: its elements are enumerated 0 .. (length - 1)
%
% See also:
%     btree_parent, btree_child_type, btree_freqzone,
%     btree_nb_scales
%

if (nb_scales < 1 || nb_scales > 10)
    error('Invalid number of scales.');
end;

scale_start = 5;
for scale = 2:11;
    if (nb_scales == (scale - 1))
        len = scale_start;
        return;
    end;
    scale_start = scale_start * 4 + 1;
end;

error('This error NEVER occurs');
