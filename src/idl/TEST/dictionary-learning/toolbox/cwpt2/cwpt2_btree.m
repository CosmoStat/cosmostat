function btree = cwpt2_btree(nb_scales, btree_type)
% function btree = cwpt2_btree(nb_scales, btree_type)
%
% Creates a specific btree structure ready to be used.
%
% nb_scales (real double scalar)
%     The number of scales : 1, 2, 3, ... 10
%
% btree_type (real double scalar)
%     1 = full decomposition
%     2 = classic or 3-way wavelet decomposition
%
% btree (B x 1 uint8 vector)
%     The btree structure ready to be used by cwpt2 and cwpt2i functions.
%
% See also:
%     cwpt2
%

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
