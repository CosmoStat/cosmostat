function reconstructed_image = c3wpti(transformed_image, btree, alf, alf_off, ahf, ahf_off, rlf, rlf_off, rhf, rhf_off)
% function reconstructed_image = c3wpti(transformed_image, btree, alf, alf_off, ahf, ahf_off, rlf, rlf_off, rhf, rhf_off)
%
% Inverse 2D 3-way Continuous Wavelet Packet Transform.
%
% transformed_image (M x N x L real double N-D matrix)
%     The transformed full-resolution subbands provided by the cwpt2 function.
%     2 <= M, N <= 16384, L must be coherent with the provided btree. See
%     c3wpt for more details on subband order.
%
% btree (B x 1 or 1 x B uint8 vector)
%     The btree structure used to obtain the transformed image.
%
% alf, ahf (double vectors) [optional]
%     Analysis low-pass and high-pass filters
%
% alf_off, ahf_off (double scalars) [optional]
%     Corresponding offsets of the center coefficients (>=0, C-like notation)
%
% rlf, rhf (double vectors) [optional]
%     Reconstruction low-pass and high-pass filters
%
% rlf_off, rhf_off (double scalars) [optional]
%     Corresponding offsets of the center coefficients (>=0, C-like notation)
%
% If these arguments are omitted, 7-9 filters are used.
%
% reconstructed_image (M x N real double matrix)
%     The reconstructed image.
%
% See also:
%     c3wpt, cwpt2_btree
%

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
