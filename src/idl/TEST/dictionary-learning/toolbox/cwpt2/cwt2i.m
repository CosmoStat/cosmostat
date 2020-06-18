function reconstructed = cwt2i(transformed, cubic_spline_filters, first_scale_shifted)
% function reconstructed = cwt2i(transformed, cubic_spline_filters, first_scale_shifted)
%
% Inverse 2D Continuous Wavelet Transform.
%
% Given a M x N x ( 1 + 3 * nb_scales ) N-D matrix of a
% 2D continuous wavelet transform provided by the cwt2 function,
% provides a reconstructed M x N image.
%
% Use cubic_spline_filters = 0 [optional] for 7-9 filters
% Use cubic_spline_filters = 1 [optional] for cubic spline filters
%
% Use first_scale_shifted = 1 [optional] if you wish to consider a transform of
% scale N as a transform of scale (N+1)
%
% This function is based on cwpt2i.
%

if nargin < 2
    cubic_spline_filters = 0;
end;

if nargin < 3
    first_scale_shifted = 0;
end;

if length(size(transformed))~=3
    error('The input argument must be a M x N x ( 1 + 3 * nb_scales ) N-D matrix.');
end;

nb_scales = size(transformed,3) - 1;

if mod(nb_scales,3)~=0
    error('Invalid data: the size along the dimension 3 must be ( 1 + 3 * nb_scales ).');
end;

nb_scales = first_scale_shifted + (nb_scales / 3);
btree = cwpt2_btree(nb_scales,2);       % 2 = wavelet

if cubic_spline_filters
    reconstructed = cwpt2i(transformed, btree, -1, -1, -1, -1);
else
    reconstructed = cwpt2i(transformed, btree);
end;

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
