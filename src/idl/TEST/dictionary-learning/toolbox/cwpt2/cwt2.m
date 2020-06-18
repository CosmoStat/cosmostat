function transformed = cwt2(signal, nb_scales, cubic_spline_filters)
% function transformed = cwt2(signal, nb_scales, cubic_spline_filters)
%
% 2D Continuous Wavelet Transform.
%
% Given a M x N image, provides a M x N x ( 1 + 3 * nb_scales )
% N-D matrix. The order of results is:
%
% HL, LH, HH (scale 1),
% HL, LH, HH (scale 2),
% ...
% LL, HL, LH, HH (scale nb_scales),
%
% where the first filtering direction is the first memory raster direction of data.
% For instance, LH means for Matlab images low-filtering along columns (first dimension)
% and high-filtering along rows (second dimension).
%
% Use cubic_spline_filters = 0 [optional] for 7-9 filters
% Use cubic_spline_filters = 1 [optional] for cubic spline filters
%
% This function is based on cwpt2.
%

if nargin < 3
    cubic_spline_filters = 0;
end;

btree = cwpt2_btree(nb_scales,2);          % 2 = wavelet

if cubic_spline_filters
    transformed = cwpt2(signal, btree, -1, -1, -1, -1);
else
    transformed = cwpt2(signal, btree);
end;

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
