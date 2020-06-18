function transformed_signal = cwt(original_signal, nb_scales)
% function transformed_signal = cwt(original_signal, nb_scales)
%
% Forward 1D Continuous Wavelet Transform.
%
% original_signal (N x 1 or 1 x N real double vector)
%     The signal to be transformed.
%
% nb_scales (real double scalar)
%     The number of scales : 1, 2, 3, ... 10
%
% transformed_signal (N x (nb_scales + 1) real double matrix)
%     The full-resolution subbands of the transformed signal. The order
%     of the subbands is:
%     high(1), high(2), ... high(nb_scales), low(nb_scales)
%
% See also:
%     cwti, plotzoom
%

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
