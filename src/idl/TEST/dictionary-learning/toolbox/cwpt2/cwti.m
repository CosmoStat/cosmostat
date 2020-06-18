function reconstructed_signal = cwti(transformed_signal)
% function reconstructed_signal = cwti(transformed_signal)
%
% Inverse 1D Continuous Wavelet Transform.
%
% transformed_signal (N x (nb_scales + 1) real double matrix)
%     The full-resolution subbands of the transformed signal obtained
%     by the cwt function. The order of the subbands is:
%     high(1), high(2), ... high(nb_scales), low(nb_scales)
%
% reconstructed_signal (N x 1 real double vector)
%     The reconstructed signal.
%
% See also:
%     cwt
%

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
