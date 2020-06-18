function sp = btree_freqzone_fft(node, style)
% function sp = btree_freqzone_fft(node, style)
%
% Same as btree_freqzone, but returns the frequency zone as:
%
% 1) an image (style = 1 or omitted). If output parameter is
%    omitted, the image is displayed.
%
% 2) an array of data ready to be passed to ifft2 function (style = 2)
%
% Can be used for generating signals containing frequencies in the
% given frequency zone only or for visualizing the corresponding zones.
%
% See also:
%     btree_freqzone
%

if nargin < 2
    style = 1;
end;

f = btree_freqzone(node);
scale = btree_scale(node);
n = 2 ^ scale;
f = f * n;

sp1 = zeros(n, n);
sp1(f(2),f(4)) = 1;

if style == 1
    sp = sp1;

    if nargout == 1
        return;
    end;

    i_show = chessboard(n,n,0,0.1);
    i_show(find(sp)) = 1;
    imagesc(i_show);
    colormap gray;
    axis square;
    axis off;
    title(sprintf('Node %d', node));

    return;
end;

% add 'nearest' parameter if you wish to widen the spectrum
sp = fpartition2fft(sp1, size(sp1)*4);      % frequency oversampling factor
