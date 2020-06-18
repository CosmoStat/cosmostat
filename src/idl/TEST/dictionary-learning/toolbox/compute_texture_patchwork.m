function [M,Id] = compute_texture_patchwork(H)

% compute_texture_patchwork - mix several textures
%
%   M = compute_texture_patchwork(H);
%
%   H can be a cell array of textures or a (n,n,3) matrix.
%   Handles up to 5 textures.
%
%   Copyright (c) 2006 Gabriel Peyr?

if iscell(H)
    H = cat(3,H{:});
end

p = size(H,3); % number of textures
n = size(H,1);

M = H(:,:,1);
Id = ones(n);
if p>1
    M(1:end/2,end/2+1:end) = H(1:end/2,end/2+1:end, mod(1,p)+1);
    Id(1:end/2,end/2+1:end) = 2;
end
if p>2
    M(end/2+1:end,1:end/2) = H(end/2+1:end,1:end/2, mod(2,p)+1);
    Id(end/2+1:end,1:end/2) = 3;
end
if p>3
    M(end/2+1:end,end/2+1:end) = H(end/2+1:end,end/2+1:end, mod(3,p)+1);
    Id(end/2+1:end,end/2+1:end) = 4;
end
if p>4
    x = linspace(-1,1,n);
    [Y,X] = meshgrid(x,x);
    r = 0.5;
    I = find( X.^2 + Y.^2 <= r^2 );
    M0 = H(:,:, mod(4,p)+1);
    M(I) = M0(I);
    Id(I) = 5;
end