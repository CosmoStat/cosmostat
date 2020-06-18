function H = compute_all_patch(M,w, options)

% compute_all_patch - compute all the patch from an image
%
% Extract the patches from an image
%   H = compute_all_patch(M,w, options);
% Reconstruct an image using averaging
%   M = compute_all_patch(H,n, options);
%
% options.bound is either 'per' or 'stop'
%   for periodic or stop at boundary.
% options.sub controls the subsampling.
%   sub=1 extracts every patch.
%   sub=w extracts non overlapping patches
% options.epsilon can add an initial jitter 
%   to the sampling.
%
%   Copyright (c) 2006 Gabriel PeyrŽ

options.null = 0;

if isfield(options, 'sub')
    sub = options.sub;
else
    sub = 1;
end
if isfield(options, 'bound')
    bound = options.bound;
else
    bound ='per';
end
if isfield(options, 'epsilon')
    epsilon = options.epsilon;
else
    epsilon = [0 0];
end

dir = 1;
if size(M,4)>1
    dir = -1;
end

if dir==1
    n = size(M,1);
else
    n = w;
    w = size(M,1);
end
s = size(M,3);

if strcmp(bound, 'per')
    x = 1:sub:n;
else
    x = 1:sub:n-w;
end

[Y,X] = meshgrid(x, x);
X = X(:) + epsilon(1);
Y = Y(:) + epsilon(2);

% number of patches
p = size(X,1);
if dir==-1 && not(p==size(M,4))
    error('Problem during patch extraction.');
end
    
%%% in this case, a fast sampling can be used %%%
[dY,dX] = meshgrid(0:w-1,0:w-1);
Xp = repmat( reshape(X,[1,1,1,p]) ,[w w s 1]) + repmat(dX,[1 1 s p]);
Yp = repmat( reshape(Y,[1,1,1,p]) ,[w w s 1]) + repmat(dY,[1 1 s p]);
Cp = repmat( reshape(1:s,[1 1 s]), [w w 1 p]);


if strcmp(bound, 'per')
    Xp = mod(Xp-1,n)+1;
    Yp = mod(Yp-1,n)+1;
end

I = sub2ind([n n s], Xp,Yp,Cp);
    
if dir==1
    H = M(I);
else
    W = zeros(n,n,s);
    H = W;
    for i=1:p
        H(I(:,:,:,i)) = H(I(:,:,:,i)) + M(:,:,:,i);
        W(I(:,:,:,i)) = W(I(:,:,:,i)) + 1;
    end
    H = H./W;
end