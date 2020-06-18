function [D,X] = perform_ksvd(Y,K,options)

% perform_ksvd - learn a dictionnary using K-SVD algorithm
%
%   D = perform_ksvd(Y,K,options)
%
%   Y is a matrix of size (n,m) containing m column examplar
%       vector in R^n.
%
%   D is a dictionnary matrix of size (n,K) of K vectors in R^n
%       that should approximate well the vectors of Y 
%       with few components.
%
%   The number of iteration of the algorithm is given in
%       options.niter_ksvd.
%
%   The algorithm uses Orthogonal Matching Pursuit coding if 
%       options.use_omp=1, otherwise it used sparse inversion.
%   The number of iteration for the sparse inversion is given 
%       options.niter_inversion
%       (see perform_iterative_thresholding).
%
%   The algorithm is described in
%       K-SVD: An Algorithm for Designing Overcomplete Dictionaries for Sparse Representation
%       Michal Aharon Michael Elad Alfred Bruckstein
%       2006
%
%   Copyright (c) 2006 Gabriel Peyr?


options.null = 0;
if isfield(options, 'niter_ksvd')
    niter = options.niter_ksvd;
else
    niter = 1;
end
if isfield(options, 'niter_inversion')
    niter_inversion = options.niter_inversion;
else
    niter_inversion = 20;
end
if isfield(options, 'enforce_orthogonality')
    enforce_orthogonality = options.enforce_orthogonality;
else
    enforce_orthogonality = 0;
end


if niter>1
    options.niter_ksvd = 1;
    T = [];
    fprintf('-> Performing K-SVD ...\n');
    for i=1:niter
        [D,X] = perform_ksvd(Y,K,options);
        options.D = D;
        options.X = X;
    end
    return;
end

[n,m] = size(Y);

if isfield(options, 'D')
    D = options.D;
else
    % D = rand(n,K);
    % intialize with random vectors from input
    sel = randperm(m); sel = sel(1:K);
    D = Y(:,sel);
end
   
if isfield(options, 'X')
    options.x = options.X;
end

  
if isfield(options, 'use_omp')
    use_omp = options.use_omp;
else
    use_omp = 0;
end

% perform sparse coding
if use_omp
    X = perform_omp(D,Y,options);
    X = full(X);
else
    options.niter = niter_inversion;
    [X,T] = perform_iterative_thresholding(D,Y,options);
    
    if 0
    for i=1:size(Y,2)
        maxIters = 10;
        % [X(:,i), iters, activationHist] = SolveMP(D, Y(:,i), K, maxIters);
        [X(:,i), iters, activationHist] = SolveOMP(D, Y(:,i), K, maxIters);
    end
    end
end

T = 1e-3;
% perform SVD update
for k=1:K
    % find the exemplar that are using dictionnary basis k
    I = find( abs(X(k,:))>T );
    if ~isempty(I)
        % compute redisual
        if 0
            D0 = D; D0(:,k) = 0;
            E0 = Y - D0*X;
            % restrict to element actually using k
            E0 = E0(:,I);
        else
            S = X(:,I);
            S(k,:) = 0;
            E = Y(:,I) - D*S;
        end
        % perform SVD
        [U,S,V] = svd(E);
        D(:,k) = U(:,1);
        X(k,I) = S(1) * V(:,1)';
    end
end

if enforce_orthogonality
    % perform the projection on the set of orthgonal matrices
    if K==n
        [U,S,V] = svd(D);
        S = diag( diag(S)./abs(diag(S)) );
        D = U*S*V';
        % recompute projection
        X = D' * Y;
    else
        warning('Orthogonality is enforced only for square matrices.');
    end
end