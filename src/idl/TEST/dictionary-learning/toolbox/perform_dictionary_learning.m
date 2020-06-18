function [D,X,E] = perform_dictionary_learning(Y,options)

% perform_dictionary_learning - learn a dictionnary using K-SVD algorithm
%
%   D = perform_dictionary_learning(Y,options)
%
%   Y is a matrix of size (n,m) containing m column examplar
%       vector in R^n.
%
%   D is a dictionnary matrix of size (n,K) of K vectors in R^n
%       that should approximate well the vectors of Y 
%       with few components.
%       K is given by options.K
%
%   The algorihtm perform the folowing optimization jointly on D and X
%       min_{D,X}  |Y-D*X|^2 + lambda * |X|_1
%           subject to columns of D having unit norm.
%   where |X|_1=\sum |X(i,j)| and lambda is given by options.lambda
%
%   It does this optimisation using block coordinate descent.
%       * Optimization of X knowing D amount to sparse coding, solved using
%       iterative thresholding:
%           X <- Thresh_lambda( X + D'*( Y - D*X ) )
%       This step is repeated options.niter_inversion times.
%       * Optimization of D knowing X amount to L2 best fit:
%           D <- Y*X^+   where   X^+ = X'*(X*X')^{-1}
%
%   This algorithm is very much inspired by the MOD algorithm
%   for dictionary design.
%
%   The number of iteration of the algorithm is given in
%       options.niter_learning.
%   The number of iteration for sparse coding during each step of learning
%       is given by options.niter_inversion.
%
%   options.options.learning_method can be set to:
%       'ksvd': in this case, the dictionary update is computed using the
%           K-SVD algorithm explained in 
%           "K-SVD: An Algorithm for Designing Overcomplete Dictionaries
%                   for Sparse Representation"
%           Michal Aharon Michael Elad Alfred Bruckstein, 2006
%       'mod': in this case, the dictionary update is computed using the L2
%           best fit as proposed by
%           "Method of optimal directions for frame design",
%           K. Engan, S. O. Aase, and J. H. Hus¿y,
%           in Proc. ICASSP, Phoenix, AZ, Mar. 1999, pp. 2443Ð2446. 
%       'randomized': apply randomly one or the other algorithm.
%
%   MOD is usualy faster but KSVD gives a better optimization of the energy.
%
%   Copyright (c) 2007 Gabriel Peyre


[n,m] = size(Y);

options.null = 0;
if isfield(options, 'niter_learning')
    niter = options.niter_learning;
else
    niter = 1;
end
if isfield(options, 'K')
    K = options.K;
else
    K = 1;
end
if isfield(options, 'niter_inversion')
    niter_inversion = options.niter_inversion;
else
    niter_inversion = 4;
end
if isfield(options, 'lambda')
    lambda = options.lambda;
else
    lambda = mean( sqrt( sum( Y.^2, 1 ) ) ) / 50;
end
if isfield(options, 'thresh_type')
    thresh_type = options.thresh_type;
else
    thresh_type = 'soft';
end
if isfield(options, 'mu')
    mu = options.mu;
else
    mu = 1/K;
end
if isfield(options, 'init_dico')
    init_dico = options.init_dico;
else
    init_dico = 'input';
end

if niter>1
    if isfield(options, 'lambda_min')
        lambda_min = options.lambda_min;
    else
        lambda_min = lambda;
    end
    if isfield(options, 'lambda_min')
        lambda_max = options.lambda_max;
    else
        lambda_max = lambda_min*4;
    end
    options.niter_learning = 1;
    E = [];
    for i=1:niter
        options.lambda = lambda_max - (i-1)/(niter-1)*(lambda_max-lambda_min);
        progressbar(i,niter);
        [D,X] = perform_dictionary_learning(Y,options);
        options.D = D;
        options.X = X;
        if strcmp(thresh_type, 'strict') || strcmp(thresh_type, 'omp')
            E(end+1) = norm( Y-D*X, 'fro');
        else
            E(end+1) = 1/2*norm( Y-D*X, 'fro')^2 + lambda_min * sum( abs(X(:)) );
        end
    end
    % sort the basis vector according to energy
    e = sum( X.^2, 2 );
    [tmp,I] = sort(e); I = I(end:-1:1);
    D = D(:,I); X = X(I,:);
    return;
end

if isfield(options, 'D') && not(isempty(options.D))
    D = options.D;
else
    % intialize dictionary
    switch init_dico
        case 'input'
            sel = randperm(m); sel = sel(1:K);
            D = Y(:,sel);
        case 'haar1d'
            D = compute_haar_matrix(n,1);
        case 'haar2d'
            D = compute_haar_matrix([sqrt(n) sqrt(n)],1);
        case 'rand'
            D = randn(n,K);
    end
    if size(D,2)<K
        D = [D; randn(n,K-size(D,2))];
    elseif size(D,2)>K
        D = D(:,1:K);
    end
    D = D - repmat( mean(D), n,1 );
    D(:,1) = 1;
    D = D ./ repmat( sqrt(sum(D.^2,1)), n,1 );
end

if isfield(options, 'mu')
    mu = options.mu;
else
    [a,s,b] = svd(D); 
    mu = 1/max(diag(s))^2;
end


if isfield(options, 'X') && not(isempty(options.X))
    X = options.X;
else
    X = zeros(size(D,2), size(Y,2));
end

if isfield(options, 'strict_sparsity')
    strict_sparsity = options.strict_sparsity;
else
    strict_sparsity = 10;
end

if isfield(options, 'mu_dampling')
    mu_dampling = options.mu_dampling;
else
    mu_dampling = 1;
end

% sparse inversion
E = zeros(niter_inversion,1); % to monitor the enery
if not(strcmp(thresh_type, 'omp'))
    for i=1:niter_inversion
        X = X + mu * mu_dampling * D'*( Y - D*X );
        if strcmp(thresh_type, 'strict')
            t = strict_sparsity;
        elseif strcmp(thresh_type, 'soft')
            t = lambda*mu_dampling*mu;
        elseif strcmp(thresh_type, 'hard')
            % for hard thresholding, the scaling is different
            t = lambda*sqrt(mu_dampling*mu);
        end
        X = perform_thresholding( X, t, thresh_type);
        if strcmp(thresh_type, 'strict')
            E(i) = norm( Y-D*X, 'fro');
        else
            E(i) = 1/2*norm( Y-D*X, 'fro')^2 + lambda * sum( abs(X(:)) );
        end 
    end
else
    if not(isfield(options, 'use_mex'))
        options.use_mex = 1;
    end
    options.nbr_max_atoms = strict_sparsity;
    options.verb = 0;
    X = perform_omp(D, Y, options);
    if issparse(X)
        X = full(X);
    end
end

if E(end)>E(1)
    warning('Iterative thresholding did not converge, you should lower options.mu_dampling.');
end

if isfield(options, 'learning_method')
    learning_method = options.learning_method;
else
    learning_method = 'mod';
end

if strcmp(learning_method, 'randomized')
    if isfield(options, 'randomized_ksvd_proportion')
        randomized_ksvd_proportion = options.randomized_ksvd_proportion;
    else
        randomized_ksvd_proportion = 0.3;
    end
    if rand<randomized_ksvd_proportion
        learning_method = 'ksvd';
    else
        learning_method = 'mod';
    end
end

if strcmp(learning_method, 'ksvd')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       KSVD algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = 1e-3;
    % the first atoms is supposed to be constant
    for k=2:K
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
elseif strcmp(learning_method, 'mod')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       MOD algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dictionary fit
    warning off;
    % D = Y * X'*(X*X')^(-1);
    D = Y * pinv(X);
    warning on;
    % Normalize
    d = sqrt(sum(D.^2,1));
    I = find(d<1e-9); d(I) = 1;
    D(:,I) = randn(size(D,1),length(I));
    d = sqrt(sum(D.^2,1));
    D = D ./ repmat( d, n,1 );
    % [a,b,s] = svd(D*sqrt(mu)); plot(diag(b));    
else
    error('Unknown learning method');
end
        
    
    
    
