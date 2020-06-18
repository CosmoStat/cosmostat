% test for segmentation + learning

path(path, 'toolbox/');
path(path, 'images/');

name = 'zebra';
name = 'jupe';
name = 'plaids';
name = 'chemise';
name = 'wood2';
name = 'leopard';
name = 'multitextures';
name = 'multitextures2';
name = 'barb';
name = 'bag';

%% load the image
n0 = [];
ndico = 2;
if strcmp(name, 'multitextures')
    ndico = 4;
end
if strcmp(name, 'bag')
    ndico = 3;
end
if strcmp(name, 'leopard')
    ndico = 5;
end
if strcmp(name, 'barb')
    ndico = 5;
    n0 = 250;
end
n = 200;
n = 150;
M = crop( load_image(name, n0), n); M = sum(M,3);
M = rescale( M );
s = 1;

%% parameters of the dictionaries 
redun = 0.8;  % redundancy, should be larger than 1
% width of the patches
w = 10;  % should be arround 10
% dimensionality of the problem
d = s*w^2;
% size of the dictionary
p = round(redun*d);
% oversampling for the training
tr = 8;     % minimum overtraining
% number of samples
m = tr*p;


use_picking = 0;
% initialize the dictionaries
for i=1:ndico
    if use_picking
        % use picking to select
        clf;
        hold on;
        imagesc(M); axis ij; axis image; axis off;
        title(['Region to initialize dictionary ' num2str(i) '/' num2str(ndico) '.']);
        colormap gray(256);
        [y,x,b] = ginput(1);
        hold off;

        % compute the patches arround picked points
        r0 = 10;
        [Yc,Xc] = meshgrid(1:n,1:n);
        I = find( (Xc-x).^2+(Yc-y).^2 <= r0^2 ); % center of patches
        pc = length(I);
        [dY,dX] = meshgrid(-w/2+1:w/2,-w/2+1:w/2);
        J = repmat(reshape(I(:),[1 1 1 pc]),[w w 1 1]) + repmat(dX,[1 1 1 pc]) + repmat(dY,[1 1 1 pc])*n;
        H0 = M(J); Y0 = reshape( H0, [s*w^2, size(H0,4)] );

        % compute close random patches
        m0 = p*10;
        H1 = compute_random_patches(M,w, m0, w);
        Y1 = reshape( H1, [s*w^2, m0] );
        dd = compute_distance_to_points(Y0,Y1);
        [dd,Ind] = min(dd,[],2);
        [tmp,I] = sort(dd); I = I(1:p);
        H1 = H1(:,:,:,I);
        D{i} = reshape( H1, d, p );
        D{i} = D{i} - repmat( mean(D{i}), [d 1] );
    else
        % just use a bunch of random patches with high contrast
        m0 = p*10;
        H1 = compute_random_patches(M,w, m0, w);
        Y1 = reshape( H1, [s*w^2, m0] );
        v = std( Y1 ); [tmp,I] = sort(v); I = I(end:-1:1);
        H1 = H1(:,:,:,I(1:p));
        D{i} = reshape( H1, d, p );
        D{i} = D{i} - repmat( mean(D{i}), [d 1] );
    end
    D{i} = D{i} ./ repmat( sqrt(sum(D{i}.^2)), [d 1] );
end

% parameters for sparse coding
lambda = 0.03;
lambda = 0.25;
options.lambda_min = lambda;
options.lambda_max = lambda*6;
options.thresh_type = 'soft';   % use iterative thresholding as L1-penalized sparse solver
options.thresh_type = 'strict'; % use iterative StOMP as sparse solver (usually the best)
options.thresh_type = 'omp';    % use OMP as strict sparse solver
options.strict_sparsity = 5;    % sparsity imposer by strict sparse solver
% parameters for learning
options.K = p;
options.niter_inversion = 10;
options.niter_learning = 25;
options.learning_method = 'ksvd';   % 'ksvd' or 'mod'
options.use_mex = 0;                % for OMP implementation (fast mex sometimes crashes, need to be fixed)

if strcmp(options.thresh_type, 'strict')
    options.mu_dampling = 5; % speed up in case of agressive sparsity (this is experimental)
end

% do dictionary optimization
niter = 6;
options.bound = 'stop';
options.sub = 2;
H = compute_all_patch(M, w, options);
q = size(H, 4); q1 = sqrt(q);
H = reshape(H, d, q);
E = zeros(ndico, q);
Bsvg = {};
for i=1:niter
    disp(['--> Iteration ' num2str(i) '/' num2str(niter) '.']);
    % sparse code the set of pixels
    for k=1:ndico
        options.niter = 60;
        options.X = [];
        fprintf('Sparse coding: ');
        [X{k},err] = perform_iterative_thresholding(D{k},H,options);
        if strcmp(options.thresh_type, 'soft')
            E(k,:) = 1/2 * sum( (D{k}*X{k}-H).^2 ) + options.lambda_min*sum( abs(X{k}) );
        else    % strict sparsity 
            E(k,:) = sum( (D{k}*X{k}-H).^2 );
        end
    end
    % best dictionary
    sigma = 4; E1 = E;
    for k=1:ndico
        B = reshape( E(k,:), q1, q1 );
        B = perform_blurring(B, sigma);
        E1(k,:) = B(:)';
    end
    [tmp,B] = min(E1); B = reshape( B, q1, q1 );
    Bsvg{end+1} = B;
    % learning of the dictionary
    for k=1:ndico
        I = find( B==k );
        if length(I)<m
            % not enough samples ...
            A = zeros(q1); A(I) = 1;
            A = perform_blurring(A,2); 
            J = find( (A>1e-3) & (A<0.999) ); 
            sel = randperm(length(J)); 
            s = min( m-length(I), length(sel) );
            sel = sel(1:s);
            I = [I; J(sel)];
        end
        options.D = D{k};
        options.X = X{k}(:,I);
        fprintf('Dictionary learning: ');
        [D{k},X{k},Etmp] = perform_dictionary_learning(H(:,I),options);
    end
end

