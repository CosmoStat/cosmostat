% perform texture synthesis using learned dictionary

% size of input image
n = 128;
% redundancy of the dictionary
redun = 2;
% strict sparsity of the decomposition
sparsity = 4;

path(path, '../images/');

if not(exist('name'))
    names{2} = 'fine_grid';
    % names{3} = 'herringbone';
    names{3} = 'warped_grid';
    names{4} = 'hair';
    names{5} = 'fur';
    names = {'simoncelli7'};
    names = {'reptilskin'};
    names = {'corral'};
else
    names = {name};
end

basestr = [];
for i=1:length(names)
    basestr = [basestr names{i}];
    if i<length(names)
        basestr = [basestr '-'];
    end
end
addstr = ['r' num2string_fixeddigit(redun,2) 's' num2string_fixeddigit(sparsity,2)];
repimg = ['results/synthesis/' basestr '/'];
repeps = ['results/synthesis/' basestr '/eps/'];
repiter = ['results/synthesis/' basestr '/iter' addstr '/'];
if not(exist(repimg))
    mkdir(repimg);
end
if not(exist(repeps))
    mkdir(repeps);
end
if not(exist(repiter))
    mkdir(repiter);
end

H = {}; Y = {};
Ma = {};
nTextures = length(names);
for k=1:nTextures
    name = names{k};
    switch lower(name)
        case 'corral'
            n0 = 130;
        case 'dunes'
            n0 = 200;
        case 'warped_grid'
            n0 = 450;
        case 'wood'
            n0 = 230;
        case 'grass-big'
            n0 = 200;
        case 'dead_leaf'
            n0 = 400;
        case 'fabric'
            n0 = 400;
        case 'simoncelli3'
            n0 = 300;
        case 'simoncelli7'
            n0 = 180;
        case 'text'
            n0 = 200;
        case 'turbulence1'
            n0 = 150;
        otherwise 
            n0 = [];
    end
    A = load_image(name, n0);
    n = min(size(A,1),n);
    Ma{k} = rescale( crop(sum(A,3),n) );
end
[M,Id] = compute_texture_patchwork(Ma);
s = size(Ma{1},3); % number of channel

warning off;
for i=1:nTextures
    imwrite(rescale(Ma{i}), [repimg basestr '-original-' num2str(i) '.png'], 'png');
end
warning on;

%% parameters of the dictionary 
% width of the patch
w = 8;
% dimensionality of the problem
d = s*w^2;
% size of the dictionary
p = round(redun*d);
% oversampling for the training
tr = 10;
% number of samples
m = tr*p;

%% computing samples for training
H = [];
for i=1:nTextures
    H = cat( 4, H, compute_random_patches(Ma{i},w, ceil(m/nTextures), w) );
end
m = size(H,4);
% column mode
Y = reshape( H, d,m );

% parameters for sparse coding
lambda = 0.03;
lambda = 0.35;
options.lambda_min = lambda;
options.lambda_max = lambda*4;
options.thresh_type = 'hard';
options.thresh_type = 'strict'; 
options.thresh_type = 'soft';
options.thresh_type = 'omp';
options.strict_sparsity = sparsity;
options.use_mex = 0;
% parameters for learning
options.K = p;
options.niter_inversion = 5;    % only used if thresh_type is not 'omp'
options.niter_learning = 20;
options.learning_method = 'ksvd';   % 'ksvd' or 'mod'
options.use_mex = 0;                % for OMP implementation (fast mex sometimes crashes, need to be fixed)

if strcmp(options.thresh_type, 'strict')
    options.mu_dampling = 1.5; % 4; % speed up in case of agressive sparsity (this is experimental)
end

use_wavelet_histo = 1;
n1 = n;
Ms = randn(n1);
niter = 15;

%% run learning
disp('--> Dictionary learning.');
options.D = [];
options.X = [];
[D1,X1,E1] = perform_dictionary_learning(Y,options);


if 0
options.thresh_type = 'omp';
[D1,X1,E1] = perform_dictionary_learning(Y,options);
end

%% display most representatives vectors
nb = [8 12];
A = display_dictionnary(D, X); %, nb, 2, s);
imageplot(A); 

warning off;
imwrite(rescale(A), [repimg basestr '-dictionary-w' num2str(w) '.png'], 'png');
warning on;

% samples to get the states
impose_statistics = 0;
if impose_statistics
    H0 = compute_random_patches(M,w, min( round(0.5*n^2), 10000), w);
    Y0 = reshape(H0, d, size(H0,4));
    options.niter = 30;
    options.X = [];
    fprintf('Samples coding: ');
    [X0,err] = perform_iterative_thresholding(D,Y0,options);
end

% perform synthesis
Ms = randn(n1);
Xs = [];
for i=1:niter
    disp(['--> iter ' num2str(i) '/' num2str(niter) '.']);
    % first histogram matching
    if use_wavelet_histo
        options.dotransform = 1; Jmin = 4;
        Ms = perform_histogram_matching_wavelet(Ms,M,Jmin,options);
    else
        % Ms = perform_histogram_matching(Ms, M, options);
        Ms = perform_histogram_equalization(Ms, Ma{1});
    end
    options.sub = 1;
    options.epsilon = floor( rand(2,1)*w );
    Hs = compute_all_patch(Ms,w, options);
    Ys = reshape( Hs, d,size(Hs,4) );    
    % sparse coding
    options.niter = max(6,20-3*i);
    options.X = Xs;
    [Xs,err] = perform_iterative_thresholding(D,Ys,options);
    % match histograms
    if impose_statistics
        for j=1:p
            Xs(j,:) = perform_histogram_equalization(Xs(j,:), X0(j,:));
        end
    end
    % reconstruct
    Ys = D*Xs;
    Hs = reshape( Ys, [w,w,s,size(Hs,4)] );
    Ms = compute_all_patch(Hs, n1, options);
    % save
    warning off;
    imwrite(rescale(Ms), [repiter basestr '-iter-' num2string_fixeddigit(i,2) '.png'], 'png');
    warning on;
end


warning off;
imwrite(rescale(Ms), [repimg basestr addstr '.png'], 'png');
warning on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return;



% do progressive learning
lambda_min = 0.05;
lambda_max = 0.5;
options.D = [];
options.X = [];
Xs = []; Xd = []; X0 = [];
niter = 15;
for i=1:niter
    disp(['--> iter ' num2str(i) '/' num2str(niter) '.']);
    lambda = lambda_max - (lambda_max-lambda_min)*(i-1)/(niter-1);
    options.lambda_min = lambda;
    options.lambda_max = lambda;
    
    % learn dictionary
    options.niter_learning = max(2, 30-4*i);
    options.niter_inversion = max(5, 30-3*i);
    options.X = Xd;
    fprintf('Dictionary learning: ');
    [D,Xd,E] = perform_dictionary_learning(Y,options);

    % sparse code some patches
    Y0 = reshape(H0, d, size(H0,4));
    options.niter = max(5,30-2*i);
    options.X = X0;
    fprintf('Samples coding:      ');
    [X0,err] = perform_iterative_thresholding(D,Y0,options);
    
    niter_synth = max(1,5-i);
    for j=1:niter_synth

        % first histogram matching
        if use_wavelet_histo
            options.dotransform = 1; Jmin = 4;
            Ms = perform_histogram_matching_wavelet(Ms,M,Jmin,options);
        else
            % Ms = perform_histogram_matching(Ms, M, options);
            Ms = perform_histogram_equalization(Ms, Ma{1});
        end
        Hs = compute_all_patch(Ms,w, options);
        Ys = reshape( Hs, d,size(Hs,4) );

        % sparse coding
        options.niter = max(5, 30-2*i);
        options.X = Xs;
        fprintf('Coefficients coding: ');
        [Xs,err] = perform_iterative_thresholding(D,Ys,options);
        % match histograms
        for j=1:p
            Xs(j,:) = perform_histogram_equalization(Xs(j,:), X0(j,:));
        end
        % reconstruct
        Ys = D*Xs;
        Hs = reshape( Ys, [w,w,s,size(Hs,4)] );
        Ms = compute_all_patch(Hs, n1, options);

    end
    % save
    warning off;
    imwrite(rescale(Ms), [repiter basestr '-iter-' num2string_fixeddigit(i,2) '.png'], 'png');
    warning on;
end


warning off;
imwrite(rescale(Ms), [repimg basestr addstr '-synth.png'], 'png');
warning on;
