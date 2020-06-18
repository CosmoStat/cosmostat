% test for MCA decomposition in a learned wavelet basis

path(path, 'toolbox/');
path(path, 'images/');

name = 'multitextures';
name = 'barb';
name = 'zebra';
name = 'jupe';
name = 'leopard';
name = 'plaids';
name = 'bag';
name = 'chemise';
name = 'wood2';
name = 'barb';
name = 'boat';

%% global options for the tests
use_transparency = 1;
use_wav = 1;
use_dct = 1;
use_dico = 1;
do_inpainting = 1;
mask_type = 'chessboard';
mask_type = 'user';
mask_type = 'rand';

if use_dico
    if use_transparency
        add_str = '-transpar';
    else
        add_str = '-locdico';
    end
else
    add_str = '-nodico';
end
if do_inpainting
    add_str = [add_str '-inpaint' mask_type];
end
repimg = ['results/learned/' name add_str '/'];
if not(exist(repimg))
    mkdir(repimg);
end


n = 256;
n = 128;


n0 = 180;
n0 = 400;
n0 = [];

M0 = load_image(name,n0);
M0 = sum(M0,3);
n = min(n,size(M0,1));
M0 = rescale( crop(M0,n) );
s = size(M0,3);

if use_transparency
    % load the texture
    namet = 'hexholes';
    namet = 'reptilskin';
    namet = 'crochetcropti';
    namet = 'herringbone';
    namet = 'simoncelli2';
    namet = 'hair';
    n0 = n; n0 = []; n0 = n; 
    if strcmp(namet,'hexholes')
        n0 = 200;
    end
    if strcmp(namet,'crochetcropti')
        n0 = 256;
    end
    MT = load_image(namet,n0);
    MT = sum(MT,3);
    MT = rescale( crop(MT,n) );
    % suppress low frequencies
    h = ones(11); h = h / sum(h(:));
    MT = MT - perform_convolution(MT,h);
    warning on;
    % add texture
    M = M0+MT*0.4;
else
    M = M0;
end


warning off;
if use_transparency
    imwrite(rescale(M0), [repimg name '-background.png'], 'png');
    imwrite(rescale(MT), [repimg name '-texture.png'], 'png');
end
imwrite(rescale(M), [repimg name '-original.png'], 'png');
warning on;


%% parameters of the dictionary 
redun = 1.4;  % redundancy, should be larger than 1
% width of the patch
w = 16;  % should be arround 10
% dimensionality of the problem
d = s*w^2;
% size of the dictionary
p = round(redun*d);
% oversampling for the training
tr = 12;     % should be arround 6 or more
% number of samples
m = tr*p;


%% inpainting mask computation
if do_inpainting
    if not(exist('U'))
        switch mask_type
            case 'user'
                options.r = 2;
                options.mode = 'line';
                U = grab_inpainting_mask(M, options);
            case 'chessboard'
                a = 5; % width of the holes
                [Y,X] = meshgrid(0:n-1,0:n-1);
                V = mod( floor(X/a)+floor(Y/a), 2);
                U = M; U(V==0) = Inf;
            case 'rand'
                a = 65; % percent of removed pixels
                sel = randperm(n^2); 
                sel = sel(1:round(a/100*n^2));
                U = M; U(sel) = Inf;
            otherwise
                error('Unknown mask');
        end
    end
else
    U = [];
end


%% creates the set of patches and learn the basis
k = 0; centers = {};
while true && use_dico
    
    if use_transparency && k==0
        % select from the image
        H1 = compute_random_patches(MT,w, m, w);
        k = k+1;
    else
        % use picking to select
        clf;
        hold on;
        imagesc(M0); axis ij; axis image; axis off;
        colormap gray(256);
        [y,x,b] = ginput(1);
        hold off;
        if not(b==1)
            break;
        end
        k = k+1;
        centers{end+1} = [x;y];

        % compute the patches arround picked points
        r0 = 10;
        [Yc,Xc] = meshgrid(1:n,1:n);
        I = find( (Xc-x).^2+(Yc-y).^2 <= r0^2 ); % center of patches
        pc = length(I);
        [dY,dX] = meshgrid(-w/2+1:w/2,-w/2+1:w/2);
        J = repmat(reshape(I(:),[1 1 1 pc]),[w w 1 1]) + repmat(dX,[1 1 1 pc]) + repmat(dY,[1 1 1 pc])*n;
        H0 = M0(J); Y0 = reshape( H0, [s*w^2, size(H0,4)] );

        % compute close random patches
        m0 = m*10;
        H1 = compute_random_patches(M0,w, m0, w);
        Y1 = reshape( H1, [s*w^2, m0] );
        dd = compute_distance_to_points(Y0,Y1);
        [dd,Ind] = min(dd,[],2);
        [tmp,I] = sort(dd); I = I(1:m);
        H1 = H1(:,:,:,I);
    end
    
    % remove the mean
    H1 = H1 - repmat( mean(mean(H1)), [w w s 1] );
    
    % output the patches
    H{k} = H1;
    Y{k} = reshape( H{k}, [s*w^2, m] );
    
    if use_transparency
        break;
    end
end

%% ksvd options
options.nbr_max_atoms = 5; % sparsity, 10 is fine
options.tol = 1e-3;
options.niter_ksvd = 30;    % 12 is fine
options.use_omp = 1;
options.enforce_orthogonality = 0;

%% learn using KSVD
components = {};
% options for transforms
opt.w = w;      % size of patches
opt.Jmin = 4;   % minimum scale for wavelets
opt.n = n;      % size of the image
opt.threshold_factor = 1;
% options for the learning
options.K = p;
options.niter_inversion = 2;
options.niter_learning = 40;
options.thresh_type = 'soft';
options.init_dico = 'rand';
lambda = mean( sqrt( sum( Y{1}.^2, 1 ) ) ) / 20;
options.lambda_min = lambda*0;
options.lambda_max = lambda*2;
D = {};
for i=1:k
    % learn the dictionary
    if exist('save_dico') && length(save_dico)>=i
        D{i} = save_dico{i};
    else
        % [D{i},X] = perform_ksvd(Y{i},p,options);
        [D{i},X,E] = perform_dictionary_learning(Y{i},options);
        save_dico{i} = D{i};
    end
    % set the options for the transform
    opt.D = D{i};
    opt.threshold_factor = 1; % reduce influency
    % set up the component
    clear cpt;
    cpt.options = opt;
    cpt.callback =  @callback_dictionary;
    cpt.name = ['dico' num2str(i)];
    % add to the list of components
    components{end+1} = cpt;
    
    %% display most representatives vectors
    nb0 = [8,12];
    nb = ceil(sqrt(p)*0.7); nb = [nb, ceil(p/nb)];
    nb = min(nb,nb0);
    A = display_dictionnary(D{i},X,nb, s);
    clf;
    imageplot(A);
    warning off;
    imwrite(rescale(A), [repimg name '-dictionary-' num2str(i) '.png'], 'png');
    warning on;
end
%% classical dictionaries
if use_wav
    opt.threshold_factor = 1; % reduce influency
    clear cpt;
    cpt.options = opt;
    cpt.callback =  @callback_atrou;
    cpt.name = 'wav';
    cpt.tv_correction = 1; % add TV penalty
    components{end+1} = cpt;
end
if use_dct
    opt.w = 32;      % size of patches
    opt.q = opt.w/4; % controls overlap
    opt.threshold_factor = 1; % reduce influency if <1
    clear cpt;
    opt.dct_type = 'orthogonal4';
    opt.dct_type = 'redundant';
    opt.remove_lowfreq = 1; % force to 0 the lowfrequency
    cpt.options = opt;
    cpt.callback =  @callback_localdct;
    cpt.name = 'dct';
    components{end+1} = cpt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stability analysis
if 0
Lambda = [];
for i=1:k
    [U,lambda,V] = svd(D{i});
    Lambda = [Lambda, diag(lambda)];
    % inner product matrix
    Q = abs( D{i}'*D{i} );
end
if k>0
    clf;
    plot(log2(Lambda));
    axis([1 size(Lambda,1) -10 max(log2(Lambda(:)))]);
    legend({dico_names{1:k}});
    saveas(gcf, [repimg name '-stability.png'], 'png');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perform the decomposition
% options for MCA
options.niter = 80;
options.Tmax = 2.7;
options.Tmin = 0;
options.InpaintingMask = U;
options.saverep = [repimg 'iter/'];
if not(exist(options.saverep))
    mkdir(options.saverep);
end
Minit = M;
if do_inpainting
    Minit(U==Inf) = 0;
end
% perform the mca
ML = perform_mca(Minit, components, options);
% number of components
s = size(ML,3);
% recovered image
M1 = sum(ML,3);


clf;
if do_inpainting
    s1 = 2;
    s2 = max(3,s);
    imageplot(M, 'Ground trust', s1,s2,1);
    subplot(s1,s2,2);
    hold on;
    imagesc(Minit); axis ij; axis image; axis off; title('Original');
    for i=1:length(centers)
        hh = plot(centers{i}(2), centers{i}(1), 'r.');
        set(hh, 'MarkerSize', 20);
    end
    hold off;
    imageplot(M1, 'Inpainted', s1,s2,3); 
    warning off;
    imwrite(rescale(Minit), [repimg name '-toinpaint.png'], 'png');
    imwrite(clamp(M1), [repimg name '-inpainted.png'], 'png');
    warning on;
else
    s1 = 1;
    s2 = s+1;
    subplot(s1,s2,2); 
    hold on;
    imagesc(M); axis ij; axis image; axis off; title('Original');
    for i=1:length(centers)
        hh = plot(centers{i}(2), centers{i}(1), 'r.');
        set(hh, 'MarkerSize', 20);
    end
    hold off;
end
for i=1:s
    A = ML(:,:,i);
    cpt = components{i};
    if strcmp(cpt.name, 'wav')
        A = clamp(A); % ensure same dynamic range for cartoon layer
    end
    imageplot(ML(:,:,i), ['Layer ' cpt.name], s1,s2,i+s2);
    warning off;
    imwrite(rescale(A), [repimg name '-' cpt.name '.png'], 'png');
    warning on;
end
saveas(gcf, [repimg name add_str '-results.png'], 'png');


return;

clf;
ax = [];
ax(1) = subplot(1,s+1,1);
hold on;
imagesc(M); axis ij; axis image; axis off;
for i=1:length(centers)
    hh = plot(centers{i}(2), centers{i}(1), 'r.');
    set(hh, 'MarkerSize', 20);
end
hold off;
title('Original');
