% test for dictionary learning
% The learning is done 
%   * for piecewise-constant 1D signals, 
%   * high contrast image patches, 
% so the learned basis should looks like a Wavelet basis.
%
%   Copyright (c) 2007 Gabriel Peyre

path(path, 'toolbox/');
path(path, '../toolbox_image_data/');

test_type = 'signal';
test_type = 'image';

redun = 2;        % redundandy of the dictionary
overtraining = 15;   % over training factor
        
switch test_type
    case 'signal'
        % load piecewise constant signals
        n = 128/2;  
        K = round( n*redun );           % number of atoms
        m = round( overtraining*K );    % number of samples

        % generate random data
        s = 10; % sparsity
        Y = zeros(n,m);
        for i=1:m
            sel = randperm(n); sel = sel(1:s);
            Y(sel,i) = randn(s,1);
        end
        Y = cumsum(Y);
    case 'image'
        w = 10; n = w^2;
        K = round( n*redun );           % number of atoms
        m = round( overtraining*K );    % number of samples
        % load random patches
        name_list = {,'barb', 'lena', 'boat', 'flinstones', 'goldhill', 'peppers'};     % 'barb', 
        nbimg = length(name_list);
        Y = [];
        for i=1:nbimg
            M = load_image(name_list{i});
            % high pass filter to whiten a little the patches
            M = perform_blurring(M,15)-M;
            H = compute_random_patches(M,w,5*m, w);
            H = reshape(H, n,5*m);
            H = H ./ repmat( sqrt(sum(H.^2,1)), n,1 );
            s = std(H); [tmp,I] = sort(s); I = I(end:-1:1);
            Y = [Y, H(:,I(1:m))];
        end
        sel = randperm(size(Y,2)); sel = sel(1:m);
        Y = Y(:,sel);
end
Y = Y ./ repmat( sqrt(sum(Y.^2,1)), n,1 );



% do the learning with several methods
options.K = K;
options.niter_inversion = 20;
options.niter_learning = 30;
lambda = 0.05/20; % mean( sqrt( sum( Y.^2, 1 ) ) ) / 20;
options.thresh_type = 'soft';
options.lambda_min = lambda;
options.lambda_max = lambda*10;
options.thresh_type = 'omp';
options.strict_sparsity = 5;

%% run learning
disp('--> Dictionary learning.');
options.D = [];
options.X = [];
options.options.init_dico = 'input';
options.learning_method = 'mod';
options.use_mex = 0;
[D1,X1,E1] = perform_dictionary_learning(Y,options);
options.learning_method = 'ksvd';
[D2,X2,E2] = perform_dictionary_learning(Y,options);
options.learning_method = 'randomized';
[D3,X3,E3] = perform_dictionary_learning(Y,options);
options.options.init_dico = 'rand';
options.learning_method = 'ksvd';
[D4,X4,E4] = perform_dictionary_learning(Y,options);

plot(1:length(E1), E1, 1:length(E1), E2, 1:length(E1), E3, 1:length(E1), E4 );
legend('MOD', 'KSVD', 'Mixed', 'KSVD+Input random');

% plot some learned basis vector
if strcmp(test_type, 'image')
    nb = [10 14]; ndim = 1;
else
    nb = [4 6]; ndim = 1;
end

clf;
display_dictionnary(D1, X1, nb,ndim );



return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the rate distortion
Dlist = {D1};

niter_coding = 30;
lambda_min = lambda/5;
lambda_max = lambda*4;
err  = zeros(niter_coding, length(Dlist)); 
bits = zeros(niter_coding, length(Dlist)); 
niter_inversion = 5;
thresh_type = 'hard';

err = zeros(niter_coding,length(Dlist));
for ii=1:length(Dlist)

    D = Dlist{ii};
    [a,s,b] = svd(D);
    mu = 1/max(diag(s))^2;
    X = zeros(K,m);
    % sparse inversion
    for i=1:niter_coding
        progressbar(i,niter_coding);
        t = lambda_max - (i-1)/(niter_coding-1)*(lambda_max-lambda_min);
        for j=1:niter_inversion
            X = X + mu * D'*( Y - D*X );
            X = perform_thresholding( X, t*mu, thresh_type);
        end
        % quantize
        [Xv, Xq] = perform_quantization(X, t, 1);
        % compute error
        err(i,ii) = norm( Y-D*Xv, 'fro' ) / (n*m);
        bits(i,ii) = compute_entropy(Xq);
    end

end
plot(bits, -log(err));
xlabel('#bits'); ylabel('-log(err)'); 
legend('MOD', 'K-SVD');
