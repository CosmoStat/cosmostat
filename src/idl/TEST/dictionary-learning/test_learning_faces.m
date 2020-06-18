% test for dictionary learning of faces


redun = 2;        % redundandy of the dictionary
overtraining = 15;   % over training factor
        
w = 10; n = w^2;
K = round( n*redun );           % number of atoms
m = round( overtraining*K );    % number of samples
Y = readpan('images_dir/Y.pan',1); % read image of patches (assumed normalized to a unit l2 norm)

subset = 3000; % put end if all patches are to be included in the learning process
Y = Y(:,1:subset);

% do the learning with several methods
options.K = K;
options.niter_inversion = 20;
options.niter_learning = 30;
lambda = 0.05/20; % mean( sqrt( sum( Y.^2, 1 ) ) ) / 20;
options.lambda_min = lambda;
options.lambda_max = lambda*10;

%% run learning
disp('--> Dictionary learning.');
options.D = [];
options.X = [];
options.options.init_dico = 'input';
options.thresh_type = 'strict';
options.learning_method = 'mod';
options.use_mex = 0;
[D1,X1,E1] = perform_dictionary_learning(Y,options);
options.thresh_type = 'omp';
options.strict_sparsity = 5;
options.learning_method = 'ksvd';
[D2,X2,E2] = perform_dictionary_learning(Y,options);
options.learning_method = 'randomized';
[D3,X3,E3] = perform_dictionary_learning(Y,options);
options.options.init_dico = 'rand';
options.learning_method = 'ksvd';
[D4,X4,E4] = perform_dictionary_learning(Y,options);

figure
plot(1:length(E1), E1, 1:length(E1), E2, 1:length(E1), E3, 1:length(E1), E4 );
legend('MOD', 'KSVD', 'Mixed', 'KSVD+Input random');

nb = [10 14]; ndim = 1;
figure
display_dictionnary(D1, X1, nb,ndim );
figure
display_dictionnary(D2, X2, nb,ndim );
figure
display_dictionnary(D3, X3, nb,ndim );
figure
display_dictionnary(D4, X4, nb,ndim );

tilefigs

