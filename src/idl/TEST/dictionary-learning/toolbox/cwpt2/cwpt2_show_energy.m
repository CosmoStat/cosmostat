function cwpt2_show_energy(transformed, btree)
% function cwpt2_show_energy(transformed, btree)
%
% Given a signal transformed by cwpt2 and an associated btree,
% shows its energy distribution over the subbands
%

leaf_order = double(cwpt2_get_leaf_order(btree, 4));        % 4 = quadtree
nb_leaves = size(transformed,3);

if nb_leaves ~= length(leaf_order)
    error('Invalid btree or transformed data.');
end;

nb_scales = btree_nb_scales(btree);
n = 2 ^ nb_scales;
nrj_image = zeros(n, n);

for nb_leaf = 1:nb_leaves               % it is better to make a loop because dt data may be huge
    node = leaf_order(nb_leaf);

    f = btree_freqzone(node);
    f = f * n;

    nrj_image((f(1)+1):f(2),(f(3)+1):f(4)) = sum(sum(transformed(:,:,nb_leaf).^2));

end;

figure(2);
imagesc(nrj_image);
colormap gray;
axis square;
axis off;
