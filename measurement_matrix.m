function A = measurement_matrix(mu, normals, SIGMA, x)

old_path = path;
addpath([pwd '/../ann_mwrapper'])

nnidx = annquery(mu', x', 10);

path(old_path)

n = size(x, 1);
p = size(mu, 1);
p_nn = size(nnidx, 1);

A_i = repmat((1:n)', [1 p_nn]);
A_j = double(nnidx');
A_values = measurement_matrix_values(mu, normals, SIGMA, x, nnidx);

A = sparse(A_i(nnidx' > 0), ...
           A_j(nnidx' > 0), ...
           A_values(:), ...
           n, p);