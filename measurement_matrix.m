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

% The sizes of passed parameters are not known in advance, so we compile
% measurement_matrix_values just before calling it. For large n or p_nn
% this still pays off.
emlmex measurement_matrix_values -eg {mu, normals, SIGMA, x, nnidx}

A_values = measurement_matrix_values(mu, normals, SIGMA, x, nnidx);

A = sparse(A_i(nnidx' > 0), ...
           A_j(nnidx' > 0), ...
           A_values(nnidx' > 0), ...
           n, p);