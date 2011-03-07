function A = measurement_matrix(mu, normals, SIGMA, x)

%% Find nearest neighbors using ann library.
% Save current path, so we can restore it later.
old_path = path;

% Change path, so that ann wrapper files are found.
addpath([pwd '/../ann_mwrapper']) % Full path needed, so prepend working dir

% Query for nearest 10 mu for each x.
nnidx = annquery(mu', x', 10);

% Clean up. (Restore previous path.)
path(old_path)

%% Determine sizes of input arguments
n = size(x, 1);
p = size(mu, 1);
p_nn = size(nnidx, 1);

%% Indices of measurement matrix' non-zero values
A_i = repmat((1:n)', [1 p_nn]);
A_j = double(nnidx');

%% (Non-zero) values of measurement matrix
% The sizes of passed parameters are not known in advance, so we compile
% measurement_matrix_values just before calling it. For large n or p_nn
% this still pays off.
emlmex measurement_matrix_values -eg {mu, normals, SIGMA, x, nnidx}

A_values = measurement_matrix_values(mu, normals, SIGMA, x, nnidx);

%% Compose measurement matrix
A = sparse(A_i(nnidx' > 0), ...
           A_j(nnidx' > 0), ...
           A_values(nnidx' > 0), ...
           n, p);
