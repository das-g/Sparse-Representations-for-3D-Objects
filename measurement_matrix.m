function A = measurement_matrix(mu, normals, SIGMA, x)

old_path = path;
addpath([pwd '/../ann_mwrapper'])

nnidx = annquery(mu', x', 10);

path(old_path)

A = measurement_matrix_values(mu, normals, SIGMA, x, nnidx);