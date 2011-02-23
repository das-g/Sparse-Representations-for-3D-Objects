function A = measurement_matrix(mu, normals, SIGMA, x)

old_path = path;
addpath([pwd '/../ann_mwrapper/extension'])

[x_kdtree, pts] = annConstructTree(x');

[n dim] = size(x);
p = size(mu,1);

A = sparse(n,p);

for j = 1:p
    circumcircle_radius = max(eig(reshape(SIGMA(:,j,:,:), [dim dim])));
    neighborXIndices = annFindKnn(x_kdtree, mu(j, :)', ...
                                  circumcircle_radius ^ 2 * 1e1, ... % FIXME: max. number of neighbors should also depend on the density of x
                                  circumcircle_radius, 0);
    % The i'th measurement point might have less than size(nnidx, 1)
    % neighbors. Eliminate empty idx entires
    neighborXIndices = neighborXIndices(neighborXIndices > 0);
    numXNeighbors = size(neighborXIndices,1);
    
    kernel_value = gauss(x(neighborXIndices, :), ...
                         repmat(mu(j, :), [numXNeighbors 1]), ...
                         reshape(SIGMA(:,j,:,:), [dim dim]));
    A(neighborXIndices,j) = dot(repmat(normals(j, :), [numXNeighbors 1]), ...
                            x(neighborXIndices, :) - repmat(mu(j, :), [numXNeighbors 1]), ...
                            2) .* kernel_value;
end

annDeleteTree(x_kdtree, pts)
path(old_path)