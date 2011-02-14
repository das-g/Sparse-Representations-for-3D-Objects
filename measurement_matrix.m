function A = measurement_matrix(mu, normals, SIGMA, x)

old_path = path;
addpath([pwd '/../ann_mwrapper'])

nnidx = annquery(mu', x', 10);

path(old_path)

[n dim] = size(x);
p = size(mu,1);

A = sparse(n,p);

for i = 1:n
    neighborMuIndices = nnidx(:, i);
    % The i'th measurement point might have less than size(nnidx, 1)
    % neighbors. Eliminate empty idx entires
    neighborMuIndices = neighborMuIndices(neighborMuIndices > 0);
    
    denominator = 0; % for normalization
    
    for j=neighborMuIndices'
        kernel_value = gauss(x(i, :), mu(j, :), reshape(SIGMA(:,j,:,:), [dim dim]));
        A(i,j) = dot(normals(j,:), ...
                                x(i,:) - mu(j,:), ...
                                2) * kernel_value;
        denominator = denominator + kernel_value;
    end
    A(i, neighborMuIndices) = A(i, neighborMuIndices) ./ denominator;
end