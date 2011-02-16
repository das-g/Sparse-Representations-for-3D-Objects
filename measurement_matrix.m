function A = measurement_matrix(mu, normals, SIGMA, x)

old_path = path;
addpath([pwd '/../ann_mwrapper'])

tic
nnidx = annquery(mu', x', 10);
toc

path(old_path)

[n dim] = size(x);
p = size(mu,1);

tic
A_i_indices = [];
A_j_indices = [];
A_values = [];

counter = 1;
for i = 1:n
    neighborMuIndices = nnidx(:, i);
    % The i'th measurement point might have less than size(nnidx, 1)
    % neighbors. Eliminate empty idx entires
    neighborMuIndices = neighborMuIndices(neighborMuIndices > 0)';
    A_i_indices = [A_i_indices repmat(i,size(neighborMuIndices))];
    A_j_indices = [A_j_indices neighborMuIndices];
    A_values = [A_values zeros(size(neighborMuIndices))];
    
    denominator = 0; % for normalization
    row_start = counter;
    
    for j=neighborMuIndices
        kernel_value = gauss(x(i, :), mu(j, :), reshape(SIGMA(:,j,:,:), [dim dim]));
        A_values (counter) = dot(normals(j,:), ...
                                 x(i,:) - mu(j,:), ...
                                 2) * kernel_value;
        denominator = denominator + kernel_value;
        counter = counter + 1;
    end
    
    A_values(row_start:end) = A_values(row_start:end) ./ denominator;
end
toc

A = sparse(A_i_indices, double(A_j_indices), A_values, n, p);
toc