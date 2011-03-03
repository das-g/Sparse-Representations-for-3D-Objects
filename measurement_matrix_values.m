function A = measurement_matrix_values(mu, normals, SIGMA, x, nnidx)

[n dim] = size(x);
p_nn = size(nnidx,1);

A = zeros(n,p_nn);

for i = 1:n
    neighborMuIndices = nnidx(:, i);
    % The i'th measurement point might have less than size(nnidx, 1)
    % neighbors. Eliminate empty idx entires
    neighborMuIndices = neighborMuIndices(neighborMuIndices > 0);
    
    denominator = 0; % for normalization
    
    j_ = 0;
    
    for j=neighborMuIndices'
        j_ = j_ + 1;
        kernel_value = gauss(x(i, :), mu(j, :), reshape(SIGMA(:,j,:,:), [dim dim]));
        A(i, j_) = dot(normals(j,:), ...
                                x(i,:) - mu(j,:), ...
                                2) * kernel_value;
        denominator = denominator + kernel_value;
    end
    A(i, :) = A(i, :) ./ denominator;
end