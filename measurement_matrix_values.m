function A = measurement_matrix_values(mu, SIGMA, x, nnidx) %#eml

[n dim] = size(x);
p_nn = size(nnidx,1);

A = zeros(n,p_nn);

for i = 1:n
    neighborMuIndices = nnidx(:, i);
    
    j_ = 0; % Index into neighborMuIndices for the neighbor we're
            % currently considering.
    
    for j=neighborMuIndices'
        j_ = j_ + 1;
        
        % The i'th measurement point might have less than size(nnidx, 1)
        % neighbors. Spare entries in nnidx will be 0. Don't perform any
        % computation for those.
        if j > 0
            kernel_value = gauss(x(i, :), mu(j, :), reshape(SIGMA(:,j,:,:), [dim dim]));
            A(i, j_) = kernel_value;
        end
    end
end
