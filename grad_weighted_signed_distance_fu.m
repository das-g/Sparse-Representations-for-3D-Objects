function [ grad_f ] = grad_weighted_signed_distance_fu( mu, normals, SIGMA, x )
%Compute the gradient of the weighted signed distance function (with gauss kernels)
%
%   value = GRAD_WEIGHTED_SIGNED_DISTANCE_FU(mu, normals, SIGMA, x)
%   Compute the gradient of the weighted signed distance function
%   (for gauss kernels at centers mu with covariance matrices SIGMA)
%   at position x
%
%       mu      is a p-by-d matrix where each of the p rows represents the
%               (d-dimensional) position of a center
%
%       normals is a p-by-d matrix where each of the p rows represents the
%               (d-dimensional) directed normal of the zero-isosurface
%               at a center
%
%       SIGMA   is a p-by-d-by-d array where SIGMA(j,:,:) is the d-by-d
%               covariance matrix corresponding to the j-th center. A valid
%               covariance matrix must be positive-semidefinite. There will
%               be no warning when the passed matrices aren't.
%
%       x       n-by-d array where each of the n rows represents a position
%               where the function shall be evaluated
%
% See also: weighted_signed_distance_fu

[p d] = size(mu);
n = size(x,1);

numerator = 0;
denominator = 0; % for normalization

f = weighted_signed_distance_fu(mu, normals, SIGMA, x);

for j=1:p
    Sigma = squeeze(SIGMA(j,:,:));
    kernel_value = gauss(x, repmat(mu(j,:), [n 1]), Sigma);
    distance = x - repmat(mu(j,:), [n 1]);
    numerator = numerator + ((repmat(dot(repmat(normals(j,:), [n 1]), ...
                                distance, ...
                                2) - f, [1 d]) .* (-2 * distance / Sigma')) ...
                             + repmat(normals(j,:), [n 1]))...
                             .* repmat(kernel_value, [1 d]);
    denominator = denominator + kernel_value;
end

grad_f = numerator ./ repmat(denominator, [1 d]);

end
