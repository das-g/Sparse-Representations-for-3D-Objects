function [ f ] = weighted_signed_distance_fu( mu, normals, SIGMA, x )
% Compute the weighted signed distance function (with gauss kernels)
%
%   value = WEIGHTED_SIGNED_DISTANCE_FU(mu, normals, SIGMA, x)
%   Compute the weighted signed distance function (for gauss kernels at
%   centers mu with covariance matrices SIGMA) at position x
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
% See also: grad_weighted_signed_distance_fu, plot_f

p = size(mu,1);
n = size(x,1);

numerator = 0;
denominator = 0; % for normalization

for j=1:p
    kernel_value = gauss(x, repmat(mu(j,:), [n 1]), squeeze(SIGMA(j,:,:)));
    numerator = numerator + dot(repmat(normals(j,:), [n 1]), ...
                                x - repmat(mu(j,:), [n 1]), ...
                                2) .* kernel_value;
    denominator = denominator + kernel_value;
end

f = numerator ./ denominator;

end
