function [ result ] = gauss_mix_eval( mu, SIGMA, x )
% Evaluate the unweighted gaussian mixture probability density function
%
%   Evaluate the gaussian mixture probability density function given by the
%   centers mu and the corresponding covariance matrices SIGMA at
%   position(s) x
%
%       mu      is a p-by-d matrix where each of the p rows represents the
%               (d-dimensional) position of a center
%
%       SIGMA   is a p-by-d-by-d array where SIGMA(j,:,:) is the d-by-d
%               covariance matrix corresponding to the j-th center
%
%       x       n-by-d array where each of the n rows represents a position
%               where the function shall be evaluated
%
%   See also: gauss, plot_gauss_mix

p = size(mu,1);
n = size(x,1);

result = 0;
for j=1:p
    result = result + gauss(x, repmat(mu(j,:), [n 1]), squeeze(SIGMA(j,:,:)));
end

end
