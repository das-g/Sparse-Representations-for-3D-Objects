function [ result ] = gauss_mix_eval( mu, SIGMA, x )
%Evaluate the gaussian mixture probability density function
%   Evaluate the gaussian mixture probability density function given by the
%   centers mu and the corresponding covariance matrices SIGMA at
%   position x
%
%       mu      is a p-by-d matrix where each of the p rows represents the
%               (d-dimensional) position of a center
%
%       SIGMA   is a p-by-d-by-d array where SIGMA(j,:,:) is the d-by-d
%               covariance matrix corresponding to the j-th center
%
%       x       1-by-d position where the function is evaluated

p = size(mu,1);

result = 0;
for j=1:p
    result = result + gauss(x, mu(j,:), squeeze(SIGMA(j,:,:)));
end

end

