function [ phi ] = gauss( x, mu, Sigma) %#eml
% Multivariate Gauss distribution
%
% phi = GAUSS(x, mu, Sigma)
%   Evaluate the multivariate gauss function centered at mu with covariance
%   matrix Sigma at position x.
%
%       mu      1-by-d vector representing the (d-dimensional) position of
%               the center
%
%       Sigma   the d-by-d matrix covariance matrix. A valid covariance
%               matrix must be positive-semidefinite. There will be no
%               warning when the passed matrix isn't.
%
%       x       1-by-d vector representing the position where the function
%               shall be evaluated
%
%   For performance, you may evaluate several Gauss distributions with the
%   same covariance matrix with a single call:
%
%   phis = GAUSS(xs, mus, Sigma)
%
%       mus     is a n-by-d array where each of the n rows represents the
%               (d-dimensional) position of the center of one of n
%               multivariate Gauss distributions
%
%       Sigma   is a d-by-d matrix covariance matrix, shared by all n
%               evaluated functions. A valid covariance matrix must be
%               positive-semidefinite. There will be no warning when the
%               passed matrix isn't.
%
%       xs      n-by-d array where the ith row represents the position
%               where the ith function shall be evaluated
%
%   The n results will be returned as the rows of the n-by-1 array phis.
%
%   Example: To evaluate the univariate gauss function centered at c with
%            unit variance at points p1, p2 and p3 (c, p1, p2 and p3 all 2D
%            row vectors) run
%
%               gauss( [p1; p2; p3], repmat(c, [3 1]), eye(2) )
%
%   See also: GAUSS_MIX_EVAL
dim = size(x,2);
phi =  (2 * pi) ^ (-dim/2) * 1 / sqrt(det(Sigma)) ...
       * exp(-0.5 * dot((x - mu) / Sigma, (x - mu), 2) );
end
