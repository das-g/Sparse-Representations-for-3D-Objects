function [ phi ] = gauss( x, mu, Sigma)
% Multivariate Gauss distribution
dim = size(x,2);
phi =  (2 * pi) ^ (-dim/2) * 1 / sqrt(det(Sigma)) ...
       * exp(-0.5 * dot((x - mu) / Sigma, (x - mu), 2) );
end