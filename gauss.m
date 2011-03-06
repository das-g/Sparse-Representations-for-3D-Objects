function [ phi ] = gauss( x, mu, Sigma) %#eml
% Multivariate Gauss distribution
dim = size(x,2);
phi =  (2 * pi) ^ (-dim/2) * 1 / sqrt(norm(Sigma)) ...
       * exp(-0.5 * dot((x - mu) / Sigma, (x - mu), 2) );
end