function [ output_args ] = EM( x )
%EM Summary of this function goes here
%   Detailed explanation goes here

max_steps = 10;
n = size(x, 1)
dim = size(x, 2)

%% Iteration Start Values
mu = x;
pi_j = ones(size(x, 1), 1);
sigma = 1;
SIGMA = repmat(eye(dim),[1 1 n]) * sigma^2;

%% EM Iteration
for step=1:max_steps
    %% E-step
    pi_ij=zeros(n);
    for j=1:n
        pi_ij(:,j) = pi_j(j) * gauss(x, repmat(mu(j,:),n,1), SIGMA(:,:,j));
    end
    pi_ij = pi_ij ./ repmat(sum(pi_ij, 2),1,n);
    
    %...
    %% M-step
end

end