function [ output_args ] = EM( x )
%EM Summary of this function goes here
%   Detailed explanation goes here

max_steps = 10;
[n dim] = size(x)

%% Iteration Start Values
mu = x;
pi_j = ones(size(x, 1), 1);

%% EM Iteration
for step=1:max_steps
    %% E-step
    pi_ij=zeros(n);
    for j=1:n
        pi_ij(:,j) = pi_j(j) * gauss(x, repmat(mu(j,:),n,1), eye(dim)); % TODO use real sigma instead of eye
    end
    pi_ij = pi_ij ./ repmat(sum(pi_ij, 2),1,n);
    
    %...
    %% M-step
end

end