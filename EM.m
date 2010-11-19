function [ mu SIGMA ] = EM( x )
%EM Summary of this function goes here
%   Detailed explanation goes here

max_steps = 10;
n = size(x, 1)
dim = size(x, 2)

S_x = cov(x);
a = 0.01;

%% Iteration Start Values
mu = x;
p = size(mu,1);
pi_j = ones(1, n); % row vector, because we treat j as the 2nd dimension
sigma = 1;
SIGMA = repmat(reshape(eye(dim),[1 1 dim dim]),[1 p 1 1]) * sigma^2;

%% EM Iteration
for step=1:max_steps
    %% E-step
    pi_ij=zeros(n);
    for j=1:n
        pi_ij(:,j) = pi_j(j) * gauss(x, repmat(mu(j,:),n,1), squeeze(SIGMA(:,j,:,:)));
    end
    pi_ij = pi_ij ./ repmat(sum(pi_ij, 2),1,n);
    
    %% M-step
    pi_j = mean(pi_ij, 1); % Mean along i; results in a row vector
    mu = (pi_ij' * x) ./ repmat(n * pi_j', [1 dim]);
    tmp = repmat(reshape(x,[n 1 dim]),[1 p 1]) ...
        - repmat(shiftdim(mu,-1),[n 1 1]);
    S_j = dot(repmat(repmat(pi_ij, [1 1 2]) .* tmp,[1 1 1 2]),repmat(reshape(tmp,[n p 1 dim]),[1 1 2 1]), 1);
    SIGMA = (2 * a * repmat(shiftdim(S_x,-2),[1 p 1 1]) + S_j) ...
            ./ (2 * a + n * repmat(pi_j,[1 1 dim dim]));
end

end