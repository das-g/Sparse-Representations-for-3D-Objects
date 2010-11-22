function [ mu SIGMA ] = EM( x )
%EM Summary of this function goes here
%   Detailed explanation goes here

max_steps = 10;
[n dim] = size(x)

S_x = cov(x);
a = 0.01;

%% Iteration Start Values
random_order = randperm(n);
p = round(n / 3)
mu = x(random_order(1:p),:);
pi_j = ones(1, p); % row vector, because we treat j as the 2nd dimension
sigma = 1;
SIGMA = repmat(reshape(eye(dim),[1 1 dim dim]),[1 p 1 1]) * sigma^2;

% Plot initial values
hold off
scatter(x(:,1),x(:,2))
hold on

%% EM Iteration
for step=1:max_steps
    %% E-step
    pi_ij=zeros(n,p);
    for j=1:p
        pi_ij(:,j) = pi_j(j) * gauss(x, repmat(mu(j,:),n,1), squeeze(SIGMA(:,j,:,:)));
    end
    pi_ij = pi_ij ./ repmat(sum(pi_ij, 2),1,p);
    
    %% M-step
    pi_j = mean(pi_ij, 1); % Mean along i; results in a row vector
    mu = (pi_ij' * x) ./ repmat(n * pi_j', [1 dim]);
    tmp = repmat(reshape(x,[n 1 dim]),[1 p 1]) ...
        - repmat(shiftdim(mu,-1),[n 1 1]);
    S_j = dot(repmat(repmat(pi_ij, [1 1 2]) .* tmp,[1 1 1 2]),repmat(reshape(tmp,[n p 1 dim]),[1 1 2 1]), 1);
    SIGMA = (2 * a * repmat(shiftdim(S_x,-2),[1 p 1 1]) + S_j) ...
            ./ (2 * a + n * repmat(pi_j,[1 1 dim dim]));

    %% Plot intermediate results
    scatter(mu(:,1),mu(:,2))
end

end