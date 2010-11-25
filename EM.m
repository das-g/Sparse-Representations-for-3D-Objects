function [ mu SIGMA ] = EM( x, varargin)
%EM Summary of this function goes here
%   Detailed explanation goes here

%% Default values for optional arguments
a = 0.01;
max_steps = 10;
centers_to_points_ratio = 1;
pre_plot = @(mu,SIGMA) []; % as close to a NOOP as we get with anonymous functions
step_plot = @(mu,SIGMA) [];
post_plot = @(mu,SIGMA) [];


%% Override values for optional arguments that have been given
optargin = size(varargin,2);
assert(mod(optargin,2) == 0, 'optional arguments have to be provided in pairs: name, value')
for optarg=1:2:optargin
    switch varargin{optarg}
        case 'a'
            a = varargin{optarg + 1};
        case 'max_steps'
            max_steps = varargin{optarg + 1};
        case 'centers_to_points_ratio'
            centers_to_points_ratio = varargin{optarg + 1};
        case 'pre_plot'
            pre_plot = varargin{optarg + 1};
        case 'step_plot'
            step_plot = varargin{optarg + 1};
        case 'post_plot'
            post_plot = varargin{optarg + 1};
        otherwise
            error(['unknown optional argument name: ' varargin{optarg}] )
    end % switch optarg
end % for optarg

%% Derived values that aren't changed during the iteration
[n dim] = size(x)

S_x = cov(x);

%% Iteration Start Values
random_order = randperm(n);
p = round(n * centers_to_points_ratio)
mu = x(random_order(1:p),:);
pi_j = ones(1, p); % row vector, because we treat j as the 2nd dimension
sigma = 1;
SIGMA = repmat(reshape(eye(dim),[1 1 dim dim]),[1 p 1 1]) * sigma^2;

% Plot initial values
pre_plot(mu,SIGMA);

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
    step_plot(mu,SIGMA);
end

%% Plot final results
post_plot(mu,SIGMA);

end
