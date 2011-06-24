function [ mu SIGMA ] = EM( x, varargin)
% Compute kernels of a multivariate guassian mixture approximating a point
% cloud, using a modified estimation maximization algorithm.
%
% [mu, SIGMA] = EM(x)
%
% Input:
%       x   n-by-d array, positions of points in point cloud
%
% Outputs:
%       mu  n-by-d array, positions of kernels
%
%       SIGMA
%           n-by-d-by-d array; SIGMA(j,:,:) is the covariance matrix
%           corresponding to the jth kernel. A valid covariance matrix must
%           be positive-semidefinite. There will be no warning when the
%           passed matrices aren't.
%
% EM also takes a set of optional arguments as 'name',value pairs.
%
% See also: plot_gauss_mix, gauss_mix_eval

%% Parse input arguments
ip = inputParser;

ip.addRequired('x');

ip.addParamValue('a', 0.01);
ip.addParamValue('max_steps', 10);
ip.addParamValue('centers_to_points_ratio', 1);
ip.addParamValue('initial_sigma', 1);
ip.addParamValue('target_sigma', 1);
ip.addParamValue('input_plot', @(x,SIGMA) []);
ip.addParamValue('pre_plot', @(x,SIGMA) []);
ip.addParamValue('step_plot', @(x,SIGMA) []);
ip.addParamValue('post_plot', @(x,SIGMA) []);

ip.parse(x, varargin{:});

%% Derived values that aren't changed during the iteration
[n dim] = size(x);

S_x = eye(dim) * ip.Results.target_sigma^2;

%% Iteration Start Values
p = round(n * ip.Results.centers_to_points_ratio);
if p < n
    % select random subset of input points as initial centers
    random_order = randperm(n);
    mu = x(random_order(1:p),:);
else % (p == n)
    % use all input points as initial centers
    mu = x;
end
pi_j = ones(1, p); % row vector, because we treat j as the 2nd dimension
SIGMA = repmat(reshape(eye(dim),[1 1 dim dim]),[1 p 1 1]) * ip.Results.initial_sigma^2;

% Plot input & initial values
ip.Results.input_plot(x,squeeze(SIGMA));
ip.Results.pre_plot(mu,squeeze(SIGMA));

%% EM Iteration
for step=1:ip.Results.max_steps
    %% E-step
    pi_ij=zeros(n,p);
    for j=1:p
        pi_ij(:,j) = pi_j(j) * gauss(x, repmat(mu(j,:),n,1), squeeze(SIGMA(:,j,:,:)));
    end
    pi_ij = pi_ij ./ repmat(sum(pi_ij, 2),1,p);
    
    %% M-step
    pi_j = mean(pi_ij, 1); % Mean along i; results in a row vector
    %mu = (pi_ij' * x) ./ repmat(n * pi_j', [1 dim]);
    tmp = repmat(reshape(x,[n 1 dim]),[1 p 1]) ...
        - repmat(shiftdim(mu,-1),[n 1 1]);
        S_j = zeros(1, p, dim, dim);
    for j=1:p
        S_j(:,j,:,:) = dot(repmat(repmat(pi_ij(:,j), [1 1 dim]) .* tmp(:,j,:),[1 1 1 dim]), ...
                           repmat(reshape(tmp(:,j,:),[n 1 1 dim]),[1 1 dim 1]), 1);
    end
    SIGMA = (2 * ip.Results.a * repmat(shiftdim(S_x,-2),[1 p 1 1]) + S_j) ...
            ./ (2 * ip.Results.a + n * repmat(pi_j,[1 1 dim dim]));

    %% Plot intermediate results
    ip.Results.step_plot(mu,squeeze(SIGMA));
end

%% Plot final results
ip.Results.post_plot(mu,squeeze(SIGMA));

end
