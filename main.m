res = 400;

data = load('../test/cartman.npoff');
x = data(:,1:2);
x_normals = data(:,3:4);

corners = [min(x); max(x)];
center = mean(corners,1);
corners = (corners - repmat(center,[2 1])) * 1.2 + repmat(center,[2 1]);

start_sigma = 7;
goal_sigma = 130;

% plot_f(x, x_normals, ...
%        repmat(reshape(eye(2) * start_sigma^2, ...
%                       [1 2 2]), ...
%               [size(x,1) 1 1]), ...
%        corners, ...
%        'res', res)
% hold
% quiver(x(:,1),x(:,2),x_normals(:,1),x_normals(:,2),'color','black')
% 
% figure

[mu SIGMA] = EM(x, ...
                'a', 0.01, ...
                'initial_sigma', start_sigma, ...
                'target_sigma', goal_sigma, ...
                'centers_to_points_ratio', 1, ...
                'max_steps', 40);

% plot_f(mu, mu_normals, squeeze(SIGMA), corners, 'res',res)
% hold
% quiver(mu(:,1),mu(:,2),mu_normals(:,1),mu_normals(:,2),'color','black')
% scatter(x(:,1),x(:,2),'.r')

% Query points MUST be measurement points, reference points MUST be kernel centers,
% NOT the other way around. Else only the nearest k measurement points
% would 'see' a given kernel.
Xq = [x - x_normals * start_sigma * 0.5;
      x;
      x + x_normals * start_sigma * 0.5]; % query points
%mu = x; % reference points, use input points for now
mu_normals = x_normals; % normals at reference points (i.e. at kernel centers)

[n dim] = size(Xq);
p = size(mu,1);

A = measurement_matrix(mu, mu_normals, SIGMA, Xq);

% right-hand side (measured f)
rhs = weighted_signed_distance_fu( mu, mu_normals, ...
                                   repmat(reshape(eye(2) * start_sigma^2, ...
                                                  [1 dim dim]), ...
                                   [p 1 1]), Xq );
old_path = path;
addpath([pwd '/../spgl1-1.7'])

tau = n/20;
coeff_spgl1 = spg_lasso(A, rhs, tau, struct('verbosity', 0));

path(old_path)

% sparsify (eliminate almost-zero entries)
threshold = 0.1;
coeff_spgl1_nonzero_idx = abs(coeff_spgl1) > threshold; % vector of booleans

coeff_spgl1_reduced = coeff_spgl1(coeff_spgl1_nonzero_idx);
mu_reduced = mu(coeff_spgl1_nonzero_idx, :);
mu_normals_reduced = mu_normals(coeff_spgl1_nonzero_idx, :);
SIGMA_reduced = SIGMA(:, coeff_spgl1_nonzero_idx, :, :);

% reconstruct the target function from the reduced data
plot_approx(mu_reduced, mu_normals_reduced, SIGMA_reduced, ...
            coeff_spgl1_reduced, corners, 'res', 200)
