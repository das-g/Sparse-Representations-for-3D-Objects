res = 400;

[ x, x_normals ] = load_apts('../Data/3D/Rabbit.apts');

[n dim] = size(x);
corners = [min(x); max(x)];
center = mean(corners,1);
corners = (corners - repmat(center,[2 1])) * 1.2 + repmat(center,[2 1]);

start_sigma = 7;
goal_sigma = 130;

% plot_f(x, x_normals, ...
%        repmat(reshape(eye(dim) * start_sigma^2, ...
%                       [1 dim dim]), ...
%               [size(x,1) 1 1]), ...
%        corners, ...
%        'res', res)
% hold
% quiver(x(:,1),x(:,2),x_normals(:,1),x_normals(:,2),'color','black')
% 
% figure

disp('EM start')
[mu SIGMA] = EM(x, ...
                'a', 0.01, ...
                'initial_sigma', start_sigma, ...
                'target_sigma', goal_sigma, ...
                'centers_to_points_ratio', 1, ...
                'max_steps', 40);
disp('EM done')

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

p = size(mu,1);

disp('Measurement Matrix start')
A = measurement_matrix(mu, mu_normals, SIGMA, Xq, ...
                       'compile', true);
disp('Measurement Matrix done')

disp('RHS start')
% right-hand side (measured f)
rhs = weighted_signed_distance_fu( x, x_normals, ...
                                   repmat(reshape(eye(dim) * start_sigma^2, ...
                                                  [1 dim dim]), ...
                                   [p 1 1]), Xq );
disp('RHS done\n')

old_path = path;
addpath([pwd '/../l1_ls_matlab'])

lambda = 0.01;
disp('L1 start')
[coeff_L1ls status] = l1_ls(A, rhs, lambda, 1e-3, true);
disp('L1 done')

path(old_path)

assert(all(status == 'Solved'))

% sparsify (eliminate almost-zero entries)
threshold = 0.1;
coeff_L1ls_nonzero_idx = abs(coeff_L1ls) > threshold; % vector of booleans

coeff_L1ls_reduced = coeff_L1ls(coeff_L1ls_nonzero_idx);
mu_reduced = mu(coeff_L1ls_nonzero_idx, :);
mu_normals_reduced = mu_normals(coeff_L1ls_nonzero_idx, :);
SIGMA_reduced = SIGMA(:, coeff_L1ls_nonzero_idx, :, :);

% reconstruct the target function from the reduced data
%plot_approx(mu_reduced, mu_normals_reduced, SIGMA_reduced, ...
%            coeff_L1ls_reduced, corners, 'res', 200)
