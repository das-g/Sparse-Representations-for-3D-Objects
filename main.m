size_range = 50:17:1750;
time_build_measurement_matrix = zeros(size(size_range));
time_l1_ls = zeros(size(size_range));
num_coeffs = zeros(size(size_range));
num_coeffs_reduced = zeros(size(size_range));

for i = 1:length(size_range);

x = unit_circle(size_range(i)); % 100 points on a unit cicle centered on origin
x_normals = x; % Yes, for an origin centered unit circle, this works.

corners = [min(x); max(x)];
center = mean(corners,1);
corners = (corners - repmat(center,[2 1])) * 1.2 + repmat(center,[2 1]);

start_sigma = 0.2;
goal_sigma = 5;

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

smoothed_normals = ...
grad_weighted_signed_distance_fu(x, x_normals, ...
                                 repmat(reshape(eye(2) * start_sigma^2, ...
                                                [1 2 2]), ...
                                        [size(x,1) 1 1]), ...
                                 mu);
% normalize normals
smoothed_normals = smoothed_normals ./ repmat(sqrt(smoothed_normals(:,1).^2+smoothed_normals(:,2).^2), [1 2]);

% plot_f(mu, mu_normals, squeeze(SIGMA), corners, 'res',res)
% hold
% quiver(mu(:,1),mu(:,2),mu_normals(:,1),mu_normals(:,2),'color','black')
% scatter(x(:,1),x(:,2),'.r')

% Query points MUST be measurement points, reference points MUST be kernel centers,
% NOT the other waz around. Else only the nearest k measurement points
% would 'see' a given kernel.
Xq = [x - smoothed_normals * start_sigma * 0.5;
      x;
      x + smoothed_normals * start_sigma * 0.5]; % query points
%mu = x; % reference points, use input points for now
mu_normals = x_normals; % normals at reference points (i.e. at kernel centers)

[n dim] = size(Xq);
p = size(mu,1);

tic
A = measurement_matrix(mu, mu_normals, SIGMA, Xq);
time_build_measurement_matrix(i) = toc;

% right-hand side (measured f)
rhs = weighted_signed_distance_fu( mu, mu_normals, ...
                                   repmat(reshape(eye(2) * start_sigma^2, ...
                                                  [1 dim dim]), ...
                                   [p 1 1]), Xq );
old_path = path;
addpath([pwd '/../l1_ls_matlab'])

lambda = 0.01;
tic
[coeff_L1ls status] = l1_ls(A, rhs, lambda, 1e-3, true);
time_l1_ls(i) = toc;

path(old_path)

assert(all(status == 'Solved'))

% sparsify (eliminate almost-zero entries)
threshold = 0.1;
coeff_L1ls_nonzero_idx = abs(coeff_L1ls) > threshold; % vector of booleans

coeff_L1ls_reduced = coeff_L1ls(coeff_L1ls_nonzero_idx);
mu_reduced = mu(coeff_L1ls_nonzero_idx, :);
mu_normals_reduced = mu_normals(coeff_L1ls_nonzero_idx, :);
SIGMA_reduced = SIGMA(:, coeff_L1ls_nonzero_idx, :, :);

num_coeffs(i) = numel(coeff_L1ls);
num_coeffs_reduced(i) = numel(coeff_L1ls_reduced);

end

save('circle_measurements_2.mat', 'size_range', 'time_build_measurement_matrix', 'time_l1_ls', 'num_coeffs', 'num_coeffs_reduced')
