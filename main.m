res = 800;
filenamesnippet = 'mimi_sub2';
serialNo = '_1';
data = load(['../test/' filenamesnippet '.npoff']);
x = data(:,1:2);
x_normals = data(:,3:4);

[n dim] = size(x);
corners = [min(x); max(x)];
center = mean(corners,1);
corners = (corners - repmat(center,[2 1])) * 1.2 + repmat(center,[2 1]);

start_sigma_factor = 1.105;
goal_sigma = min(corners(2,:) - corners(1,:)) / 20;

%% Find distances to nearest neighbor using ann library.
% Save current path, so we can restore it later.
old_path = path;

% Change path, so that ann wrapper files are found.
addpath([pwd '/../ann_mwrapper']) % Full path needed, so prepend working dir

% Query for nearest 2 x (incl. itself) for each x.
[nnidx, dists] = annquery(x', x', 2);

% Clean up. (Restore previous path.)
path(old_path)

% The first row of dists will be all-zero, as those are the distances to
% from all points to themselfs. The second row is are the distances to the
% nearest points. We use the maximum of the latter as an estimate for the
% sampling width and use it to set the initial kernel size.
assert( all( dists(1, :) == 0) )
start_sigma = max( dists(2, :) ) * start_sigma_factor;

%% plot input
figure
scatter(x(:,1), x(:,2))
axis off
axis equal
hold on
quiver(x(:,1),x(:,2),x_normals(:,1),x_normals(:,2),'color','black')
hold off
set(gcf, 'PaperPositionMode', 'auto', ...
         'units','normalized', ...
         'outerposition',[0 0 1 1]) % maximize before saving
saveas(gcf, [filenamesnippet serialNo '_input'], 'epsc')

corners = [xlim' ylim'];

%% plot surface definition (signed distance fu. based on input)
plot_f(x, x_normals, ...
       repmat(reshape(eye(dim) * start_sigma^2, ...
                      [1 dim dim]), ...
              [size(x,1) 1 1]), ...
       corners, ...
       'res', res)
colorbar
axis off
axis equal
hold on
scatter(x(:,1),x(:,2),'.k')
hold off
set(gcf, 'PaperPositionMode', 'auto', ...
         'units','normalized', ...
         'outerposition',[0 0 1 1]) % maximize before saving
saveas(gcf, [filenamesnippet serialNo '_signed_distance_fu'], 'epsc')
% quiver(x(:,1),x(:,2),x_normals(:,1),x_normals(:,2),'color','black')
% 
% figure

%% Run EM
[mu SIGMA] = EM(x, ...
                'a', 1/n, ...
                'initial_sigma', start_sigma, ...
                'target_sigma', goal_sigma, ...
                'centers_to_points_ratio', 1, ...
                'max_steps', 40);

%% plot signed distance fu with new kernels from EM
plot_f(mu, x_normals, squeeze(SIGMA), corners, 'res',res)
colorbar
axis off
axis equal
hold on
% quiver(mu(:,1),mu(:,2),mu_normals(:,1),mu_normals(:,2),'color','black')
scatter(x(:,1),x(:,2),'.k')
hold off
set(gcf, 'PaperPositionMode', 'auto', ...
         'units','normalized', ...
         'outerposition',[0 0 1 1]) % maximize before saving
saveas(gcf, [filenamesnippet serialNo '_signed_distance_fu_after_EM'], 'epsc')

%% Find approximately sparse coefficients with L1 minimization
% Query points MUST be measurement points, reference points MUST be kernel centers,
% NOT the other way around. Else only the nearest k measurement points
% would 'see' a given kernel.
Xq = [x - x_normals * start_sigma * 0.5;
      x;
      x + x_normals * start_sigma * 0.5]; % query points
%mu = x; % reference points, use input points for now
mu_normals = x_normals; % normals at reference points (i.e. at kernel centers)

p = size(mu,1);

A = measurement_matrix(mu, mu_normals, SIGMA, Xq);

% right-hand side (measured f)
rhs = weighted_signed_distance_fu( x, x_normals, ...
                                   repmat(reshape(eye(dim) * start_sigma^2, ...
                                                  [1 dim dim]), ...
                                   [p 1 1]), Xq );
old_path = path;
addpath([pwd '/../l1_ls_matlab'])

lambda = 0.0001;
[coeff_L1ls status] = l1_ls(A, rhs, lambda, 1e-3, true);

path(old_path)

assert(all(status == 'Solved'))

%% sparsify (eliminate kernels corresponding to low coefficients)
threshold = 1.4;
coeff_L1ls_nonzero_idx = abs(coeff_L1ls) > threshold; % vector of booleans

coeff_L1ls_reduced = coeff_L1ls(coeff_L1ls_nonzero_idx);
mu_reduced = mu(coeff_L1ls_nonzero_idx, :);
mu_normals_reduced = mu_normals(coeff_L1ls_nonzero_idx, :);
SIGMA_reduced = SIGMA(:, coeff_L1ls_nonzero_idx, :, :);


%% plot kept points
figure
scatter(mu_reduced(:,1),mu_reduced(:,2),'r')
axis off
axis equal
set(gcf, 'PaperPositionMode', 'auto', ...
         'units','normalized', ...
         'outerposition',[0 0 1 1]) % maximize before saving
saveas(gcf, [filenamesnippet serialNo '_kept_points'], 'epsc')

%% reconstruct the target function from the reduced data
plot_approx(mu_reduced, mu_normals_reduced, SIGMA_reduced, ...
            coeff_L1ls_reduced, corners, 'res', res)
colorbar
axis off
axis equal
hold on
scatter(x(:,1),x(:,2),'.k')
hold off
set(gcf, 'PaperPositionMode', 'auto', ...
         'units','normalized', ...
         'outerposition',[0 0 1 1]) % maximize before saving
saveas(gcf, [filenamesnippet serialNo 'reconstruction_from_reduced_data'], 'epsc')
