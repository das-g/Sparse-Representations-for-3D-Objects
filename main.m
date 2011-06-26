function [] = main(filename_snippet, serial_no, doEM)

if doEM
    disp(['=== START ' filename_snippet ' (' serial_no ') with EM  ==='])
else
    disp(['=== START ' filename_snippet ' (' serial_no ') without EM  ==='])
end

res = 400;

[ x, x_normals ] = load_apts(['../Data/3D/' filename_snippet '.apts']);

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
disp(['start_sigma=' num2str(start_sigma)])

%%

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

if doEM
    disp('EM start')
    output_filename_addition = ''; % nothing to add, this is our default
    [mu SIGMA] = EM(x, ...
                    'a', 0.01, ...
                    'initial_sigma', start_sigma, ...
                    'target_sigma', goal_sigma, ...
                    'centers_to_points_ratio', 1, ...
                    'max_steps', 40);
    disp('EM done')
else
    disp('(skipping EM, using spherical kernels)')
    output_filename_addition = '_noEM'; % so that we know EM was skipped
    mu = x;
    SIGMA = repmat(reshape(eye(dim) * start_sigma^2, ...
                           [1 1 dim dim]), ...
                   [1 size(x,1) 1 1]);
end

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

lambda = 0.00001;
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
txt_output([filename_snippet '_output' serial_no output_filename_addition '.txt'], coeff_L1ls, SIGMA, mu, mu_normals)
txt_output([filename_snippet '_output' serial_no output_filename_addition '_reduced.txt'], coeff_L1ls_reduced, SIGMA_reduced, mu_reduced, mu_normals_reduced)

if doEM
    disp(['=== END ' filename_snippet ' (' serial_no ') with EM ==='])
else
    disp(['=== END ' filename_snippet ' (' serial_no ') without EM ==='])
end
