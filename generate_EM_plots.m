res = 800;

[x x_normals] = unit_square(50,1);

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

%% plot input and selected kernels
figure
scatter(x(:,1), x(:,2))
axis off
axis equal
hold on
quiver(x(:,1),x(:,2),x_normals(:,1),x_normals(:,2),'color','black')

scatter(x([25 48],1), x([25 48],2),'r','filled')

hold off
set(gcf, 'PaperPositionMode', 'auto', ...
         'units','normalized', ...
         'outerposition',[0 0 1 1]) % maximize before saving
saveas(gcf, ['selected_kernels_02'], 'epsc')
close

%% plot initial kernel shapes
EM_step_plot( x, ...
       repmat(reshape(eye(dim) * start_sigma^2, ...
                      [1 dim dim]), ...
              [size(x,1) 1 1]), ...
        0, corners )

%% Run EM
[mu SIGMA] = EM(x, ...
                'a', 1/n, ...
                'initial_sigma', start_sigma, ...
                'target_sigma', goal_sigma, ...
                'centers_to_points_ratio', 1, ...
                'step_plot', @(mu_, SIGMA_, step_) EM_step_plot(mu_, SIGMA_, step_, corners), ...
                'max_steps', 41);