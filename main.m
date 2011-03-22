res = 400

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

% mu_normals = ...
% grad_weighted_signed_distance_fu(x, x_normals, ...
%                                  repmat(reshape(eye(2) * start_sigma^2, ...
%                                                 [1 2 2]), ...
%                                         [size(x,1) 1 1]), ...
%                                  mu);
% % normalize normals
% mu_normals = mu_normals ./ repmat(sqrt(mu_normals(:,1).^2+mu_normals(:,2).^2), [1 2]);
% 
% plot_f(mu, mu_normals, squeeze(SIGMA), corners, 'res',res)
% hold
% quiver(mu(:,1),mu(:,2),mu_normals(:,1),mu_normals(:,2),'color','black')
% scatter(x(:,1),x(:,2),'.r')

% Query points MUST be measurement points, reference points MUST be kernel centers,
% NOT the other way around. Else only the nearest k measurement points
% would 'see' a given kernel.
Xq = x; % query points, use input points for now

old_path = path;
addpath([pwd '/../ann_mwrapper'])

nnidx = annquery(mu', Xq', 10);

path(old_path)

[n dim] = size(Xq);
p = size(mu,1);

A = sparse(n,p);

for i = 1:n
    neighborMuIdx = nnidx(:, i);
    % The i'th measurement point might have less than size(nnidx, 1)
    % neighbors. Eliminate empty idx entires
    neighborMuIdx = neighborMuIdx(neighborMuIdx > 0);
    numNeighbors = length(neighborMuIdx);
    
    % FIXME: weighted_signed_distance_fu isn't the right thing to call
    % here: We need separate values per kernel, not the sum.
    A(i, neighborMuIdx) = weighted_signed_distance_fu(mu(neighborMuIdx,:), x_normals, ...
                                                      reshape(SIGMA(:,neighborMuIdx,:,:), [numNeighbors dim dim]), ...
                                                      Xq(i,:));
end