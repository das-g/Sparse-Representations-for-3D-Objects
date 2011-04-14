res = 400;

data = load('../Data/2D/arrow01.npoff');
x = data(:,1:2);
x_normals = data(:,3:4);

corners = [min(x); max(x)];
center = mean(corners,1);
corners = (corners - repmat(center,[2 1])) * 1.2 + repmat(center,[2 1]);

sigma_range = 1:0.3:15;

errors = zeros(size(sigma_range));
num_coeffs_reduced = zeros(size(sigma_range));

for i = 1:length(sigma_range)

    start_sigma = sigma_range(i);

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

    mu = x;
    SIGMA = repmat(reshape(eye(2) * start_sigma^2, ...
                           [1 1 2 2]), ...
                      [1 size(x,1) 1 1]);

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
    rhs = weighted_signed_distance_fu( x, x_normals, ...
                                       repmat(reshape(eye(2) * start_sigma^2, ...
                                                      [1 dim dim]), ...
                                       [p 1 1]), Xq );
    old_path = path;
    addpath([pwd '/../l1_ls_matlab'])

    lambda = 0.01;
    [coeff_L1ls status] = l1_ls(A, rhs, lambda, 1e-3, true);

    path(old_path)

    assert(all(status == 'Solved'))

    % sparsify (eliminate almost-zero entries)
    threshold = 0.1;
    coeff_L1ls_nonzero_idx = abs(coeff_L1ls) > threshold; % vector of booleans

    coeff_L1ls_reduced = coeff_L1ls(coeff_L1ls_nonzero_idx);
    mu_reduced = mu(coeff_L1ls_nonzero_idx, :);
    mu_normals_reduced = mu_normals(coeff_L1ls_nonzero_idx, :);
    SIGMA_reduced = SIGMA(:, coeff_L1ls_nonzero_idx, :, :);

    % Measure the error of the reconstructed distance fu. at surface points
    % (from input) and off-surface points (generated from normals).
    A_reduced = measurement_matrix(mu_reduced, ...
                                   mu_normals_reduced, ...
                                   SIGMA_reduced, ...
                                   Xq);
    error = A_reduced * coeff_L1ls_reduced - rhs;
    errors(i) = sum( error .^ 2 );

    % Record number of kept coefficients
    num_coeffs_reduced(i) = numel(coeff_L1ls_reduced);

end

save('sigma_measurements_1.mat', 'sigma_range', ...
                                 'errors', ...
                                 'num_coeffs_reduced')
