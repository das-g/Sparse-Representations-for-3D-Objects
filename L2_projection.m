%% inputs
f = @(x) gauss(x, 0, 25);

x_min = -25;
x_max = 25;
x_measurement = (x_min:1:x_max)';
x_reconstruct = (x_min:0.01:x_max)';

mu = (x_min:1:x_max)';
SIGMA = 4;

%% prepare search path
% save path before we manipulate it (for restoring it later)
old_path = path;

% addpath want absolute paths, so prepend pwd/
addpath([pwd '/../l1magic-1.1/Optimization'])

%% build measurement and reconstruction matrix
p = size(mu, 1);
n_measurement = size(x_measurement, 1);
n_reconstruct = size(x_reconstruct, 1);

K_measurement = zeros(n_measurement, p);
K_reconstruct = zeros(n_reconstruct, p);

for j=1:p
    K_measurement(:,j) = gauss(x_measurement, ...
                               repmat(mu(j,:), [n_measurement 1]), ...
                               squeeze(SIGMA));
    K_reconstruct(:,j) = gauss(x_reconstruct, ...
                               repmat(mu(j,:), [n_reconstruct 1]), ...
                               squeeze(SIGMA));
end

%% measure f
f_measured = f(x_measurement);

% L2 approximation also used as start point for L1 iterations
alpha_L2 = K_measurement \ f_measured;

alpha_L1eq = l1eq_pd(alpha_L2, K_measurement, [], f_measured);

%% reconstruct f
f_reconstructed_L2 = K_reconstruct * alpha_L2;
f_reconstructed_L1eq = K_reconstruct * alpha_L1eq;

%% measure f at reconstruction points (for plotting & comparison)
f_original = f(x_reconstruct); % actual value of f at reconstruction points

%% Plot original and reconstructed f in upper half
subplot(2, 1, 1)
plot(x_reconstruct, f_original)
hold on
plot(x_reconstruct, f_reconstructed_L2, 'color', 'red')
plot(x_reconstruct, f_reconstructed_L1eq, 'color', 'green')
legend('f', ...
       'reconstructed f (from L2)', ...
       'reconstructed f (from L1 with eq. cond.)')
hold off

%% Plot error in lower half
subplot(2, 1, 2)
plot(x_reconstruct, f_reconstructed_L2 - f_original, 'color', 'red')
hold on
plot(x_reconstruct, f_reconstructed_L1eq - f_original, 'color', 'green')
hold off
title('error')

%% restore old path
path(old_path)
