%% inputs
f = @(x) cos(7 * x) + sin(5 * x);

x_min = 0;
x_max = 5;
x_measurement = (x_min:0.01:x_max)';
x_reconstruct = (x_min:0.01:x_max)';

mu = (x_min:0.1:x_max)';
SIGMA = 1;

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

alpha = K_measurement \ f_measured;

%% reconstruct f
f_reconstructed = K_reconstruct * alpha;
%% measure f at reconstruction points (for plotting & comparison)
f_original = f(x_reconstruct); % actual value of f at reconstruction points

%% Plot original and reconstructed f in upper half
subplot(2, 1, 1)
plot(x_reconstruct, f_original)
hold on
plot(x_reconstruct, f_reconstructed, 'color', 'red')
legend('f','reconstructed f')
hold off

%% Plot error in lower half
subplot(2, 1, 2)
plot(x_reconstruct, f_reconstructed - f_original)
title('error')