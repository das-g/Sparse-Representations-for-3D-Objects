x_measurement = (0:0.01:5)';
x_reconstruct = (0:0.01:5)';
mu =(0:0.1:5)';
SIGMA = 1;

p = size(mu,1);
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

f = @(x) cos(7*x) + sin(5*x);

f_measured = f(x_measurement);

alpha = K_measurement \ f_measured;

f_reconstructed = K_reconstruct * alpha;
f_original = f(x_reconstruct); % actual value of f at reconstruction points

subplot(2,1,1)
plot(x_reconstruct, f_original)
hold on
plot(x_reconstruct, f_reconstructed, 'color', 'red')
legend('f','reconstructed f')
hold off

subplot(2,1,2)
plot(x_reconstruct, f_reconstructed - f_original)
title('error')