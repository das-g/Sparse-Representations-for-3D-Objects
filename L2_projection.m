x = (0:0.01:5)';
mu =(0:0.1:5)';
SIGMA = 1;

p = size(mu,1);
n = size(x,1);

K = zeros(n,p);

for j=1:p
    K(:,j) = gauss(x, repmat(mu(j,:), [n 1]), squeeze(SIGMA));
end

f = cos(7*x) + sin(5*x);

alpha = K\f;

f_reconstructed = K * alpha;

subplot(2,1,1)
plot(x,f)
hold on
plot(x,f_reconstructed,'color','red')
legend('f','reconstructed f')
hold off

subplot(2,1,2)
plot(x,f_reconstructed - f)
title('error')