function [ ] = plot_gauss_mix( mu, SIGMA, x_res, y_res )
%PLOT_GAUSS_MIX Plot gaussian mixture
%   given by centers mu and covariance matrices SIGMA   

corners = [min(mu); max(mu)];

center = mean(corners,1);
corners = (corners - repmat(center,[2 1])) * 1.2 + repmat(center,[2 1]);

[X Y] = meshgrid( linspace(corners(1,1), corners(2,1), x_res), ...
                  linspace(corners(1,2), corners(2,2), y_res) );

Z = gauss_mix_eval(mu, squeeze(SIGMA), [X(:) Y(:)]);

Z = reshape(Z, x_res, y_res);

figure
imagesc([X(1,1) X(1,end)], [Y(1,1) Y(end,1)], Z)
set(gca,'YDir','normal')

end

