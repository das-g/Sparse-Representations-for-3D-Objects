function [ ] = plot_gauss_mix( mu, SIGMA, varargin )
%PLOT_GAUSS_MIX Plot gaussian mixture
%   given by centers mu and covariance matrices SIGMA

%% Default values for optional arguments
x_res = 100;
y_res = 100;

%% Override values for optional arguments that have been given
optargin = size(varargin,2);
assert(mod(optargin,2) == 0, 'optional arguments have to be provided in pairs: name, value')
for optarg=1:2:optargin
    switch varargin{optarg}
        case 'res'
            x_res = varargin{optarg + 1};
            y_res = varargin{optarg + 1};
        case 'x_res'
            x_res = varargin{optarg + 1};
        case 'y_res'
            y_res = varargin{optarg + 1};
        otherwise
            error(['unknown optional argument name: ' varargin{optarg}] )
    end % switch optarg
end % for optarg

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

