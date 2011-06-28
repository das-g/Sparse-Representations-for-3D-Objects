function [ ] = plot_gauss_mix( mu, SIGMA, corners, varargin )
% Plot gaussian mixture given by centers mu and covariance matrices SIGMA
%
%   PLOT_GAUSS_MIX(mu, SIGMA, corners)
%
%       mu      is a p-by-d matrix where each of the p rows represents the
%               (d-dimensional) position of a center
%
%       SIGMA   is a p-by-d-by-d array where SIGMA(j,:,:) is the d-by-d
%               covariance matrix corresponding to the j-th center. A valid
%               covariance matrix must be positive-semidefinite. There will
%               be no warning when the passed matrices aren't.
%
%       corners matrix indicating the boundary of the area to be plotted.
%               Structure: [ <left edge> , <lower edge>;
%                            <right edge>, <upper edge>  ]
%
%   PLOT_GAUSS_MIX takes the following optional arguments as
%   'name',value pairs, after the mandatory parameters:
%
%       res     height and width resolution of the plot; default: 100
%
%       x_res   width resolution of the plot, overrides res if given
%
%       y_res   height resolution of the plot, overrides res if given
%
%       kernel_indices
%               (For debugging purposes.) If given, leave away all basis
%               functions whose index is not listed in this parameter,
%               otherwise plot the complete reconstruction.
%
% See also: gauss_mix_eval

%% Parse input arguments
ip = inputParser;

ip.addRequired('mu');
ip.addRequired('SIGMA');
ip.addRequired('corners');

ip.addParamValue('res', 100);
ip.addParamValue('x_res', false);
ip.addParamValue('y_res', false);
ip.addParamValue('kernel_indices', false);

ip.parse(mu, SIGMA, corners, varargin{:});

%% Post-process input arguments
if any(strcmpi('x_res', ip.UsingDefaults))
    x_res = ip.Results.res;
else
    x_res = ip.Results.x_res;
end

if any(strcmpi('y_res', ip.UsingDefaults))
    y_res = ip.Results.res;
else
    y_res = ip.Results.y_res;
end

if ~any(strcmpi('kernel_indices', ip.UsingDefaults))
    mu = mu(ip.Results.kernel_indices, :);
    SIGMA = SIGMA(ip.Results.kernel_indices, :, :);
end

%% Coordinates for each pixel
[X Y] = meshgrid( linspace(corners(1,1), corners(2,1), x_res), ...
                  linspace(corners(1,2), corners(2,2), y_res) );

%% Evaluate function for coordinates
Z = gauss_mix_eval(mu, SIGMA, [X(:) Y(:)]);

%% Deserialize function values
Z = reshape(Z, y_res, x_res);

%% Display

figure
imagesc([X(1,1) X(1,end)], [Y(1,1) Y(end,1)], Z)
set(gca,'YDir','normal')

end
