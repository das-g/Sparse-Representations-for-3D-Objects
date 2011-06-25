function [ ] = plot_f( x, normals, SIGMA, corners, varargin )
% Plot weighted signed distance function given by points x and covariance
% matrices SIGMA
%
%   PLOT_F(x, normals, SIGMA, corners)
%
%       x       is a n-by-d matrix where each of the n rows represents the
%               (d-dimensional) position of a center
%
%       normals is a n-by-d matrix where each of the n rows represents the
%               (d-dimensional) directed normal of the zero-isosurface
%               at a center
%
%       SIGMA   is a n-by-d-by-d array where SIGMA(i,:,:) is the d-by-d
%               covariance matrix corresponding to the i-th point. A valid
%               covariance matrix must be positive-semidifinite. There will
%               be no warning when the passed matrices aren't.
%
%       corners matrix indicating the boundary of the area to be plotted.
%               Structure: [ <left edge> , <lower edge>;
%                            <right edge>, <upper edge>  ]
%
%   PLOT_F takes the following optional arguments as 'name',value pairs,
%   after the mandatory parameters:
%
%       res     height and width resolution of the plot; default: 100
%
%       x_res   width resolution of the plot, overrides res if given
%
%       y_res   height resolution of the plot, overrides res if given
%
%       x_indices
%               (For debugging purposes.) If given, leave away all basis
%               functions whose index is not listed in this parameter,
%               otherwise plot the complete reconstruction.
%
% See also: weighted_signed_distance_fu, plot_gauss_mix, plot_approx

%% Parse input arguments
ip = inputParser;

ip.addRequired('x');
ip.addRequired('normals');
ip.addRequired('SIGMA');
ip.addRequired('corners');

ip.addParamValue('res', 100);
ip.addParamValue('x_res', false);
ip.addParamValue('y_res', false);
ip.addParamValue('x_indices', false);

ip.parse(x, normals, SIGMA, corners, varargin{:});

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

if ~any(strcmpi('x_indices', ip.UsingDefaults))
    x = x(ip.Results.x_indices, :);
    SIGMA = SIGMA(ip.Results.x_indices, :, :);
end

%% Coordinates for each pixel
[X Y] = meshgrid( linspace(corners(1,1), corners(2,1), x_res), ...
                  linspace(corners(1,2), corners(2,2), y_res) );

%% Evaluate function for coordinates
Z = weighted_signed_distance_fu(x, normals, SIGMA, [X(:) Y(:)]);

%% Deserialize function values
Z = reshape(Z, y_res, x_res);

%% Display
figure
imagesc([X(1,1) X(1,end)], [Y(1,1) Y(end,1)], Z)
set(gca,'YDir','normal')
c = caxis;
caxis(c - c([2 1]))
colormap(myjet)

end
