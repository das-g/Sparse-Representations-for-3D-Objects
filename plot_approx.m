function [ ] = plot_approx( x, normals, SIGMA, coeffs, corners, varargin )
%PLOT_APPROX plot approximated function given by coefficients coeffs
%   given by points x and covariance matrices SIGMA
%
%       x       is a n-by-d matrix where each of the n rows represents the
%               (d-dimensional) position of a center
%
%       SIGMA   is a n-by-d-by-d array where SIGMA(i,:,:) is the d-by-d
%               covariance matrix corresponding to the i-th point
%
%       coeffs  is a n-by-1 array contaning coefficients for weighting
%               each basis function
%
%       corners matrix indicating the boundary of the area to be plotted.
%               Structure: [ <left edge> , <lower edge>;
%                            <right edge>, <upper edge>  ]

%% Parse input arguments
ip = inputParser;

ip.addRequired('x');
ip.addRequired('normals');
ip.addRequired('SIGMA');
ip.addRequired('coeffs');
ip.addRequired('corners');

ip.addParamValue('res', 100);
ip.addParamValue('x_res', false);
ip.addParamValue('y_res', false);
ip.addParamValue('x_indices', false);

ip.parse(x, normals, SIGMA, coeffs, corners, varargin{:});

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

%%

[X Y] = meshgrid( linspace(corners(1,1), corners(2,1), x_res), ...
                  linspace(corners(1,2), corners(2,2), y_res) );

A = measurement_matrix(x, normals, SIGMA, [X(:) Y(:)]);
Z = A * coeffs;

Z = reshape(Z, y_res, x_res);

figure
imagesc([X(1,1) X(1,end)], [Y(1,1) Y(end,1)], Z)
set(gca,'YDir','normal')

end

