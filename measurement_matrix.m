function A = measurement_matrix(mu, normals, SIGMA, x, varargin)

%% Parse input arguments
ip = inputParser;

ip.addRequired('mu');
ip.addRequired('normals');
ip.addRequired('SIGMA');
ip.addRequired('x');

ip.addParamValue('compile', false);

ip.parse(mu, normals, SIGMA, x, varargin{:});

%% Find nearest neighbors using ann library.
% Save current path, so we can restore it later.
old_path = path;

% Change path, so that ann wrapper files are found.
addpath([pwd '/../ann_mwrapper']) % Full path needed, so prepend working dir

% Query for nearest 10 mu for each x.
nnidx = annquery(mu', x', 10);

% Clean up. (Restore previous path.)
path(old_path)

%% Determine sizes of input arguments
n = size(x, 1);
p = size(mu, 1);
p_nn = size(nnidx, 1);

%% Indices of measurement matrix' non-zero values
A_i = repmat((1:n)', [1 p_nn]);
A_j = double(nnidx');

%% (Non-zero) values of measurement matrix
% We call measurement_matrix_values to get the values. For large n and/or
% large p_nn, it pays off to compile it, first. The caller can request
% compilation by passing ..., 'compile', true.

if ip.Results.compile
    %% Caller requested compile
    % The sizes of passed parameters are not known in advance, so we compile
    % measurement_matrix_values just before calling it. For large n or p_nn
    % this still pays off.
    emlmex -o mmv_compiled measurement_matrix_values -eg {mu, normals, SIGMA, x, nnidx}
    
    % Call compiled function.
    A_values = mmv_compiled(mu, normals, SIGMA, x, nnidx);
else
    %% Caller didn't request compile
    % Make sure we'll call the *.m file, not a compiled function.
    assert(exist('measurement_matrix_values', 'file') == 2, ... Check whether *.m function.
           'Please delete the measurement_matrix_values MEX file.')
    
    % Call it.
    A_values = measurement_matrix_values(mu, normals, SIGMA, x, nnidx);
end

%% Compose measurement matrix
A = sparse(A_i(nnidx' > 0), ...
           A_j(nnidx' > 0), ...
           A_values(nnidx' > 0), ...
           n, p);
