function [ positions, normals ] = unit_square( n, rotation )
%UNIT_SQUARE uniform points and corresponding normals on unit square
%   x = UNIT_SQUARE(n) computes 4*n uniformly distributed points x on
%    a unit square.
%   
%   [x,normals] = UNIT_SQUARE(n) also computes the corresponding normals.
%   
%   x = UNIT_SQUARE(n, rotation)
%   or
%   [x,normals] = UNIT_SQUARE(n, rotation) generate a pre-rotated square
%   instead of an axis aligned one. rotation is in radian.

% n points on the right edge of a square (on the complex plane),
% inclusive the lower corner, exclusive the upper corner (which will be
% part of the upper edge)
positions = 1/2 + (-1/2:1/n:1/2-1/n) * 1i;

% Rotations by multiples of 90° to get the other edges.
rotations = [1 1i -1 -1i];

% If the user requested a rotation of the whole square, incorporate it into
% the rotations vector.
if nargin > 1
    rotations = rotations * exp(1i * rotation);
end

% Apply rotations. The result will be a 4n x 4 matrix. We'll serialize it
% in the next step.
positions = positions' * rotations;

% Go from complex pane to 2D real plane. (And serealize data.)
positions = [real(positions(:)) imag(positions(:))];

if nargout > 1
    % Normals for the right edge. The corner point gets a 45° normal. All
    % others are just perpendicular to the edge.
    normals = [exp(-pi*1i/4) ones(1, n - 1)];
    
    % Rotate and go C -> R^2. Same as for positions above.
    normals = normals' * rotations;
	normals = [real(normals(:)) imag(normals(:))];
end

end
