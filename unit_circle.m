function [ x ] = unit_circle( n )
%UNIT_CIRCLE uniform points and corresponding normals on unit circle
%   x = UNIT_CIRCLE(n) computes n uniformly distributed points x on a unit
%   circle.

angles = ( 0 : 2 * pi / n : 2 * pi * (1 - 1 / n) )';

x = [cos(angles) sin(angles)];

end
