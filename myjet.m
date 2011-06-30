function J = myjet(m)
%MYJET    Variant of HSV
%   Colormap inspired by MATLAB's jet, but with a hard transition
%   from cyan to yellow in the middle.
%
%   See also JET, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
middle = round(m/2);
third = round(m/3);
sixth = middle - third;

J = zeros(m,3);
J(middle+1:end-sixth, 1) = 1;
J(sixth:middle,3) = 1;
J(1:sixth, 3) = (sixth+1:2*sixth)/(2*sixth);
J(end:-1:end-sixth,1) = (sixth:2*sixth)/(2*sixth);
J(sixth+1:middle, 2) = (1:third)/third;
J(middle:end-sixth-1, 2) = (third:-1:1)/third;
