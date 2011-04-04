function [ points, normals ] = load_apts( filename )
%LOAD_APTS read data from a *.apts file
%   [ points, normals ] = LOAD_APTS( filename ) load position data and
%   corresponding normals from a *.apts file.
%
%      points     position data. n x d matrix, each row one point.
%
%      normals    surface normal data. n x d matrix, each row one normal
%                 vector.
%
%      filename   path to the file to be read.

file_id = fopen(filename);

textscan(file_id, '%[^\n]', 7); % skip 7 lines

data = textscan(file_id, '%f %f %f ; %f ; %f %f %f ; %f %f %f %f');

status = fclose(file_id);
assert(status == 0)

points = [data{1:3}];
% radii = data{4};
normals = [data{5:7}];
% colors = [data{8:11}];

end

