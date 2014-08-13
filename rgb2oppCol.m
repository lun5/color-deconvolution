function [ oppCol_coords, brightness, theta, rotated_coordinates] = rgb2oppCol( rgb_coords, rotation_matrix,  varargin) 
%Convert from RGB space to SIC space
if nargin < 2
    error('Need at least 2 inputs: rgb coordinates and rotation matrix')
end

if isempty(rotation_matrix)
    rotation_matrix = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3); ...
        1/sqrt(6) 1/sqrt(6) -2/sqrt(6); ...
        -1/sqrt(2) 1/sqrt(2) 0];
end
%od_coords = rgb2od(rgb_coords);
% normalize the rgb coordinates to [0 1]
%normalized_rgb = rgb_coords./255;
%rotated_coordinates = rotation_matrix*normalized_rgb;%rgb_coords;
rotated_coordinates = rotation_matrix*rgb_coords;

%rotated_coordinates = rotation_matrix*od_coords;
% Eliminate chemical saturation (black) and optical saturation (white)
brightness = rotated_coordinates(1,:);

% calculate angle to map to the unit circle
theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
x = cos(theta); y = sin(theta);
oppCol_coords = [x;y];

end

