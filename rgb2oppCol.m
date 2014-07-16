function [ oppCol_coords, indx_sat] = rgb2oppCol( rgb_coords, mu_s, sigma_s, rotation_matrix, options,  varargin) 
%Convert from RGB space to SIC space
defaultopt = struct('Normalize','on'); % flag for normalizing the coordinate

if nargin < 5
    options = [];
    if nargin < 4
          error('Need at least 3 inputs: rgb coordinates, mu, sigma of normalization function, and rotation matrix')
    end
end

% get the extreme cutoff and filter optical density and plot flag
normalizeflag = optimget(options,'Normalize',defaultopt,'fast');

if isempty(rotation_matrix)
    rotation_matrix = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3); ...
        1/sqrt(6) 1/sqrt(6) -2/sqrt(6); ...
        -1/sqrt(2) 1/sqrt(2) 0];
end
rotated_coordinates = rotation_matrix*rgb_coords;

% Eliminate chemical saturation (black) and optical saturation (white)
first_component = uint8(rotated_coordinates(1,:));
indx_chemical_sat = first_component < 5;
indx_optical_sat = first_component > 250;
indx_sat = indx_chemical_sat + indx_optical_sat;

if strcmpi(normalizeflag,'on')
    F = @(x, mu, sigma) 1-exp(- max(x - mu,0)/2/sigma^2);
    Fmatrix = F(sqrt(sum(rotated_coordinates(2:3,:).^2,1)),mu_s,sigma_s); 
    oppCol_coords = repmat(Fmatrix,2,1).* rotated_coordinates(2:3,:)./...
        repmat(1 + sqrt(sum(rotated_coordinates(2:3,:).^2,1)),2,1);
else
    oppCol_coords = rotated_coordinates(2:3,:);
end

end

