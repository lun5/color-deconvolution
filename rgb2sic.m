function [ sic_coords ] = rgb2sic( rgb_coords, mu_s, sigma_s, rotation_matrix, options,  varargin) 
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

if strcmpi(normalizeflag,'on')
    F = @(x, mu, sigma) 1-exp(- max(x - mu,0)/2/sigma^2);
    Fmatrix = F(sqrt(sum(rotated_coordinates(2:3,:).^2,1)),mu_s,sigma_s); 
    sic_coords = repmat(Fmatrix,2,1).* rotated_coordinates(2:3,:)./...
        repmat(1 + sqrt(sum(rotated_coordinates(2:3,:).^2,1)),2,1);
else
    sic_coords = rotated_coordinates(2:3,:);
end


end

