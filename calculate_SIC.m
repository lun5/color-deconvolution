% calculate_SIC
% Thesis by Elmaraghi 2003
rotation_matrix = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3); ...
    1/sqrt(6) 1/sqrt(6) -2/sqrt(6); ...
    -1/sqrt(2) 1/sqrt(2) 0];
rotated_coordinates = rotation_matrix*rgb_image;

mu_s = 4; sigma_s = 4; % transition control parameters
F = @(x, mu, sigma) 1-exp(- max(x - mu,0)/2/sigma^2);
sic_coord = F(normc(rotated_coordinate(1:2,:)),mu_s,sigma_s).* ...
    rotated_coordinate(1:2,:)./(1 + normc(rotated_coordinate(1:2,:)));
