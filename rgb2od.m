function [ od_image ] = rgb2od( rgb_image )
%rgb2od convert from rgb space to optical density space 
% INPUT:rgb_image of size 3 x (xsize x ysize), the flatten rgb image
% OUPUT:od_image of size 3x number of pixels
rgb_image(rgb_image == 0) = 1;
RGBnorm = rgb_image./255;
od_image = -log10(RGBnorm); % convert image to optical density values
%od_image = -log(RGBnorm); % convert image to optical density values

end

