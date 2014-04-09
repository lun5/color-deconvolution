function [ od_image ] = rgb2od( raw_image )
%rgb2od convert from rgb space to optical density space 
% INPUT:raw_image of size xsize x ysize x 3
% OUPUT:od_image of size 3x number of pixels
RGB = raw2rgb(raw_image); % R,G,B vector
%RGBnorm = (RGB - repmat(min(RGB),3,1))./repmat(max(RGB) - min(RGB),3,1);
%RGBnorm = normc(RGB);
RGBnorm = RGB./255;
od_image = -log10(RGBnorm); % convert image to optical density values
%od_image = -log(RGBnorm); % convert image to optical density values

end

