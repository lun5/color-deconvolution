function [ rgb_image ] = od2rgb( od_image,xsize, ysize )
%od2rgb convert from optical density space back to rgb space
% INPUT:optical density od_image, matrix 3x#pixels, sizes of original image
% OUPUT:rgb_image of size xsize x ysize x 3
normRGB = power(10,-od_image); % there are 3% of these numbers > 1
%normRGB = exp(-od_image);
RGB = uint8(normRGB*255);
r = reshape(RGB(1,:),[xsize, ysize]);
g = reshape(RGB(2,:),[xsize, ysize]);
b = reshape(RGB(3,:),[xsize, ysize]);

rgb_image = cat(3,r,g,b);

end

