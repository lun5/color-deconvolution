function [ rgb_image ] = raw2rgb( raw_image )
%raw2rgb convert
% INPUT: raw_image: m x n x 3 pixel image into
% OUTPUT: rgb_image of 3xM (M = num_pixels)
num_pixels = numel(raw_image(:,:,1));
Rvec = reshape(raw_image(:,:,1),1,num_pixels); 
Gvec = reshape(raw_image(:,:,2),1,num_pixels);  
Bvec = reshape(raw_image(:,:,3),1,num_pixels); 
rgb_image = double(cat(1,Rvec,Gvec,Bvec)); % R,G,B vector

end

