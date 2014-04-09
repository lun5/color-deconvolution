function [ color_deconv_rgb ] = stainvec2rgb( stain_vec, sat_vec, xsize, ysize )
%stainvec2rgb calculates the rgb image of a stain given
%INPUT: stain_vector (3x1), sat_vec saturation vector (1x num_pixels)
%   od_image optical density of original image, size of original image
%OUTPUT: original image deconvoled into image with stain_vec
num_pixels = xsize * ysize;
od_stain_image = repmat(sat_vec,3,1).*repmat(stain_vec,1,num_pixels);

color_deconv_rgb = od2rgb(od_stain_image,xsize, ysize);
end

