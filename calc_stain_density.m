%% this function convert from RGB image into OD density 
% using Ruifrok and Johnston
% Luong Nguyen July 25, 2016

function [stain_density, hem_im_rgb, eosin_im_rgb] = calc_stain_density(im_rgb, opts)
   [h,w,c] = size(im_rgb);
   if max(max(im_rgb)) > 1
       im_rgb = double(im_rgb)./255;
   end
   optical_density = - reshape(log(im_rgb+eps),[h*w, c])';
   stain_matrix = [0.644,0.717,0.267; 0.093,0.954,0.283]'; 
   % the stain matrix is already in od values
   stain_density = pinv(stain_matrix)*optical_density;
   stain_density(stain_density < 0) = 0;
   stain_density(stain_density == Inf) = max(stain_density(stain_density < Inf));
   % reconstruct the hematoxylin/eosin images
   if ~exist('opts','var'); opts = 0; end
   if opts == 1
       hem_im = stain_matrix*cat(1,stain_density(1,:),zeros(1,h*w));
       hem_im_rgb = uint8(reshape(exp(-hem_im')*255,[h,w,c]));
       eosin_im = stain_matrix*cat(1,zeros(1,h*w),stain_density(2,:));
       eosin_im_rgb = uint8(reshape(exp(-eosin_im')*255,[h,w,c]));
   else
       hem_im_rgb = [];
       eosin_im_rgb = [];
   end
end
