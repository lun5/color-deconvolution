im1 = reshape(source_oppCol_eq(1,:),[source_xsize, source_ysize]);
%im1 =rgb2gray(source_im);
im1 = double(im1);
%rect = getrect; im1 = imcrop(im1,rect);
[dx1,dy1] = gradient(im1);
[x1 y1] = meshgrid(1:source_xsize,1:source_ysize);
%[x1 y1] = meshgrid(1:size(im1,2),1:size(im1,1));
u1 = dx1;
v1 = dy1;
figure; imshow(im1,[]);
hold on
quiver(x1,y1,u1,v1)
gradsize1 = (u1.^2 + v1.^2).^0.5;
figure; imagesc(gradsize1);axis square
colorbar;

im2 = reshape(source_oppCol_eq(2,:),[source_xsize, source_ysize]);
%rect = getrect; im2 = imcrop(im1,rect);
%[x1 y1] = meshgrid(1:source_xsize,1:source_ysize);
[x2 y2] = meshgrid(1:size(im1,2),1:size(im1,1));
[dx2,dy2] = gradient(im2);
[x2 y2] = meshgrid(1:source_xsize,1:source_ysize);
u2 = dx2;
v2 = dy2;
figure;imshow(im2,[]);
hold on
figure; quiver(x2,y2,u2,v2)
gradsize2 = (u2.^2 + v2.^2).^0.5;
figure; imagesc(gradsize2);axis square
colorbar;

sumgrad = gradsize1 + gradsize2;
figure; imshow(sumgrad,[]);
figure; imagesc(sumgrad);axis square
colorbar;

[Gmag, Gdir] = imgradient(source_im(:,:,2),'prewitt');

figure; imshowpair(Gmag, Gdir, 'montage');
title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Prewitt method')
axis off;

im = rgb2gray(imread('coloredChips.png'));
[x, y] = meshgrid(1:size(im,2),1:size(im,1));
im = double(im);
[dx,dy] = gradient(im);
u = dx;
v = dy;
figure;imshow(im,[]);
hold on
quiver(x,y,u,v)
gradsize = (u.^2 + v.^2).^0.5;
figure; imagesc(gradsize);

