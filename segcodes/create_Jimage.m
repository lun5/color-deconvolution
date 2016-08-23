%% function Jimage = create_Jimage(im, win_size)
% Luong Nguyen 10/7/15
% in the desperate hope that this might improve the performance. i'm blue
% today

function Jimage = create_Jimage(im, win_size)

% create a disk filter
h = fspecial('disk',win_size);
h(h>0) = 1;

if win_size >= 17
    ratio = floor((win_size+1)/9);
    h1 = zeros(size(h)); h1(1:ratio:end,1:ratio:end) = 1;
    h = h.*h1;
end

% pad the image with 0;
padded_im = padarray(im,[win_size win_size]);
% for each pixel, identify the surrounding disk
[nr, nc] = size(im);
if nr < win_size || nc < win_size
    error('window size has the be smaller than image sizes');
end
Jimage = zeros(nr,nc);
parfor i = 1:nr 
    for j = 1:nc
        center_window = padded_im(i:i+2*win_size, j:j+2*win_size);
        center_window = center_window.*h;
        Jimage(i,j) = calculate_jscore(center_window);
    end
end

end
