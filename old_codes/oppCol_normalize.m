% Normalization script
% Luong Nguyen
% July 12, 2014

workdir = '/Users/lun5/Research/color_deconvolution'; 
datadir = fullfile(workdir, '20x_images','TIFF');

prompt = 'Please enter target image name ';
target_imname = input(prompt,'s');
while isempty(target_imname)
    target_imname = input(prompt,'s');
end

prompt = 'Please enter source image name ';
source_imname = input(prompt,'s');
while isempty(source_imname)
    source_imname = input(prompt,'s');
end

% get the names of the files
target_im = imread(fullfile(datadir,target_imname));
target_rgb = raw2rgb(target_im);figure;imshow(target_im);
source_im = imread(fullfile(datadir,source_imname));
source_rgb = raw2rgb(source_im);figure;imshow(source_im);

% get the training data
color_space = 'oppCol';
[ training_data, labels, rotation_matrix ] = import_training_data( color_space );
% convert the image to oppCol space
options = struct('Normalize','on');
mu_s = 4; sigma_s = 2; % values for normalization
[ target_oppCol, ~] = rgb2oppCol( target_rgb, mu_s, sigma_s, rotation_matrix, options); 
[ source_oppCol, ~] = rgb2oppCol( source_rgb, mu_s, sigma_s, rotation_matrix, options); 

% get the first component of the oppCol space
target_firstComp = rotation_matrix(1,:)*target_rgb;
source_firstComp = rotation_matrix(1,:)*source_rgb;
figure; hist(target_firstComp,100);
figure; hist(source_firstComp,100);
% calculate the angle to first component/second component
% the two are orthogonal so we can just look at one of them
nstep = 10;

figure;
scatter(target_oppCol(1,1:nstep:end),target_oppCol(2,1:nstep:end),20,target_rgb(:,1:nstep:end)'./255,'filled');
hold on;
hold off;
axis([-1 1 -1 1]);
axis equal

figure;
scatter(source_oppCol(1,1:nstep:end),source_oppCol(2,1:nstep:end),20,source_rgb(:,1:nstep:end)'./255,'filled');
hold on;
hold off;
axis([-1 1 -1 1]);
axis equal

target_angles = atan2(target_oppCol(2,:),target_oppCol(1,:));
source_angles = atan2(source_oppCol(2,:),source_oppCol(1,:));
figure; hist(target_angles,100);
figure; hist(source_angles,100);

% how to convert from source to target?
% fit distribution then convert one distribution to the other
% or fit between quartile1, median, quartile3
source_stats = [mean(source_firstComp) prctile(source_firstComp,5) prctile(source_firstComp,25) ...
    prctile(source_firstComp,50) prctile(source_firstComp,75) prctile(source_firstComp,95) ...
    mean(source_angles) prctile(source_angles, 5) prctile(source_angles, 25) prctile(source_angles, 50) ...
    prctile(source_angles, 75) prctile(source_angles, 95) ];

target_stats = [mean(target_firstComp) prctile(target_firstComp,5) prctile(target_firstComp,25) ...
    prctile(target_firstComp,50) prctile(target_firstComp,75) prctile(target_firstComp,95) ...
    mean(target_angles) prctile(target_angles, 5) prctile(target_angles, 25) prctile(target_angles, 50) ...
    prctile(target_angles, 75) prctile(target_angles, 95) ];

% fit statistic of target image as a function of source image
cs = spline(source_stats,target_stats);
xx = 0:10:300;
figure; plot(source_stats,target_stats,'o',xx,ppval(cs,xx),'-');
axis([0 300 0 300]);

% conver source image into the target image
figure; plot(source_stats,target_stats,'o');
axis([0 300 0 300]);
% show the results.

% what if we only fit the linear equation to this
p = polyfit(source_stats,target_stats,1); p
yfit = polyval(p,source_stats);
yresid = target_stats - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(target_stats) - 1) * var(target_stats);
rsq = 1 - SSresid/SStotal; % pretty decent rsq
% rsq =
% 
%     0.9943

% should try the version where I convert to OD first