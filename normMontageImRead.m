function [norm_images, avg_time] = normMontageImRead()

clf;
targ_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512\Target\';
targlist = dir(fullfile(targ_dir,'*.tif')); targlist = {targlist.name}';
source_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512';
sourcelist = dir(fullfile(source_dir,'*.tif')); sourcelist = {sourcelist.name}';

for index=1:1155
    im_dir = strcat('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\', 'Luong');
    imlist = dir(fullfile(im_dir,'*.tif')); imlist = {imlist.name}';
    imname = imlist{index};
    filename = strsplit(imname, '-');
    targetname = strcat(filename{1}, '.tif');
    list_methods = cell(4,1);
    method_names = {'Original Source', 'Original Target', 'Luong', 'Macenko', 'Reinhard',...
           'Khan', 'Vahadane', 'Vahadane Fast'};
    for i=3:6
        im_dir = strcat('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\', method_names{i});
        imlist = dir(fullfile(im_dir,'*.tif')); imlist = {imlist.name}';
        imname = imlist{index};
        list_methods{i} = imread(fullfile(im_dir,imname));
    end
   
    sourcename = filename{2};
    savedir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\Montage\';
    savename = fullfile(savedir, imname);
    if ~exist(savename,'file')
        source = imread(fullfile(source_dir,sourcename));
        list_methods{1} = source;
        target = imread(fullfile(targ_dir,targetname));
        list_methods{2} = target;
        
        montage = figure(1); montage.Position = [150, 100, 1600, 900];
        for k=1:length(list_methods)
            subplot_tight(1, 4, k, .02); imshow(list_methods{k}); title(method_names{k});
        end
        printfile('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\Montage\', montage, targetname, sourcename);
    end
end
end