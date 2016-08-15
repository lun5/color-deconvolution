%function runNorm_1

% im_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512';
% targ_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512\Target\';
% print_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\15 Targets Tiles 512';
code_dir = 'D:\Documents\Tiles_Norm\codes\';addpath(genpath(code_dir));
spams_dir = 'D:\Documents\Tiles_Norm\codes\spams-matlab';
cd(spams_dir); start_spams;addpath(genpath(spams_dir));

cd(code_dir);
data_dir = 'D:\Documents\Tiles_Norm\';
im_dir = fullfile(data_dir,'Tiles_2k');
targ_dir = fullfile(data_dir,'Target_2k');
print_dir = fullfile(data_dir,'normalized_images_2k');
if ~exist(print_dir,'dir');mkdir(print_dir);end;

imlist = dir(fullfile(im_dir,'*.tif'));
imlist = {imlist.name}';
targlist = dir(fullfile(targ_dir,'*.tif'));
targlist = {targlist.name}';
%time_matrix = cell(length(targlist),1);
time_matrix = cell(length(targlist),length(imlist));
%time_matrix = cell(5,length(imlist));

method_names = {'Luong','Macenko','Reinhard','Khan','Vahadane','VahadaneFast'};
for i = 1:length(method_names)
   outdir = fullfile(print_dir,method_names{i});
   if ~exist(outdir,'dir'); mkdir(outdir);end;
end

if ~exist(targ_dir,'dir'); mkdir(targ_dir);end;
for i=1:length(targlist)
    targetname = targlist{i};
    target = imread(fullfile(targ_dir,targetname));
    time_matrix_target = cell(length(imlist),1);
    tt = tic;
    for j=1:length(imlist)
        imname = imlist{j};imname = lower(imname);
        %time_methods = zeros(1,6);
        time_methods = zeros(1,1);
        if strcmp(targetname, imname)== 0
            source_im = imread(fullfile(im_dir,imname));
            output_name = [targetname(1:end-4) '-' imname];
%             if ~exist(fullfile(print_dir, method_names{1},output_name),'file');
%                 start = tic; normLuong = NormLuong(source_im, target);
%                 time_methods(1) = toc(start);
%                 imwrite(normLuong,fullfile(print_dir, method_names{1},output_name));
%             end
%             if ~exist(fullfile(print_dir, method_names{2},output_name),'file');
%                 start = tic; normMacenko = NormMacenko(source_im, target);
%                 time_methods(2) = toc(start);
%                 imwrite(normMacenko,fullfile(print_dir, method_names{2},output_name));
%             end
%             if ~exist(fullfile(print_dir, method_names{3},output_name),'file');
%                 start = tic; normReinhard = im2uint8(NormReinhard(source_im, target));
%                 time_methods(3) = toc(start);
%                 imwrite(normReinhard,fullfile(print_dir, method_names{3},output_name));
%             end
            if ~exist(fullfile(print_dir, method_names{4},output_name),'file');
                start = tic; normSCDLeeds = NormSCDLeeds(source_im, target);
                time_methods = toc(start);
                imwrite(normSCDLeeds,fullfile(print_dir, method_names{4},output_name));
            end
%             start = tic; normVahadane = SNMFnorm(source_im, target);
%             time_methods(5) = toc(start);
%             printfile(fullfile(print_dir, method_names{5}), normVahadane, targetname, imname);
%             start = tic; normVahadaneFast = Demo_colornorm(source_im, target);
%             time_methods(6) = toc(start);
%             printfile(fullfile(print_dir, method_names{6}), normVahadaneFast, targetname, imname);
        end
        %time_matrix{i, j} = time_methods;
        time_matrix_target{j} = time_methods;
    end
    time_matrix{i} = time_matrix_target;
    fprintf('Done with target num %d: %s in %.2f seconds\n',i, targetname,toc(tt));
end
save('total_time_Aug_8_Khan.mat', 'time_matrix');
%end