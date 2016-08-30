%% Luong Nguyen 08/10/2016
% Test the Normalization for the MITOS dataset
data_dir = 'M:\MITOS\';
dir_list = {'training_aperio','testing_aperio','training_hamamatsu','testing_hamamatsu'};
method_names = {'Luong','Macenko','Khan','Reinhard','Vahadane','VahadaneFast'};

for mm = 3%4:length(method_names)
met = method_names{mm};

t1 = tic;
fprintf('Start with method %s\n',met);
for dd = 1:2
    source_list = dir(fullfile(data_dir,dir_list{dd},'*.tar'));
    source_list = {source_list.name}';
    
    for ss = 1:length(source_list)
        source_name = source_list{ss}(1:end-4);
        t2 = tic;
        resol_levels = dir(fullfile(data_dir,dir_list{dd},source_name,'frames'));
        resol_levels = {resol_levels.name}';
        for rr = 1:length(resol_levels)
            if strcmp(resol_levels{rr},'.') || strcmp(resol_levels{rr},'..')
                continue;
            end
            source_dir = fullfile(data_dir,dir_list{dd},source_name,...
                'frames',resol_levels{rr});
            target_dir = fullfile(data_dir,dir_list{dd+2},...
                ['H' source_name(2:end)],'frames',resol_levels{rr});
            %s2t_dir = fullfile(data_dir, dir_list{dd},source_name,...
            %    'normalized_frames',resol_levels{rr});
            %s2t_dir = fullfile(data_dir, dir_list{dd},source_name,...
            %    'Macenko_normalized_frames',resol_levels{rr});
            s2t_dir = fullfile(data_dir, dir_list{dd},source_name,...
                [met '_normalized_frames'],resol_levels{rr});
            
            if ~exist(s2t_dir,'dir'); mkdir(s2t_dir); end;
            %t2s_dir = fullfile(data_dir, dir_list{dd+2},...
            %    ['H' source_name(2:end)],'normalized_frames',resol_levels{rr});
            %t2s_dir = fullfile(data_dir, dir_list{dd+2},...
            %    ['H' source_name(2:end)],'Macenko_normalized_frames',resol_levels{rr});
            t2s_dir = fullfile(data_dir, dir_list{dd+2},...
                ['H' source_name(2:end)],[met '_normalized_frames'],resol_levels{rr});
            if ~exist(t2s_dir,'dir'); mkdir(t2s_dir); end;
            imlist = dir(fullfile(source_dir,'*.tiff'));
            imlist = {imlist.name}';
            for ii = 1:length(imlist)
                imname = imlist{ii};
                source_image = imread(fullfile(source_dir,imname));
                target_image = imread(fullfile(target_dir,['H' imname(2:end)]));
                if ~exist(fullfile(s2t_dir,imname),'file')
                    %s2t_im = NormLuong(source_image, target_image);
                    %s2t_im = NormMacenko(source_image, target_image);
                    if strcmp(met,'Luong'); s2t_im = NormLuong(source_image, target_image); end
                    if strcmp(met,'Macenko'); s2t_im = NormMacenko(source_image, target_image); end
                    if strcmp(met,'Khan');  s2t_im = NormSCDLeeds(source_image, target_image); end
                    if strcmp(met,'Vahadane'); s2t_im =  SNMFnorm(source_image, target_image); end
                    if strcmp(met,'VahadaneFast'); s2t_im = Demo_colornorm(source_image, target_image); end
                    if strcmp(met,'Reinhard'); s2t_im = im2uint8(NormReinhard(source_image, target_image)); end
                    imwrite(s2t_im,fullfile(s2t_dir,imname));
                end
                if ~exist(fullfile(t2s_dir,['H' imname(2:end)]),'dir')
                    %t2s_im = NormLuong(target_image, source_image);
                    %t2s_im = NormMacenko(target_image, source_image);
                    if strcmp(met,'Luong'); t2s_im = NormLuong(target_image, source_image); end
                    if strcmp(met,'Macenko'); t2s_im = NormMacenko(target_image, source_image); end
                    if strcmp(met,'Khan');  t2s_im = NormSCDLeeds(target_image, source_image); end
                    if strcmp(met,'Vahadane'); t2s_im =  SNMFnorm(target_image, source_image); end
                    if strcmp(met,'VahadaneFast'); t2s_im = Demo_colornorm(target_image, source_image); end
                    if strcmp(met,'Reinhard'); t2s_im = im2uint8(NormReinhard(target_image, source_image)); end
                    imwrite(t2s_im,fullfile(t2s_dir,['H' imname(2:end)]));
                end
            end
            clear t2s_im; clear s2t_im;
        end
        tt2 = toc(t2);
        fprintf('Done with dir %s in %.2f seconds\n',source_name,tt2);
    end
    tt1 = toc(t1);
    fprintf('Done with method %s in %.2f seconds\n',met,tt1);
end

end
