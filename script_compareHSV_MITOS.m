%% run script for MITOS

%% STEP 1: HSV calculation
data_dir = 'M:\MITOS';
dir_list = {'training_aperio','testing_aperio','training_hamamatsu','testing_hamamatsu'};
method_names = {'Luong','Macenko','Khan','Reinhard','Vahadane','VahadaneFast'};
PDollar_toolbox_path = 'C:\Users\luong_nguyen\Documents\GitHub\toolbox';
addpath(genpath(PDollar_toolbox_path));
metric_names = {'ks','chisq','emd','kl'};

for mm = 1:length(method_names)
    met = method_names{mm};    
    t1 = tic;
    fprintf('Start with method %s\n',met);
    measures_over_traintest = cell(2,1);

    for dd = 1:2
        source_list = dir(fullfile(data_dir,dir_list{dd},'*.tar'));
        t2 = tic;
        if ~exist(fullfile(data_dir,dir_list{dd},'HistDist_HSV'),'dir')
            mkdir(fullfile(data_dir,dir_list{dd},'HistDist_HSV'));
        end
        source_list = {source_list.name}';
        measures_over_sourcelist = cell(length(source_list),1); 
        for ss = 1:length(source_list)
            source_name = source_list{ss}(1:end-4);
            t3 = tic;
            resol_levels = dir(fullfile(data_dir,dir_list{dd},source_name,'frames'));
            resol_levels = {resol_levels.name}';
            measures_over_resol = cell(length(resol_levels),1);
            for rr = 1:length(resol_levels)
                if strcmp(resol_levels{rr},'.') || strcmp(resol_levels{rr},'..')...
                        || ~isempty(strfind(resol_levels{rr},'stain'))
                    continue;
                end
                t4 = tic;
                source_dir = fullfile(data_dir,dir_list{dd},source_name,...
                    'frames',resol_levels{rr});
                target_dir = fullfile(data_dir,dir_list{dd+2},...
                    ['H' source_name(2:end)],'frames',resol_levels{rr});
                s2t_dir = fullfile(data_dir, dir_list{dd},source_name,...
                    [met '_normalized_frames'],resol_levels{rr});
                
                if ~exist(s2t_dir,'dir'); mkdir(s2t_dir); end;
                t2s_dir = fullfile(data_dir, dir_list{dd+2},...
                    ['H' source_name(2:end)],[met '_normalized_frames'],resol_levels{rr});
                if ~exist(t2s_dir,'dir'); mkdir(t2s_dir); end;
                
                imlist = dir(fullfile(source_dir,'*.tiff'));
                imlist = {imlist.name}';
                measures_over_imlist = cell(length(imlist),1);
                source_images = cell(length(imlist),1);
                target_images = cell(length(imlist),1);
                parfor ii = 1:length(imlist)
                    imname = imlist{ii};
                    source_im = imread(fullfile(source_dir,imname));
                    target_im = imread(fullfile(target_dir,['H' imname(2:end)]));
                    source_images{ii} = fullfile(source_dir,imname);
                    target_images{ii} = fullfile(target_dir,['H' imname(2:end)]);
                    
                    s2t_im = imread(fullfile(s2t_dir,imname));
                    s2t_hsv = rgb2hsv(s2t_im); % normalized hsv
                    t2s_im = imread(fullfile(t2s_dir,['H' imname(2:end)]));
                    t2s_hsv = rgb2hsv(t2s_im);
                    source_hsv = rgb2hsv(source_im);
                    target_hsv = rgb2hsv(target_im);
                    metrics = zeros(2,length(metric_names)*(size(target_hsv,3)+1));% row 1: s2t, row 2: t2s
                    channel_hist_counts = cell(4,3);
                    for cc = 1:3
                        channel = source_hsv(:,:,cc);
                        [channel_hist_counts{1,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
                        channel = target_hsv(:,:,cc);
                        [channel_hist_counts{2,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
                        channel = s2t_hsv(:,:,cc);
                        [channel_hist_counts{3,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
                        channel = t2s_hsv(:,:,cc);
                        [channel_hist_counts{4,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
                    end
                    
                    for metr = 1:length(metric_names)
                        metric = metric_names{metr};
                        start_indx = (metr-1)*(size(target_hsv,3)+1) + 1;
                        if strcmp(metric,'ks')
                            for cc = 1:3
                                sc = source_hsv(:,:,cc); tc = target_hsv(:,:,cc);
                                s2t = s2t_hsv(:,:,cc); t2s = t2s_hsv(:,:,cc);
                                [~,~,metrics(1,start_indx + cc)] = kstest2(s2t(:), tc(:));
                                [~,~,metrics(2,start_indx + cc)] = kstest2(t2s(:), sc(:));
                            end
                        else
                            for cc = 1:3
                                % target vs s2t
                                metrics(1,start_indx+cc) = pdist2(channel_hist_counts{2,cc}, ...
                                    channel_hist_counts{3,cc},metric);
                                % source vs t2s
                                metrics(2,start_indx+cc) = pdist2(channel_hist_counts{1,cc}, ...
                                    channel_hist_counts{4,cc},metric);
                            end
                        end
                        metrics(1,start_indx) = sum(metrics(1,(start_indx+1):(start_indx+3)));
                        metrics(2,start_indx) = sum(metrics(2,(start_indx+1):(start_indx+3)));
                    end
                    measures_over_imlist{ii} = cat(2,metrics(1,:), metrics(2,:));
                end
                
                measures_over_imlist = num2cell(cat(1,measures_over_imlist{:}));
                clear t2s_im; clear s2t_im;
                measures_over_resol{rr} = measures_over_imlist;
                tt4 = toc(t4);
                fprintf('Done with resolution %s in %.2f s\n',resol_levels{rr},tt4);
            end
            measures_over_resol = num2cell(cat(1,measures_over_resol{:}));
            measures_over_sourcelist{ss} = measures_over_resol;
            tt3 = toc(t3);
            fprintf('Done with dir %s in %.2f seconds\n',source_name,tt3);
        end
        measures_over_sourcelist = num2cell(cat(1,measures_over_sourcelist{:}));
        measures_over_traintest{dd} = measures_over_sourcelist;
        tt2 = toc(t2);
        fprintf('Done with dir %s in %.2f seconds\n',dir_list{dd},tt2);
    end
    measures_over_traintest = num2cell(cat(1,measures_over_traintest{:}));
    T = cell2table(cat(2,source_images,target_images,measures_over_imlist));
    
    T.Properties.VariableNames = {'SourceImages','TargetImages',...
        's2t_ks','s2t_ks_h','s2t_ks_s','s2t_ks_v',...
        's2t_chisq','s2t_chisq_h','s2t_chisq_s','s2t_chisq_v',...
        's2t_emd','s2t_emd_h','s2t_emd_s','s2t_emd_v',...
        's2t_kl','s2t_kl_h','s2t_kl_s','s2t_kl_v',...
        't2s_ks','t2s_ks_h','t2s_ks_s','t2s_ks_v',...
        't2s_chisq','t2s_chisq_h','t2s_chisq_s','t2s_chisq_v',...
        't2s_emd','t2s_emd_h','t2s_emd_s','t2s_emd_v',...
        't2s_kl','t2s_kl_h','t2s_kl_s','t2s_kl_v'};
    writetable(T,fullfile(data_dir,dir_list{dd},'HistDist_HSV',[method_names{mm} '.txt']),'Delimiter',',');
    tt1 = toc(t1);
    fprintf('Done with method %s in %.2f seconds\n',met,tt1);
end

