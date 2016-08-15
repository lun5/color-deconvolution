source_dir = 'C:\Users\tam128\Box Sync\ColorNormalizationPaper\Tiles_512 Validation Data\GUI\Source';
target_dir = 'C:\Users\tam128\Box Sync\ColorNormalizationPaper\Tiles_512 Validation Data\GUI\Target';
source_normal_dir = 'C:\Users\tam128\Box Sync\ColorNormalizationPaper\Tiles_512 Validation Data\GUI\Normalized_Target';

method_names = {'Luong', 'Macenko', 'Reinhard',...
            'Khan', 'Vahadane', 'Vahadane Fast'};
        
source_im_list = dir(fullfile(source_dir,'*.tif'));
source_im_list = {source_im_list.name}';
target_im_list = dir(fullfile(target_dir,'*.tif'));
target_im_list = {target_im_list.name}';
count = 0;
for tt = 1:length(target_im_list)
    target_imname = target_im_list{tt};
    for ss = 1:length(source_im_list)
        source_imname = source_im_list{ss};
        for mm = 1:(length(method_names)-1)
            method1 = method_names{mm};
            for mm2 = (mm+1):length(method_names)               
                count = count + 1;
                source_names{count} = fullfile('Source',source_imname);
                target_names{count} = fullfile('Target',target_imname);
                method_1_names{count} = fullfile('Normalized_source', method1,[target_imname(1:end-4) '-' source_imname]);
                method_2_names{count} = fullfile('Normalized_source',method_names{mm2}, [target_imname(1:end-4) '-' source_imname]);
                method_1_num{count} = mm;
                method_2_num{count} = mm2;
                source_num{count} =ss;
                target_num{count}=tt;
            end
        end
    end
end

T = table(source_names', target_names', source_num', target_num', method_1_names', method_2_names', method_1_num', method_2_num');
T.Properties.VariableNames = {'Source', 'Target', 'source_num', 'target_num', 'm1','m2','m1_num','m2_num'};
writetable(T,'C:\Users\tam128\Box Sync\ColorNormalizationPaper\Tiles_512 Validation Data\GUI\GUI color norm pair.txt',...
    'Delimiter',',','WriteVariableNames',true);


