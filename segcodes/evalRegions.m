function [ff_score, bb_score] = evalRegions(groundTruth,result)

% identify background and foreground in the ground truth
gto = groundTruth{1,1};

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
logical_cells = cellfun(cellfind('stroma'),gto.names);
len= length(logical_cells)+1;
stromaLabel=len;
if sum(logical_cells) > 0
    stromaLabel = find(logical_cells==1);
end
whiteLabel = len;
logical_cells = cellfun(cellfind('white'),gto.names);
if sum(logical_cells) > 0
    whiteLabel = find(logical_cells==1);
end
gt = gto.Segmentation;
gt(gt==0)=stromaLabel;
noback = gt~=stromaLabel &  gt~=whiteLabel;
noback(gto.Boundaries)=0;
L= bwlabeln(noback);
foreground_CC = bwconncomp(L); % objects in the foreground
bg_gt = ~noback;
if whiteLabel == len && stromaLabel == len
    bg_gt(:) = 0;
end
background_CC = bwconncomp(bg_gt); % connected components in the background
bg_indx_list_gt = find(bg_gt);
%% identify background and foreground in the segmentation
num_seg_result = length(unique(result)); % max(results(:))
seg_pixel_list = cell(num_seg_result,1);
result_bg = zeros(num_seg_result,1);
for i = 1:num_seg_result
    mask = result == i;
    seg_pixel_list{i} = find(mask);
    intersect_bg_gt = length(intersect(seg_pixel_list{i}, bg_indx_list_gt));
    %if intersect_bg_gt/sqrt(length(seg_pixel_list{i})*length(bg_indx_list_gt)) >= 0.6
    if intersect_bg_gt/min(length(seg_pixel_list{i}),length(bg_indx_list_gt)) >= 0.6
        result_bg(i) = 1;
    end
end
% no background
if sum(bg_gt(:)) == 0
    result_bg(:) = 0;
end 
%% calculate foreground-foregroundscore
num_fg_objs_gt = length(foreground_CC.PixelIdxList);
fg_scores = zeros(num_fg_objs_gt,1);
matched_fg_id_list = zeros(num_fg_objs_gt,1);

num_elts_fg = zeros(num_fg_objs_gt,1);
for i = 1:num_fg_objs_gt
    num_elts_fg(i) = length(foreground_CC.PixelIdxList{i});
end

[~,sort_indx] = sort(num_elts_fg,'descend');
foreground_CC.PixelIdxList = foreground_CC.PixelIdxList(sort_indx);

for i = 1:num_fg_objs_gt
   max_overlap = 0;
   max_seg_id = 0;
   for j = 1:num_seg_result
      if result_bg(j) == 1 % check if the segment is background
          continue;
      end
      % eliminate segments already matched with the groundtruth foreground
      if i > 1 && matched_fg_id_list(i-1) == j
          continue;
      end
      result_overlap = intersect(seg_pixel_list{j},foreground_CC.PixelIdxList{i});
      if length(result_overlap) > max_overlap
          max_overlap = length(result_overlap);
          max_seg_id = j;
      end      
   end   
   matched_fg_id_list(i) = max_seg_id;
   if max_overlap > 0
       % weighted by the size of the groundtruth
       fg_scores(i) = max_overlap/length(union(seg_pixel_list{max_seg_id},...
           foreground_CC.PixelIdxList{i}))*length(foreground_CC.PixelIdxList{i});
   end
   %fg_scores(i) = max_overlap;
end
%ff_score = sum(fg_scores)/sum(noback(:));
ff_score = sum(fg_scores(fg_scores> 0))/sum(noback(:));

% %% calculate background-background score
% bg_scores = zeros(num_seg_result,1);
% for j = 1:num_seg_result
%     if result_bg(j) == 0;
%         continue;
%     end
%     result_overlap = intersect(seg_pixel_list{j},bg_indx_list_gt);
%     bg_scores(j) = length(result_overlap)/length(bg_indx_list_gt);
% end
% bg_scores(result_bg == 0) = [];
% bb_score = median(bg_scores);

%% calculate background-background score
if sum(bg_gt(:)) == 0 % if there is no background
    bb_score = -1;
    return;
end

num_bg_objs_gt = length(background_CC.PixelIdxList);
bg_scores = zeros(num_bg_objs_gt,1);
matched_bg_id_list = zeros(num_bg_objs_gt,1);

num_elts_bg = zeros(num_bg_objs_gt,1);
for i = 1:num_bg_objs_gt
    num_elts_bg(i) = length(background_CC.PixelIdxList{i});
end

[~,sort_indx] = sort(num_elts_bg,'descend');
background_CC.PixelIdxList = background_CC.PixelIdxList(sort_indx);

for i = 1:num_bg_objs_gt
   max_overlap = 0;
   max_seg_id = 0;
   for j = 1:num_seg_result
      if result_bg(j) == 0 % check if the segment is foreground
          continue;
      end
      % eliminate segments already matched with the groundtruth background
      if i > 1 && matched_bg_id_list(i-1) == j
          continue;
      end
      result_overlap = intersect(seg_pixel_list{j},background_CC.PixelIdxList{i});
      if length(result_overlap) > max_overlap
          max_overlap = length(result_overlap);
          max_seg_id = j;
      end      
   end   
   matched_bg_id_list(i) = max_seg_id;
   if max_overlap > 0
       bg_scores(i) = max_overlap/length(union(seg_pixel_list{max_seg_id},...
           background_CC.PixelIdxList{i}))*length(background_CC.PixelIdxList{i}); %bg_scores(i) = max_overlap;
   end
end
%bb_score = sum(bg_scores)/sum(bg_gt(:));
if sum(bg_scores) == 0
    bb_score = 0;
else
    bb_score = sum(bg_scores(bg_scores > 0))/sum(bg_gt(:));
end
