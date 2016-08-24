function [seg] = object_proposal_by_type( img_dir,adjDela_dir, mask_dir, imname, param_string, ...
    obj_type, max_num_neighbors, num_comps )
% compute the top connected components from the object maps in Burak's
% output. 

% INPUT:
% -IMG_DIR: image and object directory
% -imname: image name
% -param_string: of form _se1_minNuc3_minStr5_minLum9_adjDela, with info
% about minimum nuclear, stroma, and lumen size
% -obj_type: type of object (N-1, S-2, L-3) for finding connected component
% -num_neighbors: number of neighbors to identify the median distances
% between objects (default 15)
% -num_comp: number of top connected component to report (default 10)
% -plot_flag: 1 then plot and save the output (default 0)

% OUTPUT:
% -seg: segmentation mask of the superpixel type
% Luong Nguyen 3/13/2016
if nargin < 5
    error('Please input: image directory, adj directory, mask dir, image name, parameter string');
elseif nargin < 6
    obj_type = 1;
    max_num_neighbors = 25;
    num_comps =20;
elseif nargin < 7
    num_comps =20;
    max_num_neighbors = 25;
elseif nargin < 8
    num_comps =20;
end

% read in image, mask, adjacency map
info = imfinfo(fullfile(img_dir,[imname '.tif']));
nrow = info.Height; ncol = info.Width; 
adj_map = dlmread(fullfile(adjDela_dir,[imname param_string '_adjDela']),',',0,0);
%circle_map = dlmread(fullfile(IMG_DIR,'circle_map',[imname param_string '_circle_map']),',',1,0);
%voronoi_map = dlmread(fullfile(IMG_DIR,'voronois',[imname param_string '_voronoiMap']),',',0,0);
%num_total_objects = adj_map(1,1);
split_name = strsplit(imname,'-');
if obj_type == 3
    whole_mask = dlmread(fullfile(mask_dir,[split_name{2} '_mask']));
    whole_mask = whole_mask(1:4:end,1:4:end);
end
adj_map = adj_map(2:end,:);
obj_types = adj_map(:,3);
obj_coords = adj_map(:,[5 4]);
obj_radii = sqrt(adj_map(:,2)./pi);


% find the centers and radii of the object of certain type
indx_type = obj_types == obj_type;
% centers = obj_coords(indx_type,:);
% radii = obj_radii(indx_type);
num_superpx = sum(indx_type);

adj_obj = adj_map(indx_type,:); 
adj_obj_coords = adj_obj(:,[5 4]);
adj_obj_radii = sqrt(adj_obj(:,2)./pi);
adj_obj_area = adj_obj(:,2);
% slow way can be fast if we calculate the index first
indx_1 = cell(num_superpx,1);
indx_2 = cell(num_superpx,1);
dist_values = cell(num_superpx,1);
first_indx = find(obj_types == obj_type,1,'first');
for i = 1:num_superpx
   id1 = adj_obj(i,1);
   num_neighbors_input = adj_obj(i,6);% number of neighbors
   % do not connect small white objects
   if obj_type == 3 && adj_obj_area(i) < 10 %100 for 2k x 2k, 35 for 512 x512
       continue;
   end
   if num_neighbors_input > 0
       % neighbors of neighbors, neighbor-ception
       %n_n = adj_map(adj_obj(i,7:(6+num_neighbors_input)),[1,7:end]);
       %n_n = n_n(:)'; n_n(n_n == 0) = []; 
       %neighbor_indx = [adj_obj(i,7:(6+num_neighbors_input)), n_n];
       %neighbor_indx(neighbor_indx == id1) = [];
       neighbor_indx = adj_obj(i,7:(6+num_neighbors_input));
       
       % check if the neighbor is actually nuclei/whatever
       neighbor_types = obj_types(neighbor_indx);
       neighbors_of_same_type = neighbor_indx(neighbor_types == obj_type);
       num_neighbors_of_same_type = length(neighbors_of_same_type);
       % limit the number of neighbor to be fewer than max_num_neighbors
       %num_neighbors_of_same_type = min(max_num_neighbors, num_neighbors_of_same_type);
       if num_neighbors_of_same_type > 0
           distances = zeros(1, num_neighbors_of_same_type);
           for j = 1:num_neighbors_of_same_type
               % CHANGE FOR FEATURES
               distances(j) = norm(obj_coords(id1,:)-obj_coords(neighbors_of_same_type(j),:));
           end
           % limit the number of neighbor to be fewer than max_num_neighbors
           num_neighbors_of_same_type = min(max_num_neighbors, num_neighbors_of_same_type);       
           [sort_dists, sort_indx] = sort(distances);
           indx_1{i} = repmat(id1-first_indx+1,[1, num_neighbors_of_same_type]);
           indx_2{i} = neighbors_of_same_type(sort_indx(1:num_neighbors_of_same_type))-first_indx+1;
           dist_values{i} = sort_dists(1:num_neighbors_of_same_type); %repmat(distances,[1, 2]);
       end
       
   end
end
%T = toc; fprintf('Indexing done in %.2f seconds\n',T);
% Elapsed time is 0.411493 seconds.
indx_1 = cat(2,indx_1{:});
indx_2 = cat(2,indx_2{:});
dist_values = cat(2,dist_values{:});
%proportion_features = cat(1, proportion_features{:});
% CHANGE FOR FEATURES
med_dist = median(unique(dist_values)); 
% save the distance values
%save(fullfile(outdir,[imname '_distances_to_nn.mat']),'dist_values');
if obj_type == 3
    indx = dist_values >= 3*med_dist;%prctile(dist_values,75);
else
    indx = dist_values >= med_dist;%prctile(dist_values,75);
end
sigma = 1.5*med_dist;
dist_new = exp(- dist_values.^2./(2*sigma^2));
dist_new(indx) = 0;
similarities = sparse([indx_1,indx_2],[indx_2 indx_1],[dist_new dist_new],num_superpx,num_superpx);
[~,C] = graphconncomp(similarities);

components_areas = zeros(max(C),1);
for i = 1:length(components_areas)
    indx_cl = C == i;
    components_areas(i) = sum(obj_radii(indx_cl).^2);
end

[~,indx] = sort(components_areas,'descend');
% make sure the area is > 500 for 2k x 2k and 35 for 512 x 512
min_area = 31;
num_comps = min(sum(components_areas > min_area),num_comps);
top_centers = cell(num_comps,1);
top_radii = cell(num_comps,1);
indx_sp_cl = cell(num_comps,1); % index of superpixel
for i = 1:num_comps
   indx_cl = C == indx(i);
   top_centers{i} = adj_obj_coords(indx_cl,:);
   top_radii{i} = adj_obj_radii(indx_cl);
   indx_sp_cl{i} = find(indx_cl > 0) + first_indx - 1;
end

%if plot_flag
%I = imread(fullfile(img_dir,'images',[imname '.tif']));
%redc = I(:,:,1); greenc = I(:,:,2); bluec = I(:,:,3);
%colors = distinguishable_colors(num_comps)*255;
%se = strel('disk',7,4);
seg = zeros(nrow,ncol);
for i = 1:num_comps
%     voronoi_mask = zeros(nrow,ncol);
%     for j = 1: length(indx_sp_cl{i})
%         voronoi_mask = voronoi_mask + (voronoi_map == indx_sp_cl{i}(j));
%     end
%     mask = imfill(voronoi_mask,8,'hole')>0;
    x = top_centers{i}(:,1); y = top_centers{i}(:,2);
    if length(x) ==1
        ang=0:0.01:2*pi;
        xp=floor(top_radii{i}*cos(ang)'+x);
        yp=floor(top_radii{i}*sin(ang)'+y);
        k = boundary(xp,yp);
        mask = poly2mask(xp(k),yp(k),nrow, ncol);
    elseif length(x) <= 10
        extended_x = [x; min(x+top_radii{i},nrow); max(x - top_radii{i},0)];
        extended_y = [y; min(y+top_radii{i},nrow); max(y - top_radii{i},0)];
        k = boundary(extended_x,extended_y,0.5);
        mask = poly2mask(extended_x(k),extended_y(k),nrow, ncol);
    else
        k = boundary(x,y,0.5);
        mask =  poly2mask(x(k),y(k),nrow, ncol);
    end
    % avoid overlapping
    %mask_obj = zeros(nrow,ncol);
    %indx_obj = sub2ind(size(mask_obj),y,x);
    %mask_obj(indx_obj) = 1;
    %mask_obj = imdilate(mask_obj,se)>0;
    areas(i) = (max(y) - min(y))*(max(x) - min(x));
    area_mask = sum(mask(:));
    %if area_mask > 8000
    %    bdry = seg2bdry(double(mask),'imageSize');
    %    bdry = imdilate(bdry,se);
    %    bdry = logical(bdry);
    %    redc(bdry) = colors(i,1);
    %    redc(mask_obj) = colors(i,1);
    %    greenc(bdry) = colors(i,2);
    %    greenc(mask_obj) = colors(i,2);
    %    bluec(bdry) = colors(i,3);
    %    bluec(mask_obj) = colors(i,3);
    %end
    if sum(seg(mask)>0)/sum(mask(:)) < 0.3
        seg(mask) = i;
    else
        mask_diff = ((seg - mask) == -1);
        seg(mask_diff) = i;
    end
    %figure; imshow(label2rgb(seg));
end
%I2 = cat(3,redc, greenc, bluec);
if obj_type == 3
    seg(whole_mask == 0) = num_comps + 1;
end

%seg = seg + 1;
%segs{1} = seg(1:4:end, 1:4:end); % downsize the input for eval
%parsave(fullfile(outdir,[imname '.mat']),segs);
%end

end
