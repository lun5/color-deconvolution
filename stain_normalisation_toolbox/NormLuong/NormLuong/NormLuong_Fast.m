function  [normalized_image, init_params] = NormLuong_Fast(source_image, target_image,init_params)

if nargin < 3
    init_params = [];
elseif nargin <2
    error('Function needs at least 2 inputs: source and target images');
end

rotation_matrix = load('rotation_matrix_tp10-867-1.mat','rotation_matrix');
rotation_matrix = rotation_matrix.rotation_matrix;
% rotation_matrix = [0.6412    0.4903    0.5903;...
%     0.1340    0.6859   -0.7152;...
%     0.7556   -0.5376   -0.3741];

numClusters = 3; % only purple and pink this time
opts_mixture.noise = 1;
%% recursive call to function for very large images
source_image_size = size(source_image);
%nrows = source_image_size(1); ncols = source_image_size(2);
if max(source_image_size)>1e4
    % cut off 500 each side, ONLY IN ADH
    cut_off = 800;
    normalized_image = source_image;
    source_image = source_image(cut_off:end-cut_off,cut_off:end-cut_off,:);   
    nrows = size(source_image,1); ncols = size(source_image,2);
    nsteps = 1e4;
    rep_rows = floor(nrows/nsteps); res_rows = nrows - nsteps*rep_rows;
    rep_cols = floor(ncols/nsteps); res_cols = ncols - nsteps*rep_cols;
    if res_rows > 0;
        dimR = [repmat(nsteps,1,rep_rows), res_rows];
    else
        dimR = repmat(nsteps,1,rep_rows);
    end
    if res_cols > 0
        dimC = [repmat(nsteps,1,rep_cols), res_cols];
    else
        dimC = [repmat(nsteps,1,rep_cols)];
    end
    tiles_rgb = mat2cell(source_image, dimR,dimC,3);
    tiles_pct_white = cell2mat(cellfun(@(x) ...
        sum(sum(mean(double(x)./255,3) > 0.9))/(size(x,1)*size(x,2)),...
        tiles_rgb, 'UniformOutput',false));
    % find the initial parameters from the most dense area
    nim2 = cell(numel(tiles_rgb),1);
    [~,indx] = min(tiles_pct_white(:));
    rgb = tiles_rgb{indx};
    %T = tic;
    [nim2{indx}, init_params] = NormLuong_Fast(rgb, target_image, init_params);   
    %fprintf('Done with getting init param in %.2f\n',toc(T));
    non_white_indx = find(tiles_pct_white < 0.93);
    %white_indx = ~ismember(1:numel(tiles_rgb), non_white_indx);
    %nim = cell(size(tiles_rgb));
    
    for ii = 1:numel(tiles_rgb)
        rgb = tiles_rgb{ii};
        %T = tic;
        %Datam = double(reshape(rgb,size(rgb,1)*size(rgb,2),size(rgb,3)));
        % how to avoid black spots
        if ismember(ii, non_white_indx) && ii ~=indx %&&( sum(mean(Datam))<650 &&  sum(var(Datam))> 500)
            [nim2{ii}, init_params] = NormLuong_Fast(rgb, target_image, init_params);
        elseif ii~=indx
            nim2{ii} = rgb;
        end
        %fprintf('Done with tile %d in %.2f\n',ii, toc(T));
        %tiles_rgb{non_white_indx(ii)} = [];
    end
    normalized_image(cut_off:end-cut_off,cut_off:end-cut_off,:) = cell2mat(reshape(nim2,size(tiles_rgb)));
%     tic;
%     nim_non_white = cellfun(@(x) NormLuong_Fast(x,target_image),...
%        tiles_rgb(non_white_indx), 'UniformOutput',false);
%     toc
%     nim{non_white_indx} = nim_non_white; clear nim_non_white;
%     nim{white_indx} = tiles_rgb{white_indx};
    clear tiles_rgb
    %normalized_image = cell2mat(nim);
    return;    
end

%% calculate features
for n = 1:2
    if n == 1
        im_rgb = double(source_image)./255;
        clear source_image
    else
        im_rgb = double(target_image)./255;
    end   
    which_features = {'hue opp', 'brightness opp','saturation opp'};
    nrows = size(im_rgb,1); ncols = size(im_rgb,2);
    X = reshape(im_rgb,[nrows*ncols,3]);
    rotated_coordinates = rotation_matrix*X';
    
    %% cluster using von Mises
    theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
    %theta_o = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
    %theta = theta_o + 3.34;
    %theta(theta>pi) = theta(theta>pi) - 2*pi;
    theta_im = reshape(theta,nrows,ncols);
    sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
    brightness = rotated_coordinates(1,:);
    % von Mises
    X_cart = [cos(theta); sin(theta)]';
    %num_pixels = nrows*ncols;
    
    if n == 1
        f_maps_source = {reshape(theta,[nrows,ncols]), reshape(brightness,[nrows,ncols]),...
            reshape(sat,[nrows,ncols])};
    else
        f_maps_target = {reshape(theta,[nrows,ncols]), reshape(brightness,[nrows,ncols]),...
            reshape(sat,[nrows,ncols])};
    end
    if nrows < 3000 &&  ncols < 3000
        [ mu_hat_polar,~, kappa_hat,posterior_probs, ~,~] =...
            moVM(X_cart,numClusters,opts_mixture, init_params);
    elseif max(nrows,ncols) <= 1e4
        if isempty(init_params)
            nsteps = 500;
            rep_rows = floor(nrows/nsteps); res_rows = nrows - nsteps*rep_rows;
            rep_cols = floor(ncols/nsteps); res_cols = ncols - nsteps*rep_cols;
            if res_rows > 0;
                dimR = [repmat(nsteps,1,rep_rows), res_rows];
            else
                dimR = repmat(nsteps,1,rep_rows);
            end
            if res_cols > 0
                dimC = [repmat(nsteps,1,rep_cols), res_cols];
            else
                dimC = [repmat(nsteps,1,rep_cols)];
            end
            tiles_rgb = mat2cell(im_rgb, dimR,dimC,[3]);
            tiles_pct_white = cell2mat(cellfun(@(x) ...
                sum(sum(mean(double(x),3) > 0.9))/(size(x,1)*size(x,2)),...
                tiles_rgb, 'UniformOutput',false));
            non_white_indx = find(tiles_pct_white < 0.6);
            tiles_theta = mat2cell(theta_im, dimR,dimC);
            numpass = 0;
            ss_indx = non_white_indx(randperm(length(non_white_indx)));
            num_samples = min(length(non_white_indx),8);
            max_num_while = min(2,ceil(length(ss_indx)/num_samples));
            num_while = 0;
            %conv = 0;
            mu_hat = []; kappa_hat = [];
            while (numpass == 0) || (num_while < max_num_while) ||...
                    isempty(mu_hat) || isempty(kappa_hat)
                num_while = num_while + 1;
                sample_indx = ss_indx(((num_while-1)*num_samples + 1):num_while*num_samples);
                %randperm(numel(tiles_rgb),min(numel(tiles_rgb),8));
                [I,J] = ind2sub([size(tiles_rgb,1),size(tiles_rgb,2)],sample_indx);
                for tt = 1:length(I)
                    %tile_rgb = tiles_rgb{I(tt), J(tt)};
                    %tile_gray = mean(double(tile_rgb),3);
                    %if sum(tile_gray(:) > 0.9)/numel(tile_gray) < 0.5
                    numpass = numpass + 1;
                    tile_theta = tiles_theta{I(tt), J(tt)};
                    X_cart_tile = [cos(tile_theta(:))'; sin(tile_theta(:))']';
                    [ mu_hat,~,kappa_hat, ~, ~, conv] =...
                        moVM(X_cart_tile,numClusters,opts_mixture);
                    if conv
                        num_while = max_num_while;
                        break;
                    end
                    %end
                end
            end
            clear theta_im tiles_rgb tiles_theta dimR dimC rep_cols rep_rows tile_rgb tile_gray X_cart_tile sample_indx
            init_params = struct('theta_hat',mu_hat,'kappa_hat',kappa_hat);
        end
        opts_mixture.maxiter = 1;       
        [ mu_hat_polar,~,kappa_hat,posterior_probs,~,~] =...
                moVM( X_cart,numClusters,opts_mixture,init_params);
            
    end
                
    if n==1
        opts_matching.source_stats = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
            'posterior_probs',posterior_probs,'prior_probs',[]);
    else
        opts_matching.target_stats = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
            'posterior_probs',posterior_probs,'prior_probs',[]);
    end
end

f_maps_source_normalized = cell(1,3);
for feature_iter = 1:length(which_features)
    f_map_source_curr = f_maps_source{feature_iter};
    f_map_target_curr = f_maps_target{feature_iter};
    %% normalization
    f_map_normalized_curr = matchingMoments(f_map_source_curr, f_map_target_curr,which_features{feature_iter}, opts_matching);
    % change this
%     if feature_iter == 1
%        f_map_normalized_curr = f_map_normalized_curr - 3.34;
%        f_map_normalized_curr(f_map_normalized_curr < -pi) = f_map_normalized_curr(f_map_normalized_curr < -pi) + 2*pi;
%     end
    f_maps_source_normalized{feature_iter} = f_map_normalized_curr;
end

%% calculate c2, c3, recover rotated coordinate
% then calculate new rgb
source_rotated_eq = zeros(3,length(f_map_normalized_curr));
source_rotated_eq(1,:) = f_maps_source_normalized{2}; % brightness normalized/equalized
source_rotated_eq(2,:) = f_maps_source_normalized{3}.*cos(f_maps_source_normalized{1}); % c2
source_rotated_eq(3,:) = f_maps_source_normalized{3}.*sin(f_maps_source_normalized{1}); % c3
source_rgb_eq = rotation_matrix\source_rotated_eq; % RGB = Rot_matrix \ rotated coordinates
source_rgb_eq(source_rgb_eq < 0) = 0;
source_rgb_eq(source_rgb_eq > 1) = 1;
source_rgb_eq_uint8 = uint8(source_rgb_eq*255); 
normalized_image = reshape(source_rgb_eq_uint8', source_image_size);
end

