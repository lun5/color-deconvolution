function  normalized_image = NormLuong_Fast(source_image, target_image)

rotation_matrix = load('rotation_matrix_tp10-867-1.mat','rotation_matrix');
rotation_matrix = rotation_matrix.rotation_matrix;
% rotation_matrix = [0.6412    0.4903    0.5903;...
%     0.1340    0.6859   -0.7152;...
%     0.7556   -0.5376   -0.3741];

numClusters = 3; % only purple and pink this time
opts_mixture.noise = 1;
%% calculate features
for n = 1:2
    if n == 1
        im_rgb = double(source_image)./255;
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
    num_pixels = nrows*ncols;
    
    if n == 1
        f_maps_source = {reshape(theta,[nrows,ncols]), reshape(brightness,[nrows,ncols]),...
            reshape(sat,[nrows,ncols])};
    else
        f_maps_target = {reshape(theta,[nrows,ncols]), reshape(brightness,[nrows,ncols]),...
            reshape(sat,[nrows,ncols])};
    end
    if nrows < 3000 &&  ncols < 3000
        [ mu_hat_polar,~, kappa_hat,posterior_probs, prior_probs] =...
            moVM(X_cart,numClusters,opts_mixture);
    else
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
        non_white_indx = find(tiles_pct_white < 0.5);
        tiles_theta = mat2cell(theta_im, dimR,dimC);
        numpass = 0;
        ss_indx = non_white_indx(randperm(length(non_white_indx)));
        num_samples = min(length(non_white_indx),8);
        max_num_while = min(5,ceil(length(ss_indx)/num_samples));
        num_while = 0;
        conv = 0;
        while (numpass == 0) && (num_while < max_num_while)
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
        mu_hat_polar_sample = mu_hat;
        kappa_hat_sample = kappa_hat;
        init_param = struct('theta_hat',mu_hat_polar_sample,'kappa_hat',kappa_hat_sample);
        if conv % does it matter if it converges?
            opts_mixture.maxiter = 1;
            %[ mu_hat_polar,~, kappa_hat,posterior_probs, prior_probs,~] =...
            %    moVM(X_cart,numClusters,opts_mixture, init_param);
        else
            opts_mixture.maxiter = 1;
        end
        
        nsteps = 15000;
        indx = 1:(nsteps*nsteps):num_pixels;
        if indx(end) < num_pixels; indx = [indx, num_pixels];end
        
        posterior_probs = cell(length(indx) -1,1);
        mu_hat_polar = cell(length(indx) -1,1);
        kappa_hat = cell(length(indx) -1,1);
        %prior_probs = cell(length(indx) -1,1);
        
        for jj = 1:(length(indx) - 1)
            %tic;
            if jj < (length(indx) - 1)
                xrange = indx(jj):(indx(jj+1)-1);
            else
                xrange = indx(jj):indx(jj+1);
            end
            [ mu_hat_polar{jj},~,kappa_hat{jj},posterior_probs{jj},~,~] =...
                moVM( X_cart(xrange,:),numClusters,opts_mixture,init_param);
            %fprintf('calculate the parameters at iter %d takes %.2f\n',jj, toc)
        end
    
        mu_hat_polar = mean(cat(1,mu_hat_polar{:}),1);
        kappa_hat = mean(cat(1,kappa_hat{:}),1);
        posterior_probs = cat(1,posterior_probs{:});
        %prior_probs = mean(cat(1,prior_probs{:}),1);
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
normalized_image = reshape(source_rgb_eq_uint8', size(source_image));
end
