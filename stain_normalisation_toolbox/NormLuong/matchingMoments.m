%% matching brightness by matching moments
%% what would we do with hue
% Initial commit: 5/27/2015
% revised: 3/3/2016
% Luong Nguyen
function source_eq = matchingMoments(source_im, target_im,which_feature, opts)
% we avoid touching white by having white_indx in the opts
if nargin < 4
    error('Please put options with indices of white pixels');
end
%indx_white_source = opts.source_stats.posterior_probs(:,3);
%indx_white_target = opts.target_stats.posterior_probs(:,3);

%% normalize the range of brightness to between 0 and 1
source_im = source_im(:); target_im = target_im(:);
%target_im_nw = target_im(~indx_white_target);
%source_im_nw = source_im(~indx_white_source);
% load from vM directores
[~, indx_membership_target] = max(opts.target_stats.posterior_probs,[],2);
[~, indx_membership_source] = max(opts.source_stats.posterior_probs,[],2); % 4 is the uniform noise
source_eq = source_im; % initialize the normalized source

for cl = 1:2 % only do purple and pink, LEAVE WHITE RED MASK ALONE(numClusters)% + 1)
    target_im_cl = target_im(indx_membership_target == cl)';
    source_im_cl = source_im(indx_membership_source == cl)';
    if (~isempty(source_im_cl)) && (~isempty(target_im_cl))
        if strcmp(which_feature,'hue opp')
            % statistics of target
            target_mean = opts.target_stats.mu_hat_polar(cl); %circ_mean(target_im_cl,[],2);
            [~,target_std] = circ_std(target_im_cl,[],[],2);
            
            % statistics of source
            source_mean = opts.source_stats.mu_hat_polar(cl); %circ_mean(source_im_cl,[],2);
            [~,source_std] = circ_std(source_im_cl,[],[],2);
            %% center to 0
            %source_eq_cl = source_im_nw - source_mean; %all of the image, not just the cluster of interest
            source_eq_cl = source_im_cl - source_mean; %just the cluster of interest
            source_eq_cl(source_eq_cl < - pi) = source_eq_cl(source_eq_cl < - pi)+ 2*pi;
            source_eq_cl(source_eq_cl > pi) = source_eq_cl(source_eq_cl > pi)- 2*pi;
            
            %% spread source to have same std as target
            source_eq_cl = source_eq_cl * target_std/source_std;
            
            %% fix mean to be equal target mean
            source_eq_cl = source_eq_cl + target_mean;
            source_eq_cl(source_eq_cl < - pi) = source_eq_cl(source_eq_cl < - pi)+ 2*pi;
            source_eq_cl(source_eq_cl > pi) = source_eq_cl(source_eq_cl > pi)- 2*pi;
            % To smooth out my choices?
            %source_eq_cl = source_eq_cl.*opts.source_stats.posterior_probs(:,cl);%??
        elseif strcmp(which_feature,'brightness opp')|| strcmp(which_feature,'saturation opp')
            % if we have empty cluster, copy the original image
            % statistics of target
            target_mean = mean(target_im_cl);
            target_std = std(target_im_cl);
            target_skewness = skewness(target_im_cl);
            target_kurtosis = kurtosis(target_im_cl);
            
            % statistics of source
            source_mean = mean(source_im_cl);
            source_std = std(source_im_cl);
            %source_skewness = skewness(source_im_cl);
            %source_kurtosis = kurtosis(source_im_cl);
            
            % fix the mean and variance
            source_eq_cl = (source_im_cl - source_mean)*target_std/source_std + target_mean;
            % fix the range
            source_eq_cl(source_eq_cl <= min(target_im_cl)) = min(target_im_cl);
            source_eq_cl(source_eq_cl >= max(target_im_cl)) = max(target_im_cl);
            % fix skewness
            %if skewness(source_eq_cl) ~= target_skewness
            %    [source_eq_cl, ~] = modskew(source_eq_cl,target_skewness);
            %end
            % fix kurtosis
            %if kurtosis(source_eq_cl) ~= target_kurtosis
            %    [source_eq_cl, ~] = modkurt(source_eq_cl,target_kurtosis);
            %end
            %source_eq_skewness = skewness(source_eq_nw);
            %source_eq_kurtosis = kurtosis(source_eq_nw);
        end
        %% put it back together
        source_eq(indx_membership_source == cl) = source_eq_cl;
    else
        source_eq(indx_membership_source == cl) = source_im_cl;
    end
end

end

% 
% if strcmp(which_feature,'brightness opp')|| strcmp(which_feature,'saturation opp')
%     % statistics of target
%     target_mean = mean(target_im_nw);
%     target_std = std(target_im_nw);
%     target_skewness = skewness(target_im_nw);
%     target_kurtosis = kurtosis(target_im_nw);
%     
%     % statistics of source
%     source_mean = mean(source_im_nw);
%     source_std = std(source_im_nw);
%     source_skewness = skewness(source_im_nw);
%     source_kurtosis = kurtosis(source_im_nw);
%     
%     % fix the mean and variance
%     source_eq_nw = (source_im_nw - source_mean)*target_std/source_std + target_mean;
%     % fix the range
%     source_eq_nw(source_eq_nw <= min(target_im_nw)) = min(target_im_nw);
%     source_eq_nw(source_eq_nw >= max(target_im_nw)) = max(target_im_nw);
%     % fix skewness
%     [source_eq_nw, ~] = modskew(source_eq_nw,target_skewness);
%     % fix kurtosis
%     [source_eq_nw, ~] = modkurt(source_eq_nw,target_kurtosis);
%     source_eq_skewness = skewness(source_eq_nw);
%     source_eq_kurtosis = kurtosis(source_eq_nw);
%     %% statistics after equalization
% %     sprintf('Feature %s \n',which_feature)
% %     sprintf('Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nMatched Source %.4f %.4f %.4f %.4f %.4f %.4f\n',...
% %         source_mean, source_std, min(source_im_nw),...
% %         max(source_im_nw), source_skewness, source_kurtosis,...
% %         target_mean, target_std, min(target_im_nw),...
% %         max(target_im_nw), target_skewness, target_kurtosis,...
% %         mean(source_eq_nw), std(source_eq_nw),...
% %         min(source_eq_nw), max(source_eq_nw),...
% %         source_eq_skewness, source_eq_kurtosis)
%     source_eq = source_im; source_eq(~indx_white_source) = source_eq_nw;
% elseif strcmp(which_feature,'hue opp')
%     % mixture of univariate von Mises distributions
%     %X_cart_target = [cos(target_im) sin(target_im)];
%     %% Call the function
%     %         numClusters = 3; opts_mixture.noise = 1;
%     %         [ mu_hat_polar_target,~, ~,posterior_probs_target, ~] =...
%     %            moVM(X_cart_target,numClusters, opts_mixture);
%     %         [~, indx_membership_target] = max(posterior_probs_target,[],2); % 4 is the uniform noise
%     %
%     %         X_cart_source = [cos(source_im) sin(source_im)];
%     % load these from precalculated
%     [~, indx_membership_target] = max(opts.target_stats.posterior_probs,[],2);
%     %kappa_hat_target = opts.kappa_hat;
%     %% Call the function
%     %[ mu_hat_polar_source,~, ~,posterior_probs_source, ~] =...
%     %   moVM(X_cart_source,numClusters, opts_mixture); % no noise version
%     source_eq = source_im;
%     [~, indx_membership_source] = max(opts.source_stats.posterior_probs,[],2); % 4 is the uniform noise
%     %% normalize for each cluster
%     for cl = 1:2 % only do purple and pink, LEAVE WHITE ALONE(numClusters)% + 1)
%         target_im_cl = target_im(indx_membership_target == cl)';
%         source_im_cl = source_im(indx_membership_source == cl)';
%         % statistics of target
%         target_mean = opts.target_stats.mu_hat_polar(cl); %circ_mean(target_im_cl,[],2);
%         [~,target_std] = circ_std(target_im_cl,[],[],2);
%         
%         % statistics of source
%         source_mean = opts.source_stats.mu_hat_polar(cl); %circ_mean(source_im_cl,[],2);
%         [~,source_std] = circ_std(source_im_cl,[],[],2);
%         %% center to 0
%         %source_eq_cl = source_im_nw - source_mean; %all of the image, not just the cluster of interest
%         source_eq_cl = source_im_cl - source_mean; %just the cluster of interest
%         source_eq_cl(source_eq_cl < - pi) = source_eq_cl(source_eq_cl < - pi)+ 2*pi;
%         source_eq_cl(source_eq_cl > pi) = source_eq_cl(source_eq_cl > pi)- 2*pi;
%         
%         %% spread source to have same std as target
%         source_eq_cl = source_eq_cl * target_std/source_std;
%         
%         %% fix mean to be equal target mean
%         source_eq_cl = source_eq_cl + target_mean;
%         source_eq_cl(source_eq_cl < - pi) = source_eq_cl(source_eq_cl < - pi)+ 2*pi;
%         source_eq_cl(source_eq_cl > pi) = source_eq_cl(source_eq_cl > pi)- 2*pi;
%         % Do I take the mean hear or something? To smooth out my choices?
%         %source_eq_cl = source_eq_cl.*opts.source_stats.posterior_probs(:,cl);%??
%         %% put it back together
%         source_eq(indx_membership_source == cl) = source_eq_cl;
%         
%         %sprintf('Feature %s Cluster number %d\n',which_feature, cl)
%         %sprintf('Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nMatched Source %.4f %.4f %.4f %.4f %.4f %.4f\n',...
%         %source_mean, source_std, min(source_im),...
%         %max(source_im), source_skewness, source_kurtosis,...
%         %target_mean, target_std, min(target_im),...
%         %max(target_im), target_skewness, target_kurtosis,...
%         %circ_mean(source_eq_cl'), circ_std(source_eq_cl'),...
%         %min(source_eq_cl), max(source_eq_cl),...
%         %circ_skewness(source_eq_cl'),circ_kurtosis(source_eq_cl'))
%     end
%     %         if opts_mixture.noise
%     %             source_eq_nw = source_eq_nw + source_im_nw.*posterior_probs_source(:,numClusters + 1);
%     %         end
%     source_eq(source_eq < - pi) = source_eq(source_eq < - pi)+ 2*pi;
%     source_eq(source_eq > pi) = source_eq(source_eq > pi)- 2*pi;
%     % do not touch the white or the noise :)
% end
% %% histogram
% A = {source_im_nw,target_im_nw,source_eq(~indx_white_source)};
% legendstr = {'source','target','matched source'};
% bincol = [0.8,0.8,0.8];
% figure; nhist(A,'legend',legendstr,'color',bincol,'xlabel',[which_feature ' values'],...
%     'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');
% end