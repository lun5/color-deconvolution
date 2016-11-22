function [ mu_hat_polar,mu_hat_cart, kappa_hat,posterior_probs, prior_probs, conv] = moVM_Fast(X_cart,k,opts,init_params)
    
    opts_default.maxiter = 100;
    opts_default.eps1 = 1e-2; % threshold for likelihood convergence
    opts_default.eps2 = 1e-2; % threshold for parameter convergence
    opts_default.noise = 1;
    if nargin < 4
        init_params = [];    
    elseif nargin <3
        opts = opts_default;
    elseif nargin <2
        error('Function needs at least 2 inputs: data, number of components');
    end
    
    if ~exist('opts','var')
        opts = opts_default;
    end
    
    if ~ isfield(opts,'maxiter')
        opts.maxiter = opts_default.maxiter;
    end
    
    if ~ isfield(opts,'eps1');
        opts.eps1 = opts_default.eps1;
    end

    if ~ isfield(opts,'eps2');
        opts.eps2 = opts_default.eps2;
    end
    
    if ~ isfield(opts,'noise');
        opts.noise = opts_default.noise;
    end

    % Set 'm' to the number of data points.
    numData = size(X_cart, 1);
    % set 'd' to the dimension
    d = size(X_cart,2); % for now it's just 2D
    X_polar = atan2(X_cart(:,2),X_cart(:,1));

    %% STEP 1: Initialization
    posterior_probs = zeros(numData,k + opts.noise);
    
    %% Loop through M-step and E-step until convergence
    % probability of each cluster -- prior
    prior_probs = ones(1,k+opts.noise)*(1/(k+ opts.noise));
    mu_hat_cart = zeros(d,k);
    mu_hat_polar = zeros(1,k);
    kappa_hat = zeros(1,k);
    
    %% randomly assign kappa
    % HOW TO AVOID HARD CODED HERE? Use peakfinder
    %mu_hat_polar(1) = -0.2; % nuclei purple
    %mu_hat_polar(2) = -1.7; % stroma pink 
    %mu_hat_polar(3) = 2.24; % lumen white
    
    %h = histogram(X_polar,100); set(gcf,'Visible','off');
    %values = h.Values; bin_centers = h.BinEdges + h.BinWidth/2;
    if isempty(init_params)
        kappa_hat(:) = 5;
        [values, edges] = histcounts(X_polar,100);
        bin_centers = edges + (edges(2) - edges(1))/2;
        bin_centers = bin_centers(1:end-1);
        pink_interval = [find(bin_centers > -2.5,1,'first') , find(bin_centers < -1,1,'last')];
        purple_interval = [find(bin_centers > -0.5,1,'first') find(bin_centers < 1.2,1,'last')];
        white_interval = find(bin_centers > 2,1,'first');
        pink_peak = peakfinder(values(pink_interval(1):pink_interval(2)));
        mu_hat_polar(2) = bin_centers(pink_interval(1) + pink_peak(1)); %stroma pink
        purple_peak = peakfinder(values(purple_interval(1)+1:purple_interval(2)));
        mu_hat_polar(1) = bin_centers(purple_interval(1) + purple_peak(1)); % purple nuclei
        white_peak = peakfinder(values(white_interval+1:end));
        mu_hat_polar(3) = bin_centers(white_interval + white_peak(1));% white
    else
        mu_hat_polar = init_params.theta_hat;
        kappa_hat = init_params.kappa_hat;
    end
    LLH = zeros(k + opts.noise, 1);
    for i = 1:k
        LLH(i) = sum(prior_probs(i) * ( - log(2*pi*besseli(0,kappa_hat(i)))+ ...
            kappa_hat(i)*cos(X_polar - mu_hat_polar(i))));
    end
    
    if opts.noise
       LLH(k+1) = (prior_probs(k+1)*log(1/(2*pi)))*length(X_polar);
    end
    % for deterministic annealing
    %mult = 1.015; kappa_threshold = 3; max_kappa = 100;
    noise_threshold = 0.05; kappa_max = 50;
for iter = 1: opts.maxiter
    
    mu_hat_old = mu_hat_polar;
    kappa_hat_old = kappa_hat;
    LLH_old = LLH;
    %% STEP 1: E-step
    for i = 1:k
        posterior_probs(:,i) = prior_probs(i)*circ_vmpdf(X_polar, mu_hat_polar(i), kappa_hat(i));
    end
    
    if opts.noise
        posterior_probs(:,k+1) = prior_probs(k+1)*repmat(1/(2*pi),numData,1);
    end
    
    posterior_probs = posterior_probs./repmat(sum(posterior_probs,2),1,k + opts.noise);

    %% STEP 2: M-step
    for i = 1:k
        prior_probs(i) = mean(posterior_probs(:,i));
        unnormalized_mean = sum(repmat(posterior_probs(:,i),1,d).*X_cart);
        mu_hat_cart(:,i) = unnormalized_mean'/norm(unnormalized_mean);
        %mu_hat_polar(i) = atan2(mu_hat_cart(2,i),mu_hat_cart(1,i));% FREEZE EVERYTHING
        rho = norm(unnormalized_mean)/sum(posterior_probs(:,i));
        rho_max = 0.99;
        if rho > rho_max % avoid singularity
            rho = rho_max - 0.001 + rand*0.001;
        end
        kappa_hat(i) = rho*(d - rho^2)/(1-rho^2);
        if kappa_hat(i) < 5; kappa_hat(i) = 5; end;
        %kappa_hat(i) = min(kappa_best,max_kappa) * (kappa_best >= threshold_kappa) + ...
        %    min(max_kappa, max(kappa_best, kappa_hat_old(i)*mult))*(kappa_best < threshold_kappa);
        LLH(i) = logLikelihood(X_polar, posterior_probs(:,i), mu_hat_polar(i), kappa_hat(i));
    end
        
    % rescale the uniform noise if it goes above noise threshold
    if opts.noise
        prior_probs(k+1) = 1 - sum(prior_probs(1:k));
        LLH(k+1) = sum(posterior_probs(:,k+1))+ (log(1/(2*pi))*length(X_polar));
    end
    
    if opts.noise && sum(prior_probs(1:k)) < 1 - noise_threshold 
        prior_probs(1:k) = prior_probs(1:k)*(1 - noise_threshold)/sum(prior_probs(1:k));
        prior_probs(k + 1) = noise_threshold;
    else
        %display('no need adjust noise level');
    end
                
    %% Stopping criteria
    % There is one very concentrated cluster of white. If the concentration
    % of this cluster is greater than 600 then we will stop the algorithm
%     if max(kappa_hat) > 300
%         break;
%         % or maybe use spkmeans at this point
%     end

    % else look at the change in posterior and parameters
    llh_change = norm(abs((LLH - LLH_old)./LLH_old));
    mu_change = norm(abs(mu_hat_polar - mu_hat_old));
    kappa_change = norm(abs(kappa_hat - kappa_hat_old));
    
    if llh_change  < opts.eps1 || (mu_change  < opts.eps2 && kappa_change  < opts.eps2)
        break;
    end
  
end

if iter < opts.maxiter
    %fprintf('The algorithm does not converge at maxiter %d',opts.maxiter);
    conv = 1;
else
    conv = 0;
end
  
end

function LLH = logLikelihood(ang, pij, mu_j, kappa_j)
    LLH = sum(pij.*(- log(2*pi*besseli(0,kappa_j))+ kappa_j*cos(ang - mu_j)));
end
