% Logistic regression model for classification of purple and pink pixels
% CV accuracy on training data
% Luong Nguyen
% July 2, 2014

% need to classify with RGB 
% Work directory and training data directory
close all;
workdir = '/Users/lun5/Research/color_deconvolution'; 
trainingdir = fullfile(workdir, 'results', '140625');
% Read in the labels and the input
training_purple = load([trainingdir filesep 'training_purple.mat'],'training_data_purple');
training_pink = load([trainingdir filesep 'training_pink.mat'],'training_data_pink');
X_purple_rgb = training_purple.training_data_purple;
X_pink_rgb = training_pink.training_data_pink;

label_purple = ones(size(X_purple_rgb,2),1);
label_pink = zeros(size(X_pink_rgb,2),1);
labels = [label_purple; label_pink];
%Convert from RGB --> OppCol space
options = struct('Normalize','on');
X_rgb = [X_purple_rgb'; X_pink_rgb'];
% training_data = [X_purple_rgb(:,1:3000) X_pink_rgb(:,1:9000)];
% [U,D,V] = svd(training_data,0);
% rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one
% mu_s = 4; sigma_s = 2; % values for normalization
% rgb_image = [X_purple_rgb X_pink_rgb];
% %[ oppCol_image, indx_sat] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options); 
% [X_purple_oppCol,~] = rgb2oppCol(X_purple_rgb, mu_s, sigma_s, rotation_matrix, options);
% [X_pink_oppCol,~] = rgb2oppCol(X_pink_rgb, mu_s, sigma_s, rotation_matrix, options);
% X_oppCol = [X_purple_oppCol';X_pink_oppCol'];
%train the classifier
%Use CV to estimate the classification error
num_folds =10; % number of folds for CV
Acc = zeros(num_folds,2);    
indices = crossvalind('Kfold',labels,num_folds);
for i = 1:num_folds
        test = (indices == i); train_indx = ~test;
        % estimate the parameters using logistic regression
        % b's are parameters
        % dev: deviance is a generalization of the residual sum of squares
        % stats: a struct with multiple fields
        %% train the model
        X_train = X_rgb(train_indx,:);
        [b,dev,stats] = glmfit(X_train,labels(train_indx),'binomial','link','logit');
        labels_prob_train = glmval(b,X_train,'logit');
        labels_pred_train = labels_prob_train > 0.5;
        acc_train = sum(labels_pred_train == labels(train_indx))./length(labels_pred_train).*100;
        % accuracy on training 
        % plot the decision boundary, where the probability of belonging to
        % each class is 0.5
        indx_bound_train = abs(labels_prob_train - 0.5) < 7*1e-2;
        h1 = figure;
        subplot(2,2,1); 
        % plot the actual classes of training data % the index sound kinda
        % wrong?
        plot(X_train(labels(train_indx) == 0,1),X_train(labels(train_indx) == 0,2),'r.', 'MarkerSize',10);
        hold on;
        plot(X_train(labels(train_indx) == 1,1),X_train(labels(train_indx) == 1,2),'b.', 'MarkerSize',10);
        hold off;
        %axis([-1 1 -1 1]);axis equal;

        % plot the predicted classes of training data
        subplot(2,2,2);% figure;
        plot(X_train(labels_pred_train == 0,1),X_train(labels_pred_train == 0,2),'r.', 'MarkerSize',10);
        hold on;
        plot(X_train(labels_pred_train == 1,1),X_train(labels_pred_train == 1,2),'b.', 'MarkerSize',10);
        plot(X_train(indx_bound_train,1),X_train(indx_bound_train,2),'k-', 'MarkerSize',30,'LineWidth',3); % dec boundary
        hold off;
        %axis([-1 1 -1 1]);axis equal;
        
        %% test the model
        X_test = X_rgb(test,:);
        labels_prob_test = glmval(b,X_test,'logit')'; % probability of label being purple (=1)
        labels_pred_test = labels_prob_test > 0.5; % threshold by 0.5; prob < 0.5 is pink, else purple
        acc_test = sum(labels_pred_test' == labels(test))./length(labels_pred_test).*100;
        % plot the actual classes of test data
        subplot(2,2,3);
        plot(X_test(labels(test) == 0,1),X_test(labels(test) == 0,2),'r.', 'MarkerSize',10);
        hold on;
        plot(X_test(labels(test) == 1,1),X_test(labels(test) == 1,2),'b.', 'MarkerSize',10);
        hold off;
        %axis([-1 1 -1 1]); axis equal;

        % plot the decision boundary, where the probability of belonging to
        % each class is 0.5
        indx_bound_test = abs(labels_prob_test - 0.5) < 9*1e-2;
        % plot the predicted classes of test data        
        subplot(2,2,4); %figure;
        plot(X_test(labels_pred_test == 0,1),X_test(labels_pred_test == 0,2),'r.', 'MarkerSize',10);
        hold on;
        plot(X_test(labels_pred_test == 1,1),X_test(labels_pred_test == 1,2),'b.', 'MarkerSize',10);
        plot(X_test(indx_bound_test,1),X_test(indx_bound_test,2),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
        hold off;
        %axis([-1 1 -1 1]);axis equal;

        %% Accuracies
        Acc(i,1) = acc_train; Acc(i,2) = acc_test; 
end

%plot the scatter plot of 
%apply the classifier on images
nstep = 100;
figure;
scatter3(X_rgb(1:nstep:end,1),X_rgb(1:nstep:end,2),X_rgb(1:nstep:end,3),20,X_rgb(1:nstep:end,:)./255,'filled');
%axis([-1 1 -1 1])

% how does pink and purple look like?
% purple
figure;
scatter3(X_purple_rgb(1,:), X_purple_rgb(2,:), X_purple_rgb(3,:), 20,X_purple_rgb'./255,'filled'); 
%axis([-1 1 -1 1])

% pink
figure;
scatter3(X_pink_rgb(1,:), X_pink_rgb(2,:), X_pink_rgb(3,:),20,X_pink_rgb'./255,'filled');
%axis([-1 1 -1 1])