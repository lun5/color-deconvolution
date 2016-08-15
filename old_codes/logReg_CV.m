% Read in the inputs
% read in the training data type
% convert the image to appropriate types
% then do the classification. 

% Logistic regression model for classification of purple and pink pixels
% CV accuracy on training data
% Luong Nguyen
% July 2, 2014

% need to classify with RGB 
% Work directory and training data directory
% for lab and hsv, the plot of separation should be a little different
function Acc = logReg_CV(num_folds, color_space, options, varargin) 

if nargin < 3
    options = [];
    if nargin < 2
        error('Need at least 2 inputs: number of folds for CV, and color space of input data');
    end
end

[ training_data, labels, rotation_matrix ] = import_training_data( color_space );

Acc = zeros(num_folds,2);    
indices = crossvalind('Kfold',labels,num_folds);

for i = 1:num_folds
        test = (indices == i); train_indx = ~test;
        % estimate the parameters using logistic regression
        % b's are parameters
        % dev: deviance is a generalization of the residual sum of squares
        % stats: a struct with multiple fields
        %% train the model
        X_train = training_data(train_indx,:);
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
        X_test = training_data(test,:);
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

% %plot the scatter plot of 
% %apply the classifier on images
% nstep = 100;
% figure;
% scatter3(X_rgb(1:nstep:end,1),X_rgb(1:nstep:end,2),X_rgb(1:nstep:end,3),20,X_rgb(1:nstep:end,:)./255,'filled');
% %axis([-1 1 -1 1])
% 
% % how does pink and purple look like?
% % purple
% figure;
% scatter3(X_purple_rgb(1,:), X_purple_rgb(2,:), X_purple_rgb(3,:), 20,X_purple_rgb'./255,'filled'); 
% %axis([-1 1 -1 1])
% 
% % pink
% figure;
% scatter3(X_pink_rgb(1,:), X_pink_rgb(2,:), X_pink_rgb(3,:),20,X_pink_rgb'./255,'filled');
% %axis([-1 1 -1 1])

end