%% function jscore = calculate_jscore(class_map)
% Luong Nguyen 10/7/15
% INPUT: local class map
% OUTPUT: jscore 

function jscore = calculate_jscore(class_map)
classes = unique(class_map(class_map > 0));
[I_all, J_all] = find(class_map > 0);
mean_all_classes = [mean(I_all), mean(J_all)];
S_T = sum(sum((([I_all, J_all] - repmat(mean_all_classes,[length(I_all),1])).^2),2));

num_classes = length(classes);
mean_classes = zeros(num_classes,2);
S_W = 0;
for i = 1:num_classes
   [I_class, J_class] = find(class_map == classes(i));
   mean_classes(i,:) = [mean(I_class), mean(J_class)];
   S_W = S_W + sum(sum((([I_class, J_class] - ...
       repmat(mean_classes(i,:),[length(I_class),1])).^2),2)); 
end

jscore = (S_T - S_W)/S_W;
end
