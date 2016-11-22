% script to process the results
% Read in the table

% Calculate for each pair of source and target
% Method from 1 to 6

count = 0;
for tt = 1:2
    for ss = 1:10
        for m1 = 1:5
            for m2 = (m1+1) : 6
                count = count +1;
                source_num{count} = ss;
                target_num{count} = tt;
                m1_num{count} = m1;
                m2_num{count} = m2;
                results{count} = randi(3) - 2;
            end
        end
    end
end

T = table(source_num',target_num',m1_num',m2_num',results');
T.Properties.VariableNames = {'source_num','target_num','m1_num','m2_num', 'Results'};
writetable(T,'C:\Users\luong_nguyen\Box Sync\ColorNormalizationPaper\Tiles_512 Validation Data\GUI\fake_outputs.txt')        

T = readtable('C:\Users\luong_nguyen\Box Sync\ColorNormalizationPaper\Tiles_512 Validation Data\GUI\fake_outputs.txt');

% test the logic
count_inconst = 0;
for tt = 1:2
    for ss = 1:10
       indx_pair = (T.source_num == ss) & (T.target_num == tt);
       pairwise_comp_mat = sparse([T.m1_num(indx_pair) T.m2_num(indx_pair)],...
          [T.m2_num(indx_pair) T.m1_num(indx_pair)], [T.Results(indx_pair) T.Results(indx_pair)]);
       pairwise_comp_mat = full(pairwise_comp_mat);
       scores = zeros(6,1);
       for mm = 1:6
           %scores(mm) = sum(pairwise_comp_mat(:,mm) == 1) + 0.5*sum(pairwise_comp_mat(:,mm) == 0);
           scores(mm) = sum(T.m1_num(indx_pair) == mm & T.Results(indx_pair) == 1) + ...
               sum(T.m2_num(indx_pair) == 1 & T.Results(indx_pair) == -1) + ...
               0.5* sum(T.m1_num(indx_pair) == mm & T.Results(indx_pair) == 0); 
       end
       [sorted_scores, sorted_indx] = sort(scores,'descend');
       for m1 = 1:5
           for m2 = (m1+1):6
               m1_num = sorted_indx(m1);
               m2_num = sorted_indx(m2);
               comp_result = pairwise_comp_mat(m1_num,m2_num);
               if (comp_result == -1) || (comp_result == 0  && ...
                       sorted_scores(m1_num) ~= sorted_scores(m2_num))
                  fprintf('method %d score is %.2f, method %d score is %.2f\n', ...
                      m1_num, sorted_scores(m1_num),m2_num,  sorted_scores(m2_num));
                  count_inconst = count_inconst + 1; 
               end
           end
       end
    end 
end