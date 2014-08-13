rng('default') % For reproducibility
X = [randn(100,2)+ones(100,2);...
     randn(100,2)-ones(100,2)];

opts = statset('Display','final');
[idx,ctrs] = kmeans(X,2,'Distance','city',...
              'Replicates',5,'Options',opts);

plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(ctrs(:,1),ctrs(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
hold off