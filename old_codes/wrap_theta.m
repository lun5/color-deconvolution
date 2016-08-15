% wrap around theta
indx_neg = source_eq < -pi;
indx_pos = source_eq > pi;
source_eq(indx_neg) = source_eq(indx_neg)+ 2*pi;
source_eq(indx_pos) = source_eq(indx_pos)- 2*pi;
% figure;circ_plot(source_eq(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% 
% % statistics
% sprintf('Source theta: mean %.4f, std %.4f, min %.4f, max %.4f, skewness %.4f, kurtosis %.4f\n',...
%     circ_mean(source_eq,[],2),circ_std(source_eq,[],[],2), min(source_eq), max(source_eq), ...
%     circ_skewness(source_eq,[],2), circ_kurtosis(source_eq,[],2))
% 

