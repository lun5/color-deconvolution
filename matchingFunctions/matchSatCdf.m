%% match saturation by cdf
target_sat_nonsat = target_sat(~target_indx_saturated);
source_sat_nonsat = source_sat(~source_indx_saturated);
% normalize the range of saturation to between 0 and 1
target_sat_norm = target_sat_nonsat./max(target_sat_nonsat);
source_sat_norm = source_sat_nonsat./max(source_sat_nonsat);
% count number of elements in each bin 
bincounts = histc(target_sat_norm,binranges);
%equalize source's brightness --> target's brightness in [0 1]
source_sat_norm_eq = histeq(source_sat_norm,bincounts);
% multiply by the range of target's brightness to complete normalization
source_sat_eq = source_sat;
source_sat_eq_nonsat = source_sat_norm_eq * max(target_sat);
source_sat_eq(~ source_indx_saturated) = source_sat_eq_nonsat;

% bincounts = histc(target_sat_norm,binranges);
% source_sat_norm_eq = histeq(source_sat_norm,bincounts);
% source_sat_eq = source_sat_norm_eq * max(target_sat);

%% statistics! 
sprintf('Saturation Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_cdf %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    mean(source_sat_nonsat), std(source_sat_nonsat), min(source_sat), max(source_sat),...
    skewness(source_sat_nonsat), kurtosis(source_sat_nonsat),...
    mean(target_sat_nonsat), std(target_sat_nonsat), min(target_sat), max(target_sat),...
    skewness(target_sat_nonsat), kurtosis(target_sat_nonsat),...
    mean(source_sat_eq_nonsat), std(source_sat_eq_nonsat),...
    min(source_sat_eq), max(source_sat_eq),...
    skewness(source_sat_eq_nonsat), kurtosis(source_sat_eq_nonsat))

%% histograms
A = {source_sat,target_sat,source_sat_eq};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');

% B = {source_sat_norm(~ source_indx_saturated),target_sat_norm(~target_indx_saturated),source_sat_norm_eq(~ source_indx_saturated)};
% %B = {source_sat_norm,target_sat_norm,source_sat_norm_eq};
% figure; nhist(B, 'legend',legendstr,'color',bincol,...
%     'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','minx',0,'maxx',1,'pdf');
