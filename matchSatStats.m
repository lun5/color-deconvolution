%% match sat by statistics
target_sat_nonsat = target_sat(~target_indx_saturated);
source_sat_nonsat = source_sat(~source_indx_saturated);

% statistics of target
target_mean = mean(target_sat_nonsat);
target_std = std(target_sat_nonsat);
target_skewness = skewness(target_sat_nonsat);
target_kurtosis = kurtosis(target_sat_nonsat);

% statistics of source
source_mean = mean(source_sat_nonsat);
source_std = std(source_sat_nonsat);
source_skewness = skewness(source_sat_nonsat);
source_kurtosis = kurtosis(source_sat_nonsat);

% fix the mean and variance
source_sat_eq_nonsat = (source_sat_nonsat-source_mean) ...
    * target_std/source_std + target_mean;
% fix the range
source_sat_eq_nonsat(source_sat_eq_nonsat <= min(target_sat_nonsat)) = min(target_sat_nonsat);
source_sat_eq_nonsat(source_sat_eq_nonsat >= max(target_sat_nonsat)) = max(target_sat_nonsat);
% fix skewness 
[source_sat_eq_nonsat, ~] = modskew(source_sat_eq_nonsat,target_skewness);
% fix kurtosis
[source_sat_eq_nonsat, ~] = modkurt(source_sat_eq_nonsat,target_kurtosis);
% put it back together
source_sat_eq = source_sat;
source_sat_eq(~source_indx_saturated) = source_sat_eq_nonsat;

%% statistics after equalization
sprintf('Saturation Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_statistics %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    source_mean, source_std, min(source_sat), max(source_sat), source_skewness, source_kurtosis,...
    target_mean, target_std, min(target_sat_nonsat), max(target_sat_nonsat), target_skewness, target_kurtosis,...
    mean(source_sat_eq_nonsat), std(source_sat_eq_nonsat),...
    min(source_sat_eq), max(source_sat_eq),...
    skewness(source_sat_eq_nonsat), kurtosis(source_sat_eq_nonsat))
    
%% histograms 
A = {source_sat_nonsat,target_sat_nonsat,source_sat_eq_nonsat};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');
