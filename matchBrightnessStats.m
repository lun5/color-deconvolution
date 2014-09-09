%% match brightness by statistics
target_brightness_nonsat = target_brightness(~target_indx_saturated);
source_brightness_nonsat = source_brightness(~source_indx_saturated);

% statistics of target
target_mean = mean(target_brightness_nonsat);
target_std = std(target_brightness_nonsat);
target_skewness = skewness(target_brightness_nonsat);
target_kurtosis = kurtosis(target_brightness_nonsat);

% statistics of source
source_mean = mean(source_brightness_nonsat);
source_std = std(source_brightness_nonsat);
source_skewness = skewness(source_brightness_nonsat);
source_kurtosis = kurtosis(source_brightness_nonsat);

% fix the mean and variance
source_brightness_eq_nonsat = (source_brightness_nonsat-source_mean) ...
    * target_std/source_std + target_mean;
% fix the range
source_brightness_eq_nonsat(source_brightness_eq_nonsat <= min(target_brightness_nonsat)) = min(target_brightness_nonsat);
source_brightness_eq_nonsat(source_brightness_eq_nonsat >= max(target_brightness_nonsat)) = max(target_brightness_nonsat);
% fix skewness 
[source_brightness_eq_nonsat, ~] = modskew(source_brightness_eq_nonsat,target_skewness);
% fix kurtosis
[source_brightness_eq_nonsat, ~] = modkurt(source_brightness_eq_nonsat,target_kurtosis);
% put it back together
source_brightness_eq = source_brightness;
source_brightness_eq(~source_indx_saturated) = source_brightness_eq_nonsat;

% statistics after equalization
sprintf('Brightness Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_statistics %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    source_mean, source_std, min(source_brightness), max(source_brightness), source_skewness, source_kurtosis,...
    target_mean, target_std, min(target_brightness), max(target_brightness), target_skewness, target_kurtosis,...
    mean(source_brightness_eq_nonsat), std(source_brightness_eq_nonsat),...
    min(source_brightness_eq), max(source_brightness_eq),...
    skewness(source_brightness_eq_nonsat), kurtosis(source_brightness_eq_nonsat))

% histograms 
bincol = [0.8,0.8,0.8];
legendstr = {'source','target','matched source'};
A = {source_brightness(~source_indx_saturated),target_brightness(~target_indx_saturated),...
    source_brightness_eq(~source_indx_saturated)};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');
