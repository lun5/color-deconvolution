%% matching brightness by cdf
target_brightness_nonsat = target_brightness(~target_indx_saturated);
source_brightness_nonsat = source_brightness(~source_indx_saturated);

% normalize the range of brightness to between 0 and 1
target_brightness_norm = target_brightness_nonsat./range(target_brightness_nonsat);
source_brightness_norm = source_brightness_nonsat./range(source_brightness_nonsat);
% count number of elements in each bin 
binranges = 0:0.01:1;
bincounts = histc(target_brightness_norm,binranges);
% equalize source's brightness --> target's brightness in [0 1]
source_brightness_norm_eq = histeq(source_brightness_norm,bincounts);
% multiply by the range of target's brightness to complete normalization
source_brightness_eq = source_brightness;
source_brightness_eq_nonsat = source_brightness_norm_eq * range(target_brightness);
source_brightness_eq(~ source_indx_saturated) = source_brightness_eq_nonsat;

%% statistics after equalization
sprintf('Brightness Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_cdf %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    mean(source_brightness_nonsat), std(source_brightness_nonsat), min(source_brightness),...
    max(source_brightness), skewness(source_brightness_nonsat), kurtosis(source_brightness_nonsat),...
    mean(target_brightness_nonsat), std(target_brightness_nonsat), min(target_brightness),...
    max(target_brightness), skewness(target_brightness_nonsat), kurtosis(target_brightness_nonsat),...
    mean(source_brightness_eq_nonsat), std(source_brightness_eq_nonsat),...
    min(source_brightness_eq), max(source_brightness_eq),...
    skewness(source_brightness_eq_nonsat), kurtosis(source_brightness_eq_nonsat))
%% histogram
A = {source_brightness,target_brightness,source_brightness_eq};
figure; nhist(A,'legend',legendstr,'color',bincol,'xlabel','brightness values',...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');
% 
% B = {source_brightness_norm,target_brightness_norm,source_brightness_norm_eq};
% figure; nhist(B, 'legend',legendstr,'color',bincol,...
%     'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','minx',0,'maxx',1,'pdf');
