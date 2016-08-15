% calculate_optical_density from raw_image, change this into 2 steps ***** 
opticalDensity = rgb2od(rgb_image); % convert image to optical density values
opticalDensity(opticalDensity <= filterOD) = 0; % threshold low stain OD
opticalDensity(isnan(opticalDensity)) = 0; % omit NAN
opticalDensity(isinf(opticalDensity)) = 0; % omit infinity