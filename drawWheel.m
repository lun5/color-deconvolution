% Adapted online at
% http://www.mathworks.com/matlabcentral/answers/105659-how-to-plot-histogram-of-hue-value-in-polar
% Luong Nguyen 04/28/2015

function [h, binValues, binEdges] = drawWheel(alpha, numberOfBins, col)

linhist = histogram(alpha, 'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'NumBins',numberOfBins);
set(gcf,'visible','off')
binValues = linhist.Values;
binEdges = linhist.BinEdges;

% if (any(r > 1) || any(r < 0))
% 	error('R must be a vector of values between 0 and 1')
% end

% if size(col,1) == 1
%     col = repmat(col,[numel(binValues),1]);
% end
% 
% if numel(binValues) ~= size(col,1)
% 	error('Length of r and cmap must be the same')
% end

n = numel(binValues);
innerRadius =  1;%80;
%outerRadius = 100;

%angles = linspace(-pi,pi,n+1);
newR = innerRadius*(1+binValues);
% Draw the hue in the annulus.
h = figure;
for k = 1:n
	%newR(k)
	%drawSpoke(innerRadius, outerRadius, angles(k), angles(k+1), cmap(k,:));
	drawSpoke(newR(k)    , innerRadius, binEdges(k), binEdges(k+1), col);
end

% Draw inner black ring.
circle_angles  = linspace(-pi,pi,n+1);
line(cos(circle_angles)*innerRadius, sin(circle_angles)*innerRadius, 'LineWidth', 3, 'Color', 'k');
axis equal;

%=============================================================
function h = drawSpoke(ri,ro,thetaStart,thetaEnd,c)
xInnerLeft  = cos(thetaStart) * ri;
xInnerRight = cos(thetaEnd)   * ri;
xOuterLeft  = cos(thetaStart) * ro;
xOuterRight = cos(thetaEnd)   * ro;

yInnerLeft  = sin(thetaStart) * ri;
yInnerRight = sin(thetaEnd)   * ri;
yOuterLeft  = sin(thetaStart) * ro;
yOuterRight = sin(thetaEnd)   * ro;

X = [xInnerLeft, xInnerRight, xOuterRight xOuterLeft];
Y = [yInnerLeft, yInnerRight, yOuterRight yOuterLeft];

h = patch(X,Y,c);
set(h,'edgeColor', 'none');