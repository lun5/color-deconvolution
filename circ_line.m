% line plot wrapped around a circle
% INPUT: angles and corresponding pdf values, color of the line
% OUPUT: the pdf line wrapped around on a circle
% LUONG NGUYEN 04/28/2015
function circ_line(angles, pdfValues, col)

x = cos(angles).*(1 + pdfValues');
y = sin(angles).*(1 + pdfValues');
%plot(x,y,'-k','LineWidth',3); hold off;
%line(x,y,'Color',col,'LineWidth',3);
PlotAxisAtOrigin(x,y,col);

