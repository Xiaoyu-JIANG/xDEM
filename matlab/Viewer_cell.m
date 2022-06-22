clear
clc
close all
fclose all

boundaryHalfLentgh = 10;
boundaryGap = 1.1;


r = importdata('output_radius.dat');
kr = importdata('output_aspectRatio.dat');
xet = importdata('output_position.dat');
zet = importdata('output_rotation.dat');
bdt = importdata('output_boundary.dat');


n = 50;
tt = 0:2*pi/n:2*pi;
xx = cos(tt);
yy = sin(tt);


xlen = bdt(2)-bdt(1);
ylen = bdt(4)-bdt(3);
x = xet(1:2:end);
y = xet(2:2:end);
z = zet;
X = nan(1,length(r)*(n+2));
Y = nan(1,length(r)*(n+2));

fig = createMyDefaultFigure('',[15,10]);
clf;
ax = gca;
set(gca, 'LooseInset', get(gca, 'TightInset'));
% ax.XLabel.String = '$x$';
% ax.YLabel.String = '$y$';
axis equal;
box on;
hold on;
xlim([-13, 13])
ylim([-8, 8])
% xlim([-6.4, 6.4])
% ylim([-4, 4])


for j = 1:length(r)
    xx1 = xx * r(j);
    yy1 = yy * r(j) / kr(j);
    xx2 = cos(z(j))*xx1 - sin(z(j))*yy1;
    yy2 = sin(z(j))*xx1 + cos(z(j))*yy1;
    xx3 = xx2 + x(j);
    yy3 = yy2 + y(j);
    X((j-1)*(n+2)+1:(j)*(n+2)-1) = xx3;
    Y((j-1)*(n+2)+1:(j)*(n+2)-1) = yy3;
end
hr = rectangle('Position',[bdt(1),bdt(3),bdt(2)-bdt(1),bdt(4)-bdt(3)]);
hl = line(X,Y,'Color','k','LineStyle','-','LineWidth',1);


cellFile = fopen('quadronAnalysis_cell.dat');
numCell = fscanf(cellFile, "%d", 1);
clearvars cell
cell(numCell) = struct();
for ii = 1:numCell
    cell(ii).order = fscanf(cellFile, "%d", 1);
    cell(ii).area = fscanf(cellFile, "%e", 1);
    fscanf(cellFile, "%e", 4);
    cell(ii).cx = fscanf(cellFile, "%e", 1);
    cell(ii).x = fscanf(cellFile, "%e", cell(ii).order);
    fscanf(cellFile, "%e", cell(ii).order);
    cell(ii).cy = fscanf(cellFile, "%e", 1);
    cell(ii).y = fscanf(cellFile, "%e", cell(ii).order);
    fscanf(cellFile, "%e", cell(ii).order);
end
% cell([cell.area]<)
maxCellOrder = max([cell.order]);
c = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    255,255,51
    166,86,40
    247,129,191
    153,153,153] / 255;

for ii = 3:maxCellOrder
    index = find([cell.order] == ii);
    cell_x = [];
    cell_y = [];
    for j = 1:length(index)
        cell_x = [cell_x; cell(index(j)).x; nan];
        cell_y = [cell_y; cell(index(j)).y; nan];
    end
    if isempty(cell_x); continue; end
    cell_ps = polyshape(cell_x,cell_y,'Simplify',false);
    if ii <= 11
        plot(cell_ps,...
            'FaceColor', c(ii-2,:),...
            'FaceAlpha', 0.6,...
            'EdgeColor', 0.5*ones(1,3));
    else
        plot(cell_ps,...
            'FaceColor', c(end,:),...
            'FaceAlpha', 0.6,...
            'EdgeColor', 0.5*ones(1,3));
    end
end

