clear
clc
close all
fclose all;

boundaryHalfLentgh = 32; 
boundaryGap = 1.1;

r = importdata('output_radius.dat');
kr = importdata('output_aspectRatio.dat');
xe = importdata('output_position.dat');
ze = importdata('output_rotation.dat');
bd = importdata('output_boundary.dat');

fig = figure('units','centimeters','position',[2 2 23 23]);

n = 200;
tt = 0:2*pi/n:2*pi;
xx = cos(tt);
yy = sin(tt);

start = 0;
step = 1;

for i = 1:step:length(xe(:,1))
    
    bdt = bd(i,:);
    xlen = bdt(2)-bdt(1);
    ylen = bdt(4)-bdt(3);
    xet = xe(i,:);
    zet = ze(i,:);
    x = xet(1:2:end);
    y = xet(2:2:end);
    z = zet;
    X = nan(1,length(r)*(n+2));
    Y = nan(1,length(r)*(n+2));
    
    clf;
    ax = gca;
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    ax.XLabel.String = [];
    ax.YLabel.String = [];
    axis equal;
    box on;
    hold on;
    xlim([-boundaryHalfLentgh, boundaryHalfLentgh])
    ylim([-boundaryHalfLentgh, boundaryHalfLentgh])
%         xlim([bdt(1)*boundaryGap, bdt(2)*boundaryGap])
%         ylim([bdt(3)*boundaryGap, bdt(4)*boundaryGap])
    
    title(['Step = ',num2str(i)]);
    
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
    hold on
    hr = rectangle('Position',[bdt(1),bdt(3),bdt(2)-bdt(1),bdt(4)-bdt(3)]);
    hl = line(X,Y,'Color','k','LineStyle','-');
    
    drawnow
%     pause
end
