clear
clc
close all
fclose all

boundaryHalfLentgh = 15;
boundaryGap = 1.1;

r = importdata('output_radius.dat');
kr = importdata('output_aspectRatio.dat');
xe = importdata('output_position.dat');
ze = importdata('output_rotation.dat');


data = importdata('clusterLabelling_correlation.dat');
[~,idx] = max(data(:,2));
data = importdata('clusterLabelling_mclus.dat');
ncl = data(idx,:);
data = importdata('clusterLabelling_nclus.dat');
par2cls = data(idx,:)+1;

fig = createMyDefaultFigure('',[15,10]);
clf;
ax = gca;
set(gca, 'LooseInset', get(gca, 'TightInset'));
ax.XLabel.String = '$x$';
ax.YLabel.String = '$y$';
axis equal;
box on;
hold on;
xlim([-32, 32])
ylim([-20, 20])

n = 50;
tt = 0:2*pi/n:2*pi;
xx = cos(tt);
yy = sin(tt);
start = 0;
step = 1;

cmap = [166,206,227
    31,120,180
    178,223,138
    51,160,44
    251,154,153
    227,26,28
    253,191,111
    255,127,0
    202,178,214
    106,61,154
    255,255,153
    177,89,40]/255;

clearvars clusters
clusters(length(r)) = struct();
clusters(length(r)).px = [];
clusters(length(r)).py = [];
for i = 1
    xet = xe(i,:);
    zet = ze(i,:);
    x = xet(1:2:end);
    y = xet(2:2:end);
    z = zet;
    X = nan(1,length(r)*(n+2));
    Y = nan(1,length(r)*(n+2));
    
    X0 = nan(1,length(r)*(n+2));
    Y0 = nan(1,length(r)*(n+2));
    
    for j = 1:length(r)
        
        xx1 = xx * r(j);
        yy1 = yy * r(j) / kr(j);
        xx2 = cos(z(j))*xx1 - sin(z(j))*yy1;
        yy2 = sin(z(j))*xx1 + cos(z(j))*yy1;
        xx3 = xx2 + x(j);
        yy3 = yy2 + y(j);
        
        if ncl(j) == 1
            X0((j-1)*(n+2)+1:(j)*(n+2)-1) = xx3;
            Y0((j-1)*(n+2)+1:(j)*(n+2)-1) = yy3;
        else
            X((j-1)*(n+2)+1:(j)*(n+2)-1) = xx3;
            Y((j-1)*(n+2)+1:(j)*(n+2)-1) = yy3;
        end
        
        if abs(x(j)) < 100 && abs(y(j)) < 100
            clusters(par2cls(j)).px = [clusters(par2cls(j)).px, xx3];
            clusters(par2cls(j)).py = [clusters(par2cls(j)).py, yy3];
        end
    end
    hl = line(X,Y,'Color','k','LineStyle','-','linewidth',1);
    hl0 = line(X0,Y0,'Color',[.5 .5 .5],'LineStyle',':','linewidth',0.5);
    for j = 1:length(clusters)
        if isempty(clusters(j).px) || ncl(j) < 2
            continue
        end
        k = boundary(clusters(j).px',clusters(j).py',0.8);
        patch(clusters(j).px(k),clusters(j).py(k),cmap(mod(j,size(cmap,1)-1)+1,:),'FaceAlpha',.8);
    end
    drawnow
    
end
