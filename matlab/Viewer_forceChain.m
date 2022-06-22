clear
clc
close all
fclose all;

boundaryHalfLentghX = 15;
boundaryHalfLentghY = 10;
boundaryGap = 1.1;

r = importdata('output_radius.dat');
kr = importdata('output_aspectRatio.dat');
xe = importdata('output_position.dat');
ze = importdata('output_rotation.dat');
bd = importdata('output_boundary.dat');

rp = importdata('quadronAnalysis_rattlerParticles.dat');



fig = createMyDefaultFigure('snapshot with force chain',[12,8]);

n = 100;
tt = 0:2*pi/n:2*pi;
xx = cos(tt);
yy = sin(tt);

start = 0;
step = 1;

for i = 1
    
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
    Xr = nan(1,length(r)*(n+2));
    Yr = nan(1,length(r)*(n+2));
    
    clf;
    ax = gca;
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    ax.XLabel.String = [];
    ax.YLabel.String = [];
    axis equal;
    box on;
    hold on;
    xlim([-boundaryHalfLentghX, boundaryHalfLentghX])
    ylim([-boundaryHalfLentghY, boundaryHalfLentghY])
    
    for j = 1:length(r)
        if abs(x(j)) > boundaryHalfLentghX * 1.1 || abs(y(j)) > boundaryHalfLentghY * 1.1
            continue
        end
        xx1 = xx * r(j);
        yy1 = yy * r(j) / kr(j);
        xx2 = cos(z(j))*xx1 - sin(z(j))*yy1;
        yy2 = sin(z(j))*xx1 + cos(z(j))*yy1;
        xx3 = xx2 + x(j);
        yy3 = yy2 + y(j);
        if rp(j) == 1
            X((j-1)*(n+2)+1:(j)*(n+2)-1) = xx3;
            Y((j-1)*(n+2)+1:(j)*(n+2)-1) = yy3;
        else 
            Xr((j-1)*(n+2)+1:(j)*(n+2)-1) = xx3;
            Yr((j-1)*(n+2)+1:(j)*(n+2)-1) = yy3;
        end
    end
    hold on
    hr = rectangle('Position',[bdt(1),bdt(3),bdt(2)-bdt(1),bdt(4)-bdt(3)]);

    particles = polyshape(X,Y,'Simplify',false);
    hParticles = plot(particles,'linewidth',1,'facecolor',[255,217,47]/255);
    hold on

    particlesR = polyshape(Xr,Yr,'Simplify',false);
    hParticlesR = plot(particlesR,'linewidth',1,'LineStyle','-','facecolor',[50,136,189]/255);
end

%%
data = importdata('output_contact.dat','\t',1);
data = data.data;
idx = abs(data(:,3)) < boundaryHalfLentghX * 1.1 & abs(data(:,4)) < boundaryHalfLentghY * 1.1;
data = data(idx,:);
[ip,jp,cx,cy,fnx,fny,ftx,fty] ...
    = deal(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),data(:,7),data(:,8));
fx = fnx + ftx;
fy = fny + fty;
f0 = sqrt(fx.^2 + fy.^2);
fmax = max(f0);
fwidth = f0 / fmax * 5;
ip = ip + 1;
jp = jp + 1;
for i = 1:size(data,1)
    line([x(ip(i)),cx(i)], [y(ip(i)),cy(i)], 'Color',[222,45,38]/255, 'linewidth',fwidth(i));
    line([cx(i),x(jp(i))], [cy(i),y(jp(i))], 'Color',[222,45,38]/255, 'linewidth',fwidth(i));
end

