clear
clc
close all
fclose all;

makeGIF = false;
if makeGIF
    gifName = 'animation.gif';
end

makePeriodicBound = false;

boundaryHalfLentgh = 20;
boundaryGap = 1.1;

frameRate = 25;

r = importdata('output_radius.dat');
kr = importdata('output_aspectRatio.dat');
xe = importdata('output_position.dat');
ze = importdata('output_rotation.dat');
bd = importdata('output_boundary.dat');


fig = createMyDefaultFigure('',[10 10]);

n = 50;
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

    if makePeriodicBound
        for ii = -1:1
            for jj = -1:1
                if ii == 0 && jj == 0; continue; end
                line(X + (bdt(2)-bdt(1))*ii, Y + (bdt(4)-bdt(3))*jj, ...
                    'Color',[0.6 0.6 0.6],'LineStyle','-');
            end
        end
    end
    
    drawnow
    %     pause(0.3)
    if makeGIF
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,gifName,'gif', 'Loopcount', inf, 'DelayTime', 1/frameRate);
        else
            imwrite(imind,cm,gifName,'gif', 'WriteMode', 'append', 'DelayTime', 1/frameRate);
        end
    end
end
