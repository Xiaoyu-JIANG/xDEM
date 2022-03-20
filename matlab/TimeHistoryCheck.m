clear
clc
close all

fig = createMyDefaultFigure('Energy and stress check',[15,8]);
tiledlayout(1,2)

nexttile
hold on
box on
E = importdata('output_kineticEnergy.dat');
plot(E)
set(gca,'yscale','log');

nexttile
hold on
box on
S = importdata('output_stress.dat');
plot(S(:,[1,4]))
% ylim([99.9,100.1]);