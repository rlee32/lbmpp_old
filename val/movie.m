clear;clc;close all;

addpath subroutines

folder = '../picset/';
middle = 'G251_M200_T100_RM3_VCM0_VCF0_Re100000';
frames = 1000;
stride = 5;
% middle = 'G501_M300_T100_RM3_VCM0_VCF0_Re1000000';
% frames = 200;
% stride = 1;

figure;
hold on;
vw = VideoWriter('test.avi');
open(vw);
uname = [folder 'u_' middle '_1.dat'];
vname = [folder 'v_' middle '_1.dat'];
u = dlmread(uname);
v = dlmread(vname);
plot_mag(u,v);
hold on;
for k = 2:stride:frames
    uname = [folder 'u_' middle '_' num2str(k) '.dat'];
    vname = [folder 'v_' middle '_' num2str(k) '.dat'];
    u = dlmread(uname);
    v = dlmread(vname);
    plot_mag(u,v);
    writeVideo (vw, getframe(gcf));
    if (mod( k, round(frames/10) == 0 ) )
        disp(['Finished frame ' num2str(k)]);
    end
end
close(vw);