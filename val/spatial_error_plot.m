clear; close all; clc;

M = 0.2;
H = 1;
Re = 100;
U = M / sqrt(3);
error = [];

x = ( linspace(0,1,10) );
v = dlmread('../results/v_G10_M200_T100_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,15) );
v = dlmread('../results/v_G15_M200_T100_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,20) );
v = dlmread('../results/v_G20_M200_T100_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,25) );
v = dlmread('../results/v_G25_M200_T20_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,50) );
v = dlmread('../results/v_G50_M200_T20_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,75) );
v = dlmread('../results/v_G75_M200_T20_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,100) );
v = dlmread('../results/v_G100_M200_T20_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,128) );
v = dlmread('../results/v_G128_M200_T20_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,150) );
v = dlmread('../results/v_G150_M200_T40_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

x = ( linspace(0,1,200) );
v = dlmread('../results/v_G200_M200_T30_RM1_VCM0_VCF0_Re100.dat')/U;
error = [error, reference_error(x,v,Re)];

ginv = 1./[10,15,20,25,50,75,100,128,150,200];
a1 = log(ginv);
a2 = log(error);
figure;
plot(a1, a2);
title('Error Convergence Plot');
xlabel('Logarithm of Normalized Grid Spacing');
ylabel('Logarithm of L2 Norm Error');

slope = ( a2(4) - a2(1) ) / ( a1(4) - a1(1) );
slope2 = ( a2(end) - a2(1) ) / ( a1(end) - a1(1) );