clear; close all; clc;

addpath subroutines

% We judge accuracy based on the deviation in the 
%   x vs. v velocity profile, because it is harder to match than
%   the y vs. u velocity profile.


figure;
hold on;
title('Spatial Error Convergence Plot for Re = 100, M = 0.1');
xlabel('Logarithm of Normalized Grid Spacing');
ylabel('Logarithm of L2 Norm Error');

error = []; % We build this vector up.
M = 0.1;
Re = 100;
U = M / sqrt(3);

v = dlmread('../results/fields/v_G11_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G25_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G51_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G75_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G101_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G125_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G151_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G175_M100_T100_RM1_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

ginv = 1./[ 11, 25, 51, 75, 101, 125, 151, 175 ];
a1 = log10(ginv);
a2 = log10(error);
plot(a1, a2, 'x-','DisplayName','SRT, Re = 100');
coefficients = polyfit(a1, a2, 1);
disp(['Overall SRT Re = 100 order of accuracy: ' ...
    num2str(coefficients(1))]);

error = []; % We build this vector up.
M = 0.1;
Re = 100;
U = M / sqrt(3);

v = dlmread('../results/fields/v_G11_M100_T100_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G25_M100_T100_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G51_M100_T100_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G75_M100_T50_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G101_M100_T100_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G125_M100_T100_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G151_M100_T100_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G175_M100_T100_RM3_VCM0_VCF0_Re100.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

ginv = 1./[ 11, 25, 51, 75, 101, 125, 151, 175 ];
a1 = log10(ginv);
a2 = log10(error);
plot(a1, a2, 'x-','DisplayName','MRT, Re = 100');
coefficients = polyfit(a1, a2, 1);
disp(['Overall MRT Re = 100 order of accuracy: ' ...
    num2str(coefficients(1))]);
x = [a1(1), a1(end)];
y = [a2(1), a2(1)+diff(x)];
plot(x,y,':','DisplayName','Ideal 1st-Order Model','LineWidth',2);
x = [a1(1), a1(end)];
y = [a2(1), a2(1)+2*diff(x)];
plot(x,y,':','DisplayName','Ideal 2nd-Order Model','LineWidth',2);

plot(log10([1 1]/129), [-1 -3], ':', ...
    'DisplayName', 'Ghia et al. resolution');
legend('show');



figure;
hold on;
title('Spatial Error Convergence Plot for Re = 1000, M = 0.1');
xlabel('Logarithm of Normalized Grid Spacing');
ylabel('Logarithm of L2 Norm Error');


error = []; % We build this vector up.
M = 0.1;
Re = 1000;
U = M / sqrt(3);

v = dlmread('../results/fields/v_G51_M100_T100_RM1_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G75_M100_T100_RM1_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G101_M100_T100_RM1_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G125_M100_T100_RM1_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G151_M100_T100_RM1_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

ginv = 1./[ 51, 75, 101, 125, 151 ];
a1 = log10(ginv);
a2 = log10(error);
plot(a1, a2, 'x-','DisplayName','SRT');
coefficients = polyfit(a1, a2, 1);
disp(['Overall SRT Re = 1000 order of accuracy: ' ...
    num2str(coefficients(1))]);



error = []; % We build this vector up.
M = 0.1;
Re = 1000;
U = M / sqrt(3);

v = dlmread('../results/fields/v_G25_M100_T100_RM3_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G51_M100_T100_RM3_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G75_M100_T100_RM3_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G101_M100_T150_RM3_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G125_M100_T100_RM3_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

v = dlmread('../results/fields/v_G151_M100_T100_RM3_VCM0_VCF0_Re1000.dat')/U;
[x, vc] = get_centerline_v(v);
error = [error, reference_error(x,vc,Re)];

ginv = 1./[25, 51, 75, 101, 125, 151];
a1 = log10(ginv);
a2 = log10(error);
plot(a1, a2, 'x-','DisplayName','MRT');
coefficients = polyfit(a1, a2, 1);
disp(['Overall MRT Re = 1000 order of accuracy: ' ...
    num2str(coefficients(1))]);
x = [a1(1), a1(end)];
y = [a2(1), a2(1)+diff(x)];
plot(x,y,':','DisplayName','Ideal 1st-Order Model','LineWidth',2);
x = [a1(1), a1(end)];
y = [a2(1), a2(1)+2*diff(x)];
plot(x,y,':','DisplayName','Ideal 2nd-Order Model','LineWidth',2);
plot(log10([1 1]/129), [-0.8 -2], ':', ...
    'DisplayName', 'Ghia et al. grid spacing');


legend('show');
