clear; close all; clc;

om = 1; % relaxation frequency.

% This script is to view the MRT matrices.
% This was used to reduce the matrix equations to arithmetic operations 
%   for efficient implementation.

% Constants
M = [ones(1,9);...
    -4, -ones(1,4), 2*ones(1,4);
    4, -2*ones(1,4), ones(1,4);
    0, 1, 0, -1, 0, 1, -1, -1, 1;
    0, -2, 0, 2, 0, 1, -1, -1, 1;
    0, 0, 1, 0, -1, 1, 1, -1, -1;
    0, 0, -2, 0, 2, 1, 1, -1, -1;
    0, 1, -1, 1, -1, zeros(1,4);
    zeros(1,5), 1, -1, 1 -1];
Minv = 1/36*[4*ones(9,1), M(2,:)', M(3,:)', 6*M(4,:)', 3*M(5,:)',...
    6*M(6,:)', 3*M(7,:)', 9*M(8,:)', 9*M(9,:)'];
    
% S_vec = [1, 1.4, 1.4, 1, 1.2, 1, 1.2, om, om]'; % Mohamad.
S_vec = [1, 1.2, 1, 1, 1.2, 1, 1.2, om, om]'; % 2015 Zhang et al.
MinvS = Minv*diag(S_vec);

% Different orientation of vector indices
reorder = [1,2,6,3,7,4,8,5,9];
M_ = M(:,reorder);
Minv_coeff = 36*inv(M_);
% Minv_ = Minv_coeff / 36;
% S_vec_ = S_vec(reorder);
% MinvS_ = Minv_*diag(S_vec_);
% MinvS_ = Minv_*diag(S_vec);
% MinvS_coeff = 36*(Minv_*diag(S_vec));
S_vec2 = S_vec;
S_vec2([8,9]) = 1;
mdiff = Minv*diag(S_vec2) - MinvS;
mdiff_coeff = mdiff*36;
