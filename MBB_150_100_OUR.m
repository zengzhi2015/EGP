%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Zhi Zeng and Fulei Ma,                   %
% School of Mechano-electronics Engineering,                               %
% Xidian University,                                                       %
% No.2 South Taibai Road, Xi'an, Shaanxi, China.                           %
% Please sent your comments to:                                            %
% zhizeng@mail.xidian.edu.cn, fuleima@xidian.edu.cn                        %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "An Efficient Gradient Projection Method for Structural Topology         %
%  Optimization",                                                          %
% Zhi Zeng and Fulei Ma, 2020                                              %
%                                                                          %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127. and                                                %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
%                                                                          %
% The code a can be downloaded from the web-site:                          %
% http://https://github.com/zengzhi2015/EGP-preview-version                %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program. The code may be distributed and used   %
% for educational purposes. For commercial usage, please contact the       %
% authors.                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
clear
clc
%%
nelx=150;
nely=50;
total = nely*nelx;
volfrac=0.5;
Node_Number_Force_Applied_x=1;
Node_Number_Force_Applied_y=1;
Angle_of_Force=-pi/2;
loop_limit=50;
xPhys_init = ones(nely,nelx)*volfrac;

Node_Number_Force_Applied = (Node_Number_Force_Applied_x)*(nely+1)+Node_Number_Force_Applied_y;
X_Force = cos(Angle_of_Force);
Y_Force =sin(Angle_of_Force);

penal=3.0;
rmin=2.5;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse([2*Node_Number_Force_Applied-1 2*Node_Number_Force_Applied],[1 1],[X_Force Y_Force],2*(nely+1)*(nelx+1),1);
fixeddofs = [1:2:2*(nely+1),2*(nelx+1)*(nely+1)];
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% INITIALIZE ITERATION
x = xPhys_init;
xPhys = x;
loop = 0;
change = ones(nelx,nely);
%% START ITERATION
fig = figure;
while ~(sum(abs(change),'all') <= 1e-2 || loop>=loop_limit)
    loop = loop + 1;
    xPhys = imgaussfilt(xPhys,rmin*0.5);
    [c,dc] = FEA(KE,xPhys,Emin,E0,nelx,nely,iK,jK,freedofs,edofMat,penal,F);
    dc = sign(dc).*min(mean(abs(dc),'all')*5,abs(dc));
    dc = imgaussfilt(dc,rmin);
    %% Modify large and small values to avoid the zero step problem
    xPhys = non_zeros_max_step_guarrentee(xPhys,volfrac,Emin);
    %% Check the confliction (Recursive) (Algorithm 1)
    dc_modified = dc;
    confliction_flag = (xPhys(:)==1)>10;
    cnt_add = 1;
    while true
        % dc may not conflict with some boundary. but the projected dc MAY
        % conflict with the boundary!!!
        neg_confliction_flag = ~confliction_flag;
        num_left = sum(neg_confliction_flag);
        dc_modified(neg_confliction_flag) = smart_projection(dc_modified(neg_confliction_flag),num_left);
        dc_modified(confliction_flag) = 0;
       %% Check whether new conflictions occur
        new_confliction_flag = (xPhys(:)==1 & dc_modified(:)>=0) | (xPhys(:)==Emin & dc_modified(:)<=0);
        if sum(confliction_flag==(confliction_flag|new_confliction_flag),'all')==length(confliction_flag)
            break
        else
            confliction_flag = confliction_flag|new_confliction_flag;
            dc_modified = dc;
        end
        cnt_add = cnt_add + 1;
    end
    %% Calculate the max step
    flat_xPhys = xPhys(:);
    pos_grad_flag = dc_modified>0;
    neg_grad_flag = dc_modified<0;
    if sum(pos_grad_flag,'all')==0 && sum(neg_grad_flag,'all')==0
        max_steps = 0;
        opt_step = max_steps;
        break
    else
        if sum(pos_grad_flag,'all')~=0
            pos_min_step = min((1-flat_xPhys(pos_grad_flag))./dc_modified(pos_grad_flag),[],'all');
        else
            pos_min_step = 999999999;
        end
        if sum(neg_grad_flag,'all')~=0
            neg_min_step = min((Emin-flat_xPhys(neg_grad_flag))./dc_modified(neg_grad_flag),[],'all');
        else
            neg_min_step = 999999999;
        end
            max_steps = min(pos_min_step,neg_min_step);
    end
    opt_step = max_steps;
    %% Update
    reshape_dc_modified = reshape(dc_modified,nely,nelx);
    modification = opt_step*reshape_dc_modified;
    xPhys = xPhys+modification;
    change = max(abs(xPhys(:)-x(:)));
    x = xPhys;
    %% PRINT RESULTS
    fprintf('It:%3i|Obj:%8.4f|Vol:%5.3f|dc:%8.3f|dcm:%8.3f|ch:%6.3f|ms:%e|os:%e\n', ...
      loop,full(c),mean(xPhys(:)),sum(full(abs(dc)),'all'),sum(full(abs(dc_modified)),'all'),sum(full(abs(change)),'all'),max_steps,opt_step);
    %% PLOT DENSITIES
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end

%% FEA
function [c,dc] = FEA(KE,xPhys,Emin,E0,nelx,nely,iK,jK,freedofs,edofMat,penal,F)
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U = zeros(2*(nely+1)*(nelx+1),1);
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
end
