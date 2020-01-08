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
% "Highly Efficient Feasible Direction Method (HEFDiM) for Structural      %
%  Topology Optimization",                                                 %
% Zhi Zeng and Fulei Ma, submitting to Struct Multidisc Optim, 2020        %
%                                                                          %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127. and                                                %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site:                                            %
% https://arxiv.org/abs/2001.01896                                         %
% https://github.com/zengzhi2015/HEFDiM-preview-version                    %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program. The code may be distributed and used   %
% for educational purposes. For commercial usage, please contact the       %
% authors.                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
%%
nelx=100;
nely=100;
total = nely*nelx;
volfrac=0.2;

Node_Number_Force_Applied_x=nelx/2;
Node_Number_Force_Applied_y=1;
Node_Number_Force_Applied1 = (Node_Number_Force_Applied_x)*(nely+1)+Node_Number_Force_Applied_y;
Node_Number_Force_Applied_x=nelx/2;
Node_Number_Force_Applied_y=nely+1;
Node_Number_Force_Applied = (Node_Number_Force_Applied_x)*(nely+1)+Node_Number_Force_Applied_y;
%%
Angle_of_Force=pi/2;
Angle_of_Force1=pi/2;

loop_limit=100;
xPhys_init = ones(nely,nelx)*volfrac;
F2 = -0.1;

X_Force = cos(Angle_of_Force);
Y_Force =sin(Angle_of_Force);
X_Force1 = cos(Angle_of_Force1);
Y_Force1 =sin(Angle_of_Force1);

penal=3.0;
rmin=1.6;
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

F1 = sparse([2*Node_Number_Force_Applied1-1 2*Node_Number_Force_Applied1],[1 1],[F2*X_Force1 F2*Y_Force1],2*(nely+1)*(nelx+1),1);
Num_fix = 5;
fixeddofs1 = [2*(nely+1)-1:2*(nely+1):2*(nely+1)*Num_fix-1,2*(nely+1):2*(nely+1):2*(nely+1)*Num_fix];
fixeddofs2 = [2*(nelx-Num_fix+2)*(nely+1)-1:2*(nely+1):2*(nelx+1)*(nely+1)-1,2*(nelx-Num_fix+2)*(nely+1):2*(nely+1):2*(nelx+1)*(nely+1)];
fixeddofs=[fixeddofs1 fixeddofs2];

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
  xPhys = imgaussfilt(xPhys,rmin/2);
  [c,dc] = FEA_Sigmud(KE,xPhys,Emin,E0,nelx,nely,iK,jK,freedofs,edofMat,penal,F,F1);
  dc = imgaussfilt(dc,rmin);
  %% Modify large and small values to avoid the zero step problem
  xPhys = non_zeros_max_step_guarrentee(xPhys,volfrac,Emin);
  %% Check the confliction (Recursive)  (Algorithm 1)
  dc_modified = dc;
  confliction_flag = (xPhys(:)==1)>10;
  cnt_add = 1;
  while true
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
  fprintf('It:%3i|Obj:%8.4f|Vol:%5.3f|dc:%8.3f|dcm:%8.3f|ch:%6.3f|ms:%7.7f|os:%7.7f\n', ...
      loop,full(c),mean(xPhys(:)),sum(full(abs(dc)),'all'),sum(full(abs(dc_modified)),'all'),sum(full(abs(change)),'all'),max_steps,opt_step);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end

%%
function [c, dc] = FEA_Sigmud(KE,xPhys,Emin,E0,nelx,nely,iK,jK,freedofs,edofMat,penal,F,F1)
  %% FE-ANALYSIS
  F2=-0.1;
  kout = 0.01;
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U = zeros(2*(nely+1)*(nelx+1),1);
  V = zeros(2*(nely+1)*(nelx+1),1);
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  V(freedofs) = K(freedofs,freedofs)\F1(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce_SE = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  ce_SE2 = reshape(sum((V(edofMat)*KE).*V(edofMat),2),nely,nelx);
  ce_MSE = reshape(sum((V(edofMat)*KE).*U(edofMat),2),nely,nelx);
  Delta_11 = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce_SE));
  Delta_22 = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce_SE2))/F2;
  Delta_21 = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce_MSE))/F2;
  Delta_12 = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce_MSE));
  
  dDelta_11 = penal*(E0-Emin)*xPhys.^(penal-1).*ce_SE; 
  dDelta_22 = penal*(E0-Emin)*xPhys.^(penal-1).*ce_SE2/F2;  
  dDelta_21 = penal*(E0-Emin)*xPhys.^(penal-1).*ce_MSE/F2; 
  dDelta_12 = penal*(E0-Emin)*xPhys.^(penal-1).*ce_MSE; 
  
  MM = (F2*Delta_11+Delta_11*Delta_22*kout-Delta_21*Delta_12*kout);
  dMM = F2*dDelta_11+kout*(dDelta_11*Delta_22+Delta_11*dDelta_22-dDelta_12*Delta_21-Delta_12*dDelta_21);
  c = F2*Delta_21 / MM;
  dc = F2/MM^2*(dDelta_21*MM-Delta_21*dMM);
end
