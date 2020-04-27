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
close all
clear
clc
%%
nelx = 60;
nely = 20;
nelz = 10;
volfrac = 0.4;
penal = 3;
rmin = 2.0;
%%
maxloop = 100;    % Maximum number of iterations
tolx = 0.01;      % Terminarion criterion
displayflag = 0;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E0 = 1;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
%% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);
xPhys = x; 
loop = 0; 
change = ones(nely,nelx,nelz);
%% START ITERATION
fig = figure;
while ~(mean(abs(change),'all') <= 2e-3 || loop>=maxloop)
    loop = loop + 1;
    xPhys = imgaussfilt3(xPhys,rmin/2);
    [c,dc] = FEA_revised(KE,xPhys,Emin,E0,nelx,nely,nelz,iK,jK,freedofs,edofMat,penal,F);
    dc = sign(dc).*min(mean(abs(dc),'all')*5,abs(dc));
    dc = imgaussfilt3(dc,rmin);
    %% Modify large and small values to avoid the zero step problem
    xPhys = non_zeros_max_step_guarrentee_3D(xPhys,volfrac,Emin);
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
    reshape_dc_modified = reshape(dc_modified,[nely,nelx,nelz]);
    modification = opt_step*reshape_dc_modified;
    xPhys = xPhys+modification;
    change = xPhys(:)-x(:);
    x = xPhys;
    %% PRINT RESULTS
    fprintf('It:%3i|Obj:%8.4f|Vol:%5.3f|dc:%8.3f|dcm:%8.3f|ch:%6.3f|ms:%e|os:%e\n', ...
      loop,full(c),mean(xPhys(:)),sum(full(abs(dc)),'all'),sum(full(abs(dc_modified)),'all'),mean(full(abs(change)),'all'),max_steps,opt_step);
end
clf; display_3D(xPhys);

%% FEA
function [c, dc] = FEA_revised(KE,xPhys,Emin,E0,nelx,nely,nelz,iK,jK,freedofs,edofMat,penal,F)
    %% FE-ANALYSIS
    nele = nelx*nely*nelz;
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U = zeros(length(K),1);
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
end

%% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end

%% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
