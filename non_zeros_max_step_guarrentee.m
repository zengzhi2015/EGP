%% Modify the xPhys to solve the split-stepping numerical problem
function new_xPhys=non_zeros_max_step_guarrentee(xPhys,volfrac,Emin)
    gap_min = 0.3; % This should not be larger than the fraction; Reduce this value if the convergence is not stable.
    if gap_min>=volfrac
        gap_min = volfrac*0.75;
    end
    new_xPhys = xPhys;
    new_xPhys(xPhys<Emin+gap_min)=Emin;
    new_xPhys(xPhys>1.0-gap_min)=1.0;
    if length(size(xPhys)) == 2
        [nelx,nely] = size(xPhys);
        vol_diff = sum(new_xPhys,'all')-volfrac*nelx*nely;
    else
        [nely,nelx,nelz] = size(xPhys);
        vol_diff = sum(new_xPhys,'all')-volfrac*nelx*nely*nelz;
    end
    changable_flag = (new_xPhys>=Emin+gap_min) & (new_xPhys<=1.0-gap_min);
    new_xPhys(changable_flag) = xPhys(changable_flag)-vol_diff/sum(changable_flag,'all');
end
