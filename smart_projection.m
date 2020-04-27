%% Analytical null-space projection
function projection = smart_projection(d,n)
    m = -1/sqrt(n);
    x = (-1+1/sqrt(n))/(n-1);
    sum_d = sum(d);
    dAAT = (n-1)*x^2*repmat(sum_d,n,1);
    dABT = repmat(sum_d,n,1);
    dABT(1) = (n-1)*(m-x)*sum_d;
    dABT = x*dABT;
    sumd_nm1d1 = d(1)*(n-1)*(m-x)+sum_d-d(1);
    dBAT =  x*repmat(sumd_nm1d1,n,1);
    mxd1 = (m-x)*d(1);
    dBBT = repmat(mxd1,n,1);
    dBBT = dBBT+d;
    dBBT(1) = (m-x)*sum_d+((n-1)*(m-x)^2-(m-x))*d(1);
    projection = dAAT+dABT+dBAT+dBBT;
end