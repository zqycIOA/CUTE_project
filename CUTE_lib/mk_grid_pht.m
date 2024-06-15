function[positions , amp] = mk_grid_pht(dx,dz,xnum,znum,zstart)
    %% Define a small phantom with scatterers
    
    %为保证信号处理时间余量打的边界反射点
    x1=-20/1000;
    x2=20/1000;
    z1=5/1000;
    z2=50/1000;
    
    %用于仿真成像的主反射点
    boundpos=mk_bound(x1,x2,z1,z2);
    xvec = (1:xnum) * dx; xvec = xvec - mean(xvec);
    zvec = (0:znum-1) * dz; zvec = zvec + zstart;
    [xs,zs]=meshgrid(xvec,zvec);
    xs=xs(:);
    zs=zs(:);
    % xs=(-5:2.5:5)/1000;
    % zs=(20:5:40)/1000;
    pointnum=length(xs);
    positions=[[xs,zeros(pointnum,1),zs];boundpos];
    amp=[ones(pointnum,1);ones(8,1)*0.00000001];
    % positions=[[0,0,40/1000];boundpos];
    % amp=[1;ones(8,1)*0.01];
    
    % save SinglePointPha positions amp;
end