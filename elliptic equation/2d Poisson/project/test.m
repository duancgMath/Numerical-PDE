
%% -------------------------------------------------
% Test problem 01

clear, clc

F = @(x,y)-exp(x+y);
bc_func = @(x,y)-exp(1-x-y);
bc_xl = @(y)bc_func(0,y);
bc_xr = @(y)bc_func(1,y);
bc_yl = @(x)bc_func(x,0);
bc_yu = @(x)bc_func(x,1);

model = poisson_2d([0 1 0 1],F,bc_xl,bc_xr,bc_yl,bc_yu,20,30);
           
sol1 = model.itersol('jacobi');
model.plot_sol(sol1);

sol2 = model.itersol('gs');
model.plot_sol(sol2);

sol3 = model.solution;
model.plot_sol(sol3);


%% --------------------------------------------
% Test problem 02

clear, clc

F = @(x,y)0.0*x+0.0*y;
bc_xl = @(y)1-y.^2;
bc_xr = @(y)-y.^2;
bc_yl = @(x)1-x.^2;
bc_yu = @(x)-x.^2;

model = poisson_2d([0 1 0 1],F,bc_xl,bc_xr,bc_yl,bc_yu,20,30);
      

sol1 = model.itersol('jacobi');
model.plot_sol(sol1);

sol2 = model.itersol('gs');
model.plot_sol(sol2);

sol3 = model.solution;
model.plot_sol(sol3);


%% --------------------------------------------
% Test problem 03

clear, clc

F = @(x,y)-0.5*sin(pi*x).*cos(2.0*pi*y)+0.25*cos(pi*x/50).*sin(2.0*pi*y);
bc_xl = @(y)0.0;
bc_xr = @(y)0.01;
bc_yl = @(x)0.01*sin(pi*x/2.0);
bc_yu = @(x)-0.01*sin(3.0*pi*x/2.0);

model = poisson_2d([0 1 0 1],F,bc_xl,bc_xr,bc_yl,bc_yu,20,30);
           
sol1 = model.itersol('jacobi');
model.plot_sol(sol1);

sol2 = model.itersol('gs');
model.plot_sol(sol2);

sol3 = model.solution;
model.plot_sol(sol3);

