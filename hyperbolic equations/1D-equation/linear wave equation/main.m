% =================================================
% Solve 1D linear wave equation: ut+a*ux=0
% =================================================

clear, clc
% ------------------------------------------------

% define the problem
model = linear_wave_problem([],@init_cond_normal);

n_grid_vec = [2048,1024,512,256,128,64,32];
n_grid_vec = n_grid_vec(end:-1:1);

figure;
hold on, grid on
axis([-1.05 1.05 -0.05 1.05])

dx_vec = zeros(length(n_grid_vec),1);
err_vec = zeros(length(n_grid_vec),1);
for i = 1:length(n_grid_vec)
    
    model.n_grid = n_grid_vec(i);
    % initial conditions
    u = model.u_init(model.xgrid);
    
    % Solve Problem
    t = 0.0;
    while (t < model.t_total)
        u = update_solution(model,u,t,@laxfriedrich_flux);
        t = t+model.dt;
    end

    % error
    [error,dx] = error_analysis(model,u,inf);
    err_vec(i) = error;
    dx_vec(i) = dx;
    
    % plot solution
    plot(model.xgrid,u,'.-');
    
end

% print error
print_error(model,dx_vec,err_vec);

