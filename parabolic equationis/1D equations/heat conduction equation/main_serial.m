% =====================================================
% 1D heat conduction equation -- serial
% =====================================================

clear, clc

% select parameters
scheme_kind = 3;
err_kind = inf;


% ----------------------------------------------------
% define the problem
dt = 1e-5; % for any dt if scheme_kind == 2 or 3
model = heat_conduct_1d([],dt);


% ----------------------------------------------------
% 
figure;
hold on, grid on
axis([-0.05 1.05 -0.05 0.05]);

n_grid_vec = [8,16,32,64,128,256,512];
dt_vec = zeros(length(n_grid_vec),1);
dx_vec = zeros(length(n_grid_vec),1);
err_vec = zeros(length(n_grid_vec),1);

for i = 1:length(n_grid_vec)
    
    model.n_grid = n_grid_vec(i);
    
    % solve
    tic
    u = solver(model,scheme_kind);
    toc
    
    % calculate error
    [err,dx,dt] = err_analysis(u,model,err_kind);
    dt_vec(i) = dt;
    dx_vec(i) = dx;
    err_vec(i) = err;
    
    % plot solution
    plot(model.xgrid,u,'.-');
end

plot(-1:0.01:1,model.exact_solution(-1:0.01:1),...
    '-','LineWidth',1.5)

% -----------------------------------------------------
% print error
print_error(dt_vec,dx_vec,err_vec,n_grid_vec);

