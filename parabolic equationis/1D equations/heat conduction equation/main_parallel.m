% =====================================================
% 1D heat conduction equation -- parallel
% =====================================================

clear, clc

% select parameters
scheme_kind = 2;
err_kind = inf;
    

% ----------------------------------------------------
% 
figure;
hold on, grid on
axis([-0.05 1.05 -0.05 0.05]);


n_grid_vec = [8,16,32,64,128,256,512];
dt_vec = zeros(length(n_grid_vec),1);
dx_vec = zeros(length(n_grid_vec),1);
err_vec = zeros(length(n_grid_vec),1);

parfor i = 1:length(n_grid_vec)
    
    % define the problem
    dt = 1e-6;
    model = heat_conduct_1d(n_grid_vec(i),dt);
    
    % solve
    tic
    u = solver(model,scheme_kind);
    toc
    
    % calculate error
    [err,dx,dt] = err_analysis(u,model,err_kind);
    dt_vec(i) = dt;
    dx_vec(i) = dx;
    err_vec(i) = err;
    
end % parfor


% -----------------------------------------------------
% print error
print_error(dt_vec,dx_vec,err_vec,n_grid_vec);

