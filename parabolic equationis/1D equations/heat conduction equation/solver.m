function result = solver(model)
% SOLVER

% check r
fprintf('r = %4.2e\t',model.r);
if (model.r > 0.5)
    fprintf('More than 1/2.\n');
else
    fprintf('Less than 1/2.\n');
end

% solve
t = 0.0;
u = model.u_init;
while (t < model.t_total)
    u = update_solution(u,model,t,1);
    t = t+model.dt;
end
result = u;

end % solver

