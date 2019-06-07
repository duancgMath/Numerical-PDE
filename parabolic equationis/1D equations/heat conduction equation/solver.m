function result = solver(model,scheme_kind)
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
    t = t+model.dt;
    u = update_solution(u,model,t,scheme_kind);
end
result = u;

end % solver

