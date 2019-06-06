function [error,dx] = error_analysis(model,u,norm_kind)
% ERROR_ANALYSIS

dx = model.dx;
error = 0.0;
exact_sol = model.u_exact(model.xgrid);
error_vec = exact_sol-u;
if (norm_kind == inf)
    error = max(abs(error_vec));
elseif (norm_kind == 2)
    error = norm(error_vec)*sqrt(model.dx);
elseif (norm_kind == 1)
    error = sum(abs(error_vec))*dx;
end

end % function error_analysis

