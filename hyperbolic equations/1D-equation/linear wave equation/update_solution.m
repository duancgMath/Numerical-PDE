function u = update_solution(model,u,t,num_flux)
% UPDATE_SOLUTION
u(2:end-1) = u(2:end-1)-(model.dt/model.dx)* ...
    (num_flux(model,u(2:end-1),u(3:end))-...
    num_flux(model,u(1:end-2),u(2:end-1)));
u(1) = model.u_bcl(t);
u(end) = model.u_bcr(t);
end

