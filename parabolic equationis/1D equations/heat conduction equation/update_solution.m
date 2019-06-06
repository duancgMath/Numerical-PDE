function result = update_solution(u,model,t,update_kind)
% UPDATE_SOLUTION

if (update_kind == 1) % forward
    u_old = u;
    u(2:end-1) = u_old(2:end-1)+model.r* ... 
        (u_old(3:end)-2.0*u_old(2:end-1)+u_old(1:end-2));
    u(1) = model.bcl(t);
    u(end) = model.bcr(t);
elseif (update_kind == 2) % backward
    
elseif (update_kind == 3) % Crank-Nickson
    
end

result = u;

end % function update_solution

