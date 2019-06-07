function result = update_solution(u,model,t,update_kind)
% UPDATE_SOLUTION

if (update_kind == 1) % forward
    u(2:end-1) = u(2:end-1)+model.r* ... 
        (u(3:end)-2.0*u(2:end-1)+u(1:end-2));
    u(1) = model.bcl(t);
    u(end) = model.bcr(t);
elseif (update_kind == 2) % backward
    n = length(u)-2;
    A = (1+2.0*model.r)*speye(n)+ ...
        sparse(2:n,1:n-1,-model.r*ones(n-1,1),n,n)+ ...
        sparse(1:n-1,2:n,-model.r*ones(n-1,1),n,n);
    b = u(2:end-1)'; 
    b(1) = b(1)+model.r*model.bcl(t);
    b(end) = b(end)+model.r*model.bcr(t);
    % update
    u(2:end-1) = A\b;
    % bc
    u(1) = model.bcl(t);
    u(end) = model.bcr(t);
elseif (update_kind == 3) % Crank-Nickson
    n = length(u)-2;
    A = (1.0+model.r)*speye(n)+ ...
        sparse(2:n,1:n-1,-0.5*model.r*ones(n-1,1),n,n)+ ...
        sparse(1:n-1,2:n,-0.5*model.r*ones(n-1,1),n,n);
    B = sparse(1:n,2:n+1,(1.0-model.r)*ones(n,1),n,n+2)+ ...
        sparse(1:n,  1:n,  0.5*model.r*ones(n,1),n,n+2)+ ...
        sparse(1:n,3:n+2,  0.5*model.r*ones(n,1),n,n+2);
    b = B*u';
    b(1) = b(1)+0.5*model.r*model.bcl(t);
    b(end) = b(end)+0.5*model.r*model.bcr(t);
    % update
    u(2:end-1) = A\b;
    % bc
    u(1) = model.bcl(t);
    u(end) = model.bcr(t);
end

result = u;

end % function update_solution

