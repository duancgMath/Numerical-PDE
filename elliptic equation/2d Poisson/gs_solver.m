function result = gs_solver(A,b)
% GS_SOLVER Gauss-Seidel Scheme

fprintf('-----------------------------------\n')
fprintf('Solver == Gauss-Seidel Scheme. \n');

% parameters
n_iter = 1e+4;
tol = 1e-6;
sol = zeros(length(b),1);

%
i = 0;
res_vec = zeros(n_iter,1);
tic
while (i <= n_iter)
    i = i+1;
    % -------------------------------------------
    
    
    % -------------------------------------------
    res = norm(A*sol-b);
    res_vec(i) = res;
    if (res<tol)
        break
    end
end
toc

res_vec = res_vec(1:i);

% print info
fprintf('Max iteration steps = %d. \n',n_iter);
fprintf('Iteration stops at %d th step. \n',i);
fprintf('Relative residual = %6.4e . \n',res/norm(b));
fprintf('-----------------------------------\n')

figure
semilogy(1:i,res_vec/norm(b),'o-','markersize',5)
xlabel('Iteration number');
ylabel('Relative residual');
grid on

% return
result = sol;

end % function 

