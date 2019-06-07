function result = linear_solver(A,b,solver_kind)
% LINEAR_SOLVER

fprintf('-------------- Solve --------------\n')
fprintf('Size of A = %d * %d .\n',size(A,1),size(A,2));
fprintf('Size of b = %d * %d .\n',size(b,1),size(b,2));
result = zeros(length(b),1);

if (solver_kind == 1)
    tic
    fprintf('-----------------------------------\n')
    fprintf('Solver == LU \n');
    result = A\b;
    fprintf('Relative residual = %6.4e . \n',norm(A*result-b)/norm(b));
    toc
    fprintf('-----------------------------------\n')
elseif (solver_kind == 2)
    result = jr_solver(A,b);
elseif (solver_kind == 3)
    result = gs_solver(A,b);
elseif (solver_kind == 4)
    result = sor_solver(A,b);
elseif (solver_kind == 5)
    result = bicg_solver(A,b);
elseif (solver_kind == 6)
    result = cgs_solver(A,b);
end

% https://ww2.mathworks.cn/help/matlab/math/sparse-matrix-operations.html

end % function linear_solver

