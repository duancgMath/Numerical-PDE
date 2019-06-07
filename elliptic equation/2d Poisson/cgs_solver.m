function result = cgs_solver(A,b)
% CGS_SOLVER Jacobi Relaxation Scheme

fprintf('-----------------------------------\n')
fprintf('Solver == Conjugate gradients squared method. \n');

% parameters
n_iter = 1e+4;
tol = 1e-6;
sol = zeros(length(b),1);

%
tic
[sol,~,relres,i,res_vec] = cgs(A,b,tol,n_iter,[],[],sol);
toc

% print info
fprintf('Max iteration steps = %d. \n',n_iter);
fprintf('Iteration stops at %d th step. \n',i);
fprintf('Relative residual = %6.4e . \n',relres);
fprintf('-----------------------------------\n')

figure
semilogy(1:i,res_vec(2:end)/norm(b),'o-','markersize',5)
xlabel('Iteration number');
ylabel('Relative residual');
grid on

% return
result = sol;

end % function 

