
clear, clc

% --------------------------------------------- 
dx = 0.05;
dy = 0.02;
nx = 100;
ny = 100;

B = sparse(1:nx-1,1:nx-1,2.0*(1.0/dx^2+1.0/dy^2)*ones(nx-1,1),nx-1,nx-1)+ ...
    sparse(2:nx-1,1:nx-2,-1.0/dx^2*ones(nx-2,1),nx-1,nx-1)+...
    sparse(1:nx-2,2:nx-1,-1.0/dx^2*ones(nx-2,1),nx-1,nx-1);

BB = kron(speye(ny-1),B);

C = sparse(2:ny-1,1:ny-2,-1.0/dy^2*ones(nx-2,1),ny-1,ny-1)+...
    sparse(1:ny-2,2:ny-1,-1.0/dy^2*ones(nx-2,1),ny-1,ny-1);

CC = kron(C,speye(nx-1));

A = BB+CC;


b = ones((nx-1)*(ny-1),1);

% -------------------------------------------

for i = [1,2,5,6]
    x = linear_solver(A,b,i);
end

