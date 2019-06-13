classdef poisson_2d
    properties
        xl
        xr
        yl
        yu
        F % -\laplace u = F
        bc_xl
        bc_xr
        bc_yl
        bc_yu
        nx
        ny
    end % properties
    
    methods 
        
        function result = dx(obj)
            result = (obj.xr-obj.xl)/(obj.nx-1);
        end
        
        function result = dy(obj)
            result = (obj.yu-obj.yl)/(obj.ny-1);
        end
        
        function result = xgrid(obj)
            result = obj.xl:obj.dx:obj.xr;
        end
        
        function result = ygrid(obj)
            result = obj.yl:obj.dy:obj.yu;
        end
        
        function result = mesh_x(obj)
            [~,result] = meshgrid(obj.ygrid,obj.xgrid);
        end
        
        function result = mesh_y(obj)
            [result,] = meshgrid(obj.ygrid,obj.xgrid);
        end
        
        function obj = poisson_2d(domain,F_func, ...
                bc_xl,bc_xr,bc_yl,bc_yu,n_xgrid,n_ygrid)
            obj.xl = domain(1);
            obj.xr = domain(2);
            obj.yl = domain(3);
            obj.yu = domain(4);
            obj.F = F_func;
            obj.bc_xl = bc_xl;
            obj.bc_xr = bc_xr;
            obj.bc_yl = bc_yl;
            obj.bc_yu = bc_yu;    
            obj.nx = n_xgrid;
            obj.ny = n_ygrid;
        end
        
        function obj = set.nx(obj,n_xgrid)
            obj.nx = n_xgrid;
        end
        
        function obj = set.ny(obj,n_ygrid)
            obj.ny = n_ygrid;
        end
        
        function result = A(obj)

            Tnx = sparse(1:obj.nx-2,1:obj.nx-2,(1.0/obj.dx^2+1.0/obj.dy^2)*ones(obj.nx-2,1),obj.nx-2,obj.nx-2)+ ...
                  sparse(2:obj.nx-2,1:obj.nx-3,-1.0/obj.dx^2*ones(obj.nx-3,1),obj.nx-2,obj.nx-2)+...
                  sparse(1:obj.nx-3,2:obj.nx-2,-1.0/obj.dx^2*ones(obj.nx-3,1),obj.nx-2,obj.nx-2);

            Tny = sparse(1:obj.ny-2,1:obj.ny-2,(1.0/obj.dx^2+1.0/obj.dy^2)*ones(obj.ny-2,1),obj.ny-2,obj.ny-2)+ ...
                  sparse(2:obj.ny-2,1:obj.ny-3,-1.0/obj.dy^2*ones(obj.ny-3,1),obj.ny-2,obj.ny-2)+...
                  sparse(1:obj.ny-3,2:obj.ny-2,-1.0/obj.dy^2*ones(obj.ny-3,1),obj.ny-2,obj.ny-2);

            result = kron(speye(obj.ny-2),Tnx)+kron(Tny,speye(obj.nx-2));
            
        end % matrix
        
        function result = b(obj)
            result = zeros((obj.nx-2)*(obj.ny-2),1);
            fmat = obj.F(obj.mesh_x,obj.mesh_y);
            for i_col = 2:obj.ny-1
                % Inhomogeneous term
                temp = fmat(2:end-1,i_col);
                % B.C. x
                ygrid = obj.ygrid;
                temp(1) = temp(1)+(1.0/obj.dx^2)*obj.bc_xl(ygrid(i_col));
                temp(end) = temp(end)+(1.0/obj.dx^2)*obj.bc_xr(ygrid(i_col));
                result((obj.nx-2)*(i_col-2)+1:(obj.nx-2)*(i_col-1)) = temp;
            end 
            % B.C. y
            xgrid = obj.xgrid;
            bc_yl_vec = obj.bc_yl(xgrid(2:end-1))';
            bc_yu_vec = obj.bc_yu(xgrid(2:end-1))';
            result(1:obj.nx-2) = result(1:obj.nx-2)+(1.0/obj.dy^2)*bc_yl_vec;
            result(end-obj.nx+3:end) = result(end-obj.nx+3:end)+(1.0/obj.dy^2)*bc_yu_vec;
        end % vector
        
        function result = xvec(obj)
            result = obj.A\obj.b;
        end 
        
        function result = vec2sol(obj)
            result = zeros(obj.ny,obj.nx);
            % xl & xr
            result(:,1) = obj.bc_xl(obj.ygrid);
            result(:,end) = obj.bc_xr(obj.ygrid);
            % yl & yu
            xgrid = obj.xgrid;
            result(1,2:end-1) = obj.bc_yl(xgrid(2:end-1));
            result(end,2:end-1) = obj.bc_yu(xgrid(2:end-1));
            % solution
            x = obj.xvec;
            temp = zeros(obj.nx-2,obj.ny-2);
            for i_col = 2:obj.ny-1
                temp(:,i_col-1) = x((obj.nx-2)*(i_col-2)+1:(obj.nx-2)*(i_col-1));
            end 
            temp = temp';
            result(2:end-1,2:end-1) = temp;
        end % function solution
        
        function result = solution(obj)
            result = obj.vec2sol;
        end 
            
        function result = itersol(obj,solver_kind)            
            sol = zeros(obj.ny,obj.nx);
            % xl & xr
            sol(:,1) = obj.bc_xl(obj.ygrid);
            sol(:,end) = obj.bc_xr(obj.ygrid);
            % yl & yu
            xgrid = obj.xgrid;
            sol(1,2:end-1) = obj.bc_yl(xgrid(2:end-1));
            sol(end,2:end-1) = obj.bc_yu(xgrid(2:end-1));
            
            fprintf('-----------------------------------\n')
            % parameters
            maxit = 2000;
            tol = 1e-6;
            sol(2:end-1,2:end-1) = zeros(obj.ny-2,obj.nx-2);
            
            if strcmp('jacobi',solver_kind)
                fprintf('Solver = Jacobi methods. \n');
                tic
                [sol,flag,iter,err] = jacobi_solver(obj,tol,maxit,sol);
                toc
            elseif strcmp('gs',solver_kind)
                fprintf('Solver = Gauss-Seidel methods. \n');
                tic
                [sol,flag,iter,err] = gs_solver(obj,tol,maxit,sol);
                toc
            elseif strcmp('sor',solver_kind)
                
            end
            
            % print info
            if (flag)
                fprintf('** NOT Convergence. **\n');
            else
                fprintf('** Convergence. **\n');
            end
            fprintf('||x1-x2|| = %6.4e. \n',err);
            fprintf('Max iteration steps = %d. \n',maxit);
            fprintf('Iteration stops at %d th step. \n',iter);
            fprintf('-----------------------------------\n')

            % return
            result = sol;
            
        end % function itersol
            
        function plot_sol(obj,sol)
            figure
            s = surf(obj.mesh_x',obj.mesh_y',sol,...
                'FaceColor','interp',...
                'FaceAlpha',0.9,...
                'FaceLighting','gouraud');
            s.EdgeColor = 'none';
            %light
            %lighting phong
            view(2);
            xlabel('x');
            ylabel('y');
        end % plot_sol

    end % methods
    
end % classdef
   

% ---------------------------------------------------------------------
% functions
% ---------------------------------------------------------------------

function [sol,flag,iter,err] = jacobi_solver(obj,tol,maxit,sol)
flag = 0;
d = 2.0*(1.0/obj.dx^2+1.0/obj.dy^2);
xgrid = obj.xgrid;
ygrid = obj.ygrid;
for iter = 1:maxit
    sol_old = sol;
    for j = 2:obj.nx-1
        for i = 2:obj.ny-1
            sol(i,j) = 1.0/d*(obj.F(ygrid(i),xgrid(j))+ ... 
                1.0/obj.dx^2*(sol_old(i+1,j)+sol_old(i-1,j))+ ...
                1.0/obj.dy^2*(sol_old(i,j+1)+sol_old(i,j-1)));
        end    
    end
    
    err = norm(sol-sol_old,inf);
    if (err < tol)
        break;
    end
end 

if (err >= tol)
    flag = 1;
end

end


function [sol,flag,iter,err] = gs_solver(obj,tol,maxit,sol)
flag = 0;
d = 2.0*(1.0/obj.dx^2+1.0/obj.dy^2);
xgrid = obj.xgrid;
ygrid = obj.ygrid;
for iter = 1:maxit
    sol_old = sol;
    for j = 2:obj.nx-1
        for i = 2:2:obj.ny-1
            sol(i,j) = 1.0/d*(obj.F(ygrid(i),xgrid(j))+ ... 
                1.0/obj.dx^2*(sol_old(i+1,j)+sol_old(i-1,j))+ ...
                1.0/obj.dy^2*(sol_old(i,j+1)+sol_old(i,j-1)));
        end    
        for i = 3:2:obj.ny-1
            sol(i,j) = 1.0/d*(obj.F(ygrid(i),xgrid(j))+ ... 
                1.0/obj.dx^2*(sol(i+1,j)+sol(i-1,j))+ ...
                1.0/obj.dy^2*(sol(i,j+1)+sol(i,j-1)));
        end 
    end
    
    err = norm(sol-sol_old,inf);
    if (err < tol)
        break;
    end
end 

if (err >= tol)
    flag = 1;
end

end
    





