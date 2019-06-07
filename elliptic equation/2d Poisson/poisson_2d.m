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
            result = (obj.xr-obj.xl)/obj.nx;
        end
        
        function result = dy(obj)
            result = (obj.yu-obj.yl)/obj.ny;
        end
        
        function result = xgrid(obj)
            result = obj.xl:obj.dx:obj.xr;
        end
        
        function result = ygrid(obj)
            result = obj.yl:obj.dy:obj.yl;
        end
        
        function result = mesh_generate(obj)
            result = meshgrid(obj.xgrid,obj.ygrid);
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
        
            B = sparse(1:obj.nx-1,1:obj.nx-1,2.0*(1.0/dx^2+1.0/dy^2)*ones(obj.nx-1,1),obj.nx-1,obj.nx-1)+ ...
                sparse(2:obj.nx-1,1:obj.nx-2,-1.0/dx^2*ones(obj.nx-2,1),obj.nx-1,obj.nx-1)+...
                sparse(1:obj.nx-2,2:obj.nx-1,-1.0/dx^2*ones(obj.nx-2,1),obj.nx-1,obj.nx-1);
            BB = kron(speye(obj.ny-1),B);

            C = sparse(2:obj.ny-1,1:obj.ny-2,-1.0/dy^2*ones(obj.nx-2,1),obj.ny-1,obj.ny-1)+...
                sparse(1:obj.ny-2,2:obj.ny-1,-1.0/dy^2*ones(obj.nx-2,1),obj.ny-1,obj.ny-1);
            CC = kron(C,speye(obj.nx-1));

            result = BB+CC;
            
        end % matrix
        
        
        function result = b(obj)
            result = ones((obj.nx-1)*(obj.ny-1),1);
        end % vector

    end % methods
    
end % classdef
    
    
