classdef heat_conduct_1d
% u_t=nu*u_xx
    properties (Constant)
        nu = 1.0/6.0;
        xl = 0.0;
        xr = 1.0;
        t_total = 0.5;
        init_cond = @(x)sin(2.0*pi*x);
        bcl = @(x)0.0;
        bcr = @(x)0.0;
    end 
    properties
        n_grid;
        dt;
    end % properties
    
    methods
        
        function result = dx(obj)
            result = (obj.xr-obj.xl)/obj.n_grid;
        end 
        
        function result = r(obj)
            result = obj.nu*obj.dt/(obj.dx)^2;
        end 
        
        function result = xgrid(obj)
            result = obj.xl:obj.dx:obj.xr;
        end
        
        function result = u_init(obj)
            result = obj.init_cond(obj.xgrid);
        end
        
        function obj = heat_conduct_1d(n_grid,dt)
            obj.n_grid = n_grid;
            obj.dt = dt;
        end 

        function obj = set.n_grid(obj,n)
            obj.n_grid = n;
        end 
        
        function result = exact_solution(obj,x)
            result = exp(-obj.nu*(2.0*pi)^2*obj.t_total).*...
                sin(2.0*pi*x);
        end

    end % methods
      
end % classdef
 
