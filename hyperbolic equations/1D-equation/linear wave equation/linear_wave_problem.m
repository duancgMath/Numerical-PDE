classdef linear_wave_problem
    properties (Constant)
        a = 1.0;
        CFL = 0.9;
        xl = -1.0;
        xr = 1.0;
        t_total = 0.5;
        u_bcl = @(t)0.0;
        u_bcr = @(t)0.0;
    end
    properties 
        u_init
        n_grid
    end % properties

    methods
        function result = dx(obj)
            result = (obj.xr-obj.xl)/obj.n_grid;
        end % function dx

        function result = xgrid(obj)
            result = obj.xl:obj.dx:obj.xr;
        end % function xgrid

        function result = dt(obj)
            result = obj.CFL*obj.dx/obj.a;
        end % function dt
        
        function result = u_exact(obj,x)
            result = obj.u_init(x-obj.a*obj.t_total);
        end

        function obj = linear_wave_problem(n_grid,u_init)
            obj.n_grid = n_grid;
            obj.u_init = u_init;
        end % end function linear_wave_problem

        function obj = set.n_grid(obj,nGrid)
            obj.n_grid = nGrid;
        end % function set_ngrid

    end % methods
    
end % classdef problem

