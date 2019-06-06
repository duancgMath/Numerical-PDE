function u_init = init_cond_riemann_shock(x)
% INIT_COND_RIEMANN
u_init = zeros(1,length(x));
u_init(x<0.0) = 1.0;
end % function init_cond_riemann

