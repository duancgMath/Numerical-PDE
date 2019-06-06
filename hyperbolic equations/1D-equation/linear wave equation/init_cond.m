function u_init = init_cond(x)
% INIT_COND initial conditions for the linear wave test case

n = length(x);
u_init = zeros(1,n);

a = 0.5;
z = -0.7;
delta = 0.005;
alpha = 10.0;
beta = log(2.0)/(36*delta^2);

for i = 1:n
    if (x(i)>=-0.8 && x(i)<=-0.6)
        u_init(i) = (1.0/6.0)*(G(x(i),beta,z-delta)+ ... 
            G(x(i),beta,z+delta)+4.0*G(x(i),beta,z));
    elseif (x(i)>=-0.4 && x(i)<=-0.2)
        u_init(i) = 1.0;
    elseif (x(i)>=0.0 && x(i)<=0.2)
        u_init(i) = 1.0-abs(10.0*(x(i)-0.1));
    elseif (x(i)>=0.4 && x(i)<=0.6)
        u_init(i) = (1.0/6.0)*(F(x(i),alpha,a-delta)+ ... 
            F(x(i),alpha,a+delta)+4.0*F(x(i),alpha,a));
    end
end
        

end % function init_cond 


% ----------------------------------------
function result = G(x,beta,z)
result = exp(-beta*(x-z)^2);
end % function G

function result = F(x,alpha,a)
result = sqrt(max(0.0,1.0-alpha^2*(x-a)^2));
end % function F





