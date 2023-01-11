function [u,p,gradp,omegaz] = mod_lamb_oseen(r,r_v,lambda_inf)





zeta = 1.2526;
beta = 1;
rho_l = 1000;

%lambda_inf = u_theta_max*2*pi*r_v/(1 - exp(-zeta)) ;

Const = (lambda_inf)/(2*pi);
k2=r_v^2;


u = (lambda_inf/(2*pi*r))* beta*(1 - exp(-zeta*(r^2)/k2)); 
        

IInt = (0.5 - (beta*exp(-zeta*r^2/k2)) + (0.5*beta^2*(exp(-2*zeta*r^2/k2))) + ((beta*zeta*r^2/k2)* expint(zeta* r^2/k2)) - ((beta^2*zeta*r^2/k2)*expint(2*zeta*r^2/k2))) ;
p =  (-rho_l*Const^2*IInt/r^2) ;

gradp = rho_l*u^2/r;

omegaz = (lambda_inf/(4*pi*r_v^2))*exp(-r^2/r_v^2);

end