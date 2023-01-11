function dydt = odefun(t,y,param)



%% Equations of motions of bubble
% Position of bubble
pos(1) = y(1); 
pos(2) = y(2);

% Velocity of bubble
u_bub(1) = y(3);
u_bub(2) = y(4);

if 0
    r_bub = y(5);
    rdot_bub = y(6);
else
    r_bub = param(3);
end



r_v = param(1);
lambda_inf=param(2);

nu_l = param(4);
rho_bub = param(5);
rho_l=param(6);
pgt0 = param(7);
r_bubt0 = param(8);


r = sqrt(pos(1)^2+pos(2)^2);
theta = atan2(pos(2),pos(1));

[u_tht,p_vort,gradp,omega_z] = mod_lamb_oseen(r,r_v,lambda_inf);

% gradp = [-gradp*sin(theta) gradp*cos(theta)];
% u_flow     = [-u_tht*sin(theta) u_tht*cos(theta)];


gradp = [pos(1)*gradp/r pos(2)*gradp/r];
u_flow = [-pos(2)*u_tht/r pos(1)*u_tht/r];



Re_bub = 2*r_bub*vecnorm(u_flow-u_bub)/nu_l;

Cd = (24/Re_bub) * (1 + 0.197*Re_bub^0.63 + 2.6e-4*Re_bub^1.38); % from paper
alpha = abs(omega_z)*r_bub/vecnorm(u_flow-u_bub); % from paper
Cl = 4*alpha/3;%5.82*(Re_bub^-0.5)/(alpha^0.5); % from paper

crossprod = cross([(u_flow-u_bub) 0],[0 0 omega_z]);


% Gather up force components in x and y direction for equation of momentum
Fd = (3/4).*Cd./r_bub.*(u_flow - u_bub).*vecnorm(u_flow - u_bub); %drag force
Fl = (3/8).*Cl.*([crossprod(1) crossprod(2)])./alpha; %lift force
Pgrad = (-3/rho_l).*(gradp); %pressure gradient induced force


%% Assemble Rayleigh Plesset equations and polytropic compression equation
if 0
    
    p_vap = 2000;
    T = 0.071;
    % % Radius of bubble
    mu_l = nu_l*rho_l;
    
    
    % % Gas pressure inside the bubble
    pg_bub = pgt0*(r_bubt0^3/r_bub^3);
    
    drdt = rdot_bub ;
    drdt2 = ( 1/(rho_l*r_bub) * (p_vap + pg_bub - p_vort - 2*T/r_bub - 4*mu_l*rdot_bub/(r_bub^2))  ) + ...
        dot((u_flow-u_bub),(u_flow-u_bub))/4*y(5) - ...
        3*y(6)^2/(2*y(5));
    
    
    
    dxdt = [ drdt ; drdt2]'; 

end




%% Consolidate vector of all equations

if 0
    dydt = [u_bub (Fd   +  Fl + Pgrad) dxdt]';
else
    dydt = [u_bub (Fd   +  Fl + Pgrad)]';% dxdt]';
end















end
