

clear all
close all
clc

%% declare parameters
param(1) = 2e-3; %r_v , lamb-oseed vortex viscous core, m
param(2) = 3;%u_theta_max, maximum tangential velocity of vortex
param(3) = 20e-5;%r_bub, bubble radius
param(4) = 1e-6;%nu_l, liquid viscosity
param(5) = 1.2;%rho_bub, bubble density
param(6) = 1000;%rho_l, fluid density


%% some more parameters for drawing the background flow, not used in
% trajectory evaluation
zeta = 1.2526; % constant term in Lamb oseen to make maximum u_theta at r_v
beta = 1; % beta is 1 for lamb-oseen, some expression for bosscher's vortex model... not updated in mod_lamb_oseen function yet
r_v  = param(1);
pinf = 20000;
uinf = 5.5;
rho_bub = param(5);
u_theta_max = param(2);
r_bub = param(3);
rho_l = param(6);
nu_l = param(4);


%% Use this plot a Lamb-oseen profile for pressure and velocity. 
% r = linspace(0,5e-3,100);
% 
% for i=1:length(r)
%     
%     [u(i),p(i)] = mod_lamb_oseen(r_v,r(i),u_theta_max,pinf);
%     
% end
% 
% plot(r,u)

%% calculating the velocity field for plotting purposes only
[x_grid,y_grid] = meshgrid(linspace(-15e-3,15e-3,100),linspace(-15e-3,15e-3,100));

r_grid = (x_grid.^2 + y_grid.^2).^0.5;
theta = atan2(y_grid,x_grid);

for i=1:size(r_grid,1)
    for j=1:size(r_grid,2)
        [u_theta(i,j),p(i,j),gradp(i,j)] = mod_lamb_oseen(r_grid(i,j),r_v,u_theta_max);
    end
end

p_dom = p+pinf;

grad_pres_r = rho_l.*u_theta.^2./r_grid; 


[gradx_p,grady_p] = gradient(p,x_grid(1,2)-x_grid(1,1));
% u = u_theta.*sin(theta);
% v = u_theta.*cos(theta);
% w = 5.5.*ones(size(u));
u = -u_theta.*sin(theta);
v = u_theta.*cos(theta);


%% using ode solver to solve simultaneous equations (use this to understand bubble trajectory behavior)



% Parameter declaration... hopefully they are the same as previous,
% apologies for messy program
param(1) = 1e-3; %r_v
param(2) = 1;%u_theta_max
param(3) = 0.5e-5;% initial bubble radius
param(4) = 1e-6;%nu_l
param(5) = 1.2;%rho_bub
param(6) = 1000;%rho_l


% Initial Conditions
post0 = [0 1e-3]; % initial position
u_bub_t0= [0 0];
r_bub_t0 = [param(3) 0];
p_g_0 = 2000;

param(7) = p_g_0;
param(8) = r_bub_t0(1);

% Solve ODE
[t,y] = ode15s(@(t,y) odefun(t,y,param), [0 55*(1/5000)], [post0 u_bub_t0 r_bub_t0]');

% Plot
figure(1)
hold all
quiver(x_grid./param(1),y_grid./param(1),u,v)
plot(y(:,1)./param(1),y(:,2)./param(1),'linewidth',2)
xlim([-2 2])
ylim([-2 2])
axis equal
% axis tight
colormap gray

% 
% figure(2)
% plot(t,(y(:,3).^2+y(:,4).^2).^0.5)% velocity

%% first order solution, not needed if you can do with the above.
post0 = [2e-3 0]; % initial position
dt = 0.5e-4; %seconds
time = (0:100).*dt;
u_bub_t0= [0 0];

pos = [];
u_bub = [];

figure(1)
hold all
contourf(x_grid,y_grid,p_dom)
xlim([-15e-3 15e-3])
ylim([-15e-3 15e-3])
axis equal
colormap gray


for i=1:length(time)
    
    if i==1
        u_bub = u_bub_t0;
        pos = post0;
    end
    
    r = sqrt(pos(1)^2+pos(2)^2);
    theta = atan2(pos(2),pos(1));
    
    [u_tht,p,gradp,omega_z] = mod_lamb_oseen(r,r_v,u_theta_max);
    
%     gradp = [-gradp*sin(theta) gradp*cos(theta)];
%     u_flow     = [-u_tht*sin(theta) u_tht*cos(theta)];
%     gradp = [interp2(x_grid,y_grid,gradx_p,pos(1),pos(2)) interp2(x_grid,y_grid,grady_p,pos(1),pos(2))];
%     u_flow = [interp2(x_grid,y_grid,v,pos(1),pos(2)) interp2(x_grid,y_grid,u,pos(1),pos(2))];
    

    gradp = [pos(1)*gradp/r pos(2)*gradp/r];
    u_flow = [-pos(2)*u_tht/r pos(1)*u_tht/r];
%     pause



    
    figure(1)
    hold all
    plot(pos(1),pos(2),'r*')
    quiver(pos(1),pos(2),u_bub(1),u_bub(2),0.001,'b')
    xlim([-0.5e-2 0.5e-2])
    ylim([-0.5e-2 0.5e-2])
    drawnow
    
    
    Re_bub = 2*r_bub*vecnorm(u_flow-u_bub)/nu_l;
    
    Cd = (24/Re_bub) * (1 + 0.197*Re_bub^0.63 + 2.6e-4*Re_bub^1.38);
    alpha = abs(omega_z)*r_bub/vecnorm(u_flow-u_bub);
    Cl = 4*alpha/3;%5.82*(Re_bub^-0.5)/(alpha^0.5);
    
    crossprod = cross([(u_flow-u_bub) 0],[0 0 omega_z]);
    
    u_bub = u_bub + dt.*((-3/rho_l).*(gradp)  +  (3/4).*Cd./r_bub.*(u_flow - u_bub).*vecnorm(u_flow - u_bub)   +  (3/8).*Cl.*([crossprod(1) crossprod(2)])./alpha    ) ;


    pos = pos + (u_bub.*dt);
    
    
%     pause
end
















