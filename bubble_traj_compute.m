%clear all
%close all
%clc



%% using ode solver to solve simultaneous equations (use this to understand bubble trajectory behavior)


% Parameter declaration... hopefully they are the same as previous,
% apologies for messy program
param(1) = 6.5e-3; %r_v
%param(2) = 1;%u_theta_max
param(2) = 0.1; % lambda_inf
param(3) = 100e-6;% initial bubble radius
param(4) = 1e-6;%nu_l
param(5) = 1.2;%rho_bub
param(6) = 1000;%rho_l

t = [];
u_vort_bubpos = [];
u_bubpos = [];

% Initial Conditions
r_post0=15e-3;
post0 = [r_post0 0]; % initial position
u_bub_t0= [0 0];
r_bub_t0 = [param(3) 0];
p_g_0 = 30000;

param(7) = p_g_0;
param(8) = r_bub_t0(1);

r_v = param(1);
lambda_inf = param(2);

%%

[x_grid,y_grid] = meshgrid(linspace(-15e-3,15e-3,25),linspace(-15e-3,15e-3,25));
r_grid = (x_grid.^2 + y_grid.^2).^0.5;
theta = atan2(y_grid,x_grid);

for i=1:size(r_grid,1)
    for j=1:size(r_grid,2)
        [u_theta(i,j),p(i,j),gradp(i,j)] = mod_lamb_oseen(r_grid(i,j),r_v,lambda_inf);
    end
end

u = -u_theta.*sin(theta);
v = u_theta.*cos(theta);

%%
if 0
    % Solve ODE considering Rayleigh Plesset
    [t,y] = ode45(@(t,y) odefun(t,y,param), [0 20*(1/5000)], [post0 u_bub_t0 r_bub_t0]');
end
% Solve ODE not considering Rayleigh Plesset
[t,y] = ode15s(@(t,y) odefun(t,y,param), [0 1000*(1/5000)], [post0 u_bub_t0]');



% Calculate theta position
theta_bubpos = atan2(y(:,2),y(:,1));
% Calculate tangential velocity
u_theta_bubpos = y(:,4).*cos(theta_bubpos) - y(:,3).*sin(theta_bubpos);


%% Plotting

figure(1)
subplot(1,2,1)
hold all
quiver(x_grid./param(1),y_grid./param(1),10.*u,10.*v,'k')
plot(y(:,1)./param(1),y(:,2)./param(1),'linewidth',2)
xlim([-2.5 2.5])
ylim([-2.5 2.5])
xlabel('$x/r_v$','interpreter','latex')
ylabel('$x/r_v$','interpreter','latex')
%axis equal
% axis tight
colormap gray
axis square
box on
set(gca,'linewidth',1,'fontsize',20)


% Estimate of the local velocity of the vortex
r_bubpos = (y(:,1).^2+y(:,2).^2).^0.5;
for i=1:length(r_bubpos)
    [u_vort_bubpos(i),~,~,~] = mod_lamb_oseen(r_bubpos(i),r_v,lambda_inf);
end
%speed_bubpos = (y(:,3).^2+y(:,4).^2).^0.5;

% Velocity of the bubble along its path
figure(1)
%subplot(1,2,2)
%hold all
%%plot(t,u_theta_bubpos-u_vort_bubpos','k.-','linewidth',1)% velocity
%plot(t,u_theta_bubpos,'r','linewidth',1.5)
%plot(t,u_vort_bubpos,'k','linewidth',1.5)
%legend('Bubble Velocity','Lamb Oseen Bubble Position')
%%axis square
%box on
%set(gca,'linewidth',1,'fontsize',20)

% Velocity difference of the bubble along its path against background flow
figure(1)
subplot(1,2,2)
hold all
plot(t,u_theta_bubpos'-u_vort_bubpos,'.-','linewidth',1.5)
xlabel('$t$ [sec]','interpreter','latex')
ylabel('$u_\theta^{b}-u_\theta^{LO}$ [m/s]','interpreter','latex')
ylim([0 0.5])
xlim([0 0.15])
axis square
box on
set(gca,'linewidth',1,'fontsize',20)




