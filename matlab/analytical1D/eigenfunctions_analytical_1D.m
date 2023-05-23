% Path integral Eigenfunctions
% Anaytical Example
clear; clc; %close all;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])
use_euler = false; % flag to use euler integration or not

%% set up the system for positive lambda2
% sys 1 stable
% lambda1 = -1;
% f = @(t, x) [-x(1,:)+x(1,:).^3];
% % true eig functions
% x = sym('x', 'real');
% phi = @(x) -lambda1*x(1,:)/sqrt(1-x(1,:).^2);

%sys 2 unstable works for alpha = 0.1 but not accurate
alpha=1;
f = @(t, x) alpha*[(x(1,:)-(x(1,:)).^3)];
% true eig functions
x = sym('x', 'real');
phi = @(x) x(1,:)/sqrt(1-x(1,:).^2);

%sys 3 
% alpha=1;
% f = @(t, x) alpha*[cos(x(1,:))];
% % true eig functions
% % possible eigenfunction
% % x = sym('x', 'real');
% % phi = @(x) -(1/sin(x(1,:))+1/tan((x(1,:))));
% x = sym('x', 'real');
% phi = @(x) x(1,:)/sqrt(1-x(1,:).^2);

%sys 4 works with cosine with alpha = 1
% alpha=1;
% f = @(t, x) alpha*[x-x.*cos(x(1,:))];
% % true eig functions
% x = sym('x', 'real');
% phi = @(x) x(1,:)/sqrt(1-x(1,:).^2);

%sys 5 works with sine with alpha = 0.1
% alpha=0.1;
% f = @(t, x) alpha*[x-x.*sin(x(1,:))];
% % true eig functions
% x = sym('x', 'real');
% phi = @(x) x(1,:)/sqrt(1-x(1,:).^2);

%% linearize the system
A = double(subs(jacobian(f(0,x),x),x,[0]));
fn = @(x) f(0,x)-A*x;
[V,D,W] = eig(A); 
l = D(1,1);
w = W(:,1); 
g = @(x) w'*fn(x);

%% Set up path Integral
dim=1;
Dom = [-0.9 0.9]; ds = 0.001;
% Dom = [-4 4]; ds = 0.05;
grid_pts = Dom(1):ds:Dom(2); %define grid where eigenfunction is well defined
x_0 = [grid_pts(:)]; phi_est_nonlin=[];phi_est_lin=[]; phi_est_direct=[];

%options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
options = odeset('RelTol',1e-9,'AbsTol',1e-300);

% t_span = [-1 0];
t_span = [0 -10];

if(use_euler)
    x_euler=[];t_euler=[];
    N = 1000;dt = 0.01;
end

for i = 1:length(x_0)
    % using euler
    if(use_euler)
        x_euler = x_0(i); t_euler=0;
        x = x_euler; t = t_euler;
        for n=1:N
            t_euler = t_euler + dt;
            x_euler = x_euler + dt*f(t_euler,x_euler); 
            x = [x;x_euler];
            t = [t;t_euler];
        end
    else
        % usde ode45 to test eig funs
        [t,x] = ode45(@(t,x)f(t,x),t_span,x_0(i,:),options);
    end

    %for forward time
    phi_est_direct = [phi_est_direct, w'*x_0(i,:)+trapz(t,exp(-l*t).*g(x')',dim)];
    phi_est_nonlin = [phi_est_nonlin, trapz(t,exp(-l*t).*g(x')',dim)];
    phi_est_lin = [phi_est_lin, w'*x_0(i,:)];

    % for negative time
    % phi_est = [phi_est, w'*x_0(i,:)' - trapz(t,exp(-l*t).*g(x')',dim)];
    % phi_est = [phi_est, trapz(t,exp(-l*t).*g(x')',dim)];
end
phi_est = phi_est_lin + phi_est_nonlin;

%% get true eig funs
Phi = zeros(size(grid_pts));
for i = 1:length(grid_pts)
    Phi(i) = phi(grid_pts(i)) ;
end

%% Eigenfunctions
figure(1)
subplot(2,4,2)
plot(grid_pts,Phi);
% title('true $\phi(x) = \frac{x}{\sqrt{1-x^2}}$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$\phi(x)$','interpreter','latex');
axes = gca;
set(gca,'fontsize',20)
grid on
% xlim([-0.9,0.9])
ylim([-2,2])
box on
axes.LineWidth=2;
axis square

subplot(2,4,2)
hold on;
plot(grid_pts,phi_est,'--'); hold on;
% plot(grid_pts,phi_est_lin,'--'); hold on;
% title('estimated $\phi(x)$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$\phi(x)$','interpreter','latex');
axes = gca;
set(gca,'fontsize',20)
grid on
% xlim([-0.9,0.9])
ylim([-2,2])
box on
axes.LineWidth=2;
axis square

% subplot(2,2,3)
% plot(grid_pts,phi_est_lin);
% % title('true $\phi(x) = \frac{x}{\sqrt{1-x^2}}$','Interpreter','latex')
% xlabel('$x$','interpreter','latex');
% ylabel('$\phi(x)$','interpreter','latex');
% axes = gca;
% set(gca,'fontsize',20)
% grid on
% % xlim([-0.9,0.9])
% ylim([-2,2])
% box on
% axes.LineWidth=2;
% axis square
% 
% subplot(2,2,4)
% plot(grid_pts,phi_est_nonlin); hold on;
% % title('estimated $\phi(x)$','Interpreter','latex')
% xlabel('$x$','interpreter','latex');
% ylabel('$\phi(x)$','interpreter','latex');
% axes = gca;
% set(gca,'fontsize',20)
% grid on
% % xlim([-0.9,0.9])
% ylim([-2,2])
% box on
% axes.LineWidth=2;
% axis square

%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end
