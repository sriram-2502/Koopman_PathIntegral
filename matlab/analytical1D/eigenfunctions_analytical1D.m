% Path integral Eigenfunctions
% Anaytical Example
clear; clc; close all;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% set up the system for positive lambda2
% sys 1
f = @(t, x) [-x(1,:)+x(1,:).^3];

% sys 2
% f = @(t, x) [lambda1*x(1,:); lambda2*(x(2,:)+x(1,:).^2)];
% beta=-lambda2/(2*lambda1-lambda2);

% get quiver
% [X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
% u = (lambda1.*X);
% v = (lambda2.*Y-X.^2);

% true eig functions
x = sym('x', 'real');
phi = @(x) x(1,:)/sqrt(1-x(1,:).^2);

%% linearize the system
A = double(subs(jacobian(f(0,x),x),x,[0]));
fn = @(x) f(0,x)-A*x;
[V,D,W] = eig(A); 
l = D(1,1);
w = W(:,1); 
g = @(x) w'*fn(x);

%% Set up path Integral
Dom = [-0.9 0.9]; ds = 0.005;
% Dom = [-2 2]; ds = 0.05;
grid = Dom(1):ds:Dom(2); %define grid where eigenfunction is well defined
x_0 = [grid(:)]; phi_est=[];

%options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
options = odeset('RelTol',1e-9,'AbsTol',1e-300);

t_span = [-200 0];
parfor i = 1:length(x_0)
    % usde ode45 to test eig funs
    [t,x] = ode45(@(t,x)f(t,x),t_span,x_0(i,:),options);

    % for negative time
    phi_est = [phi_est, w'*x_0(i,:)' - trapz(t,exp(-l*t).*g(x')')];
end

%% get true eig funs
Phi = zeros(size(grid));
for i = 1:length(grid)
    Phi(i) = phi(grid(i)) ;
end
phi_est = reshape((phi_est),size(grid));

%% Eigenfunctions
figure(1)
subplot(1,2,1)
plot(grid,Phi);
title('true $\phi(x) = \frac{-x}{\sqrt{1-x^2}}$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$\phi(x)$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(1,2,2)
plot(grid,phi_est); hold on;
title('estimated $\phi(x)$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$\phi(x)$','interpreter','latex');
set(gca,'fontsize',20)
axis square

% subplot(1,3,3)
% Error = abs((phi2_est - Phi22)./(Phi22));
% error = norm(phi2_est - Phi22);
% surf(q1,q2, Error)
% title('Error = $||\phi_2(x)-\hat \phi_2(x)||$','Interpreter','latex')
% xlabel('$x_1$','interpreter','latex');
% ylabel('$x_2$','interpreter','latex');
% set(gca,'fontsize',20)
% axis square


%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end

