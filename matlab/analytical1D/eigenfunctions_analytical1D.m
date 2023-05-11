% Path integral Eigenfunctions
% Anaytical Example
clear; clc; close all;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% set up the system for positive lambda2
% sys 1
f = @(t, x) [-x(1,:)+x(1,:).^3];

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
dim=1;
% Dom = [-0.9 0.9]; ds = 0.005;
Dom = [-2 2]; ds = 0.001;
grid = Dom(1):ds:Dom(2); %define grid where eigenfunction is well defined
grid(grid==1)=nan;
grid(grid==-1)=nan;
x_0 = [grid(:)]; phi_est=[];

%options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
options = odeset('RelTol',1e-9,'AbsTol',1e-300);

t_span = [-10 0];
parfor i = 1:length(x_0)
    % usde ode45 to test eig funs
    [t,x] = ode45(@(t,x)f(t,x),t_span,x_0(i,:),options);

    % for negative time
%     phi_est = [phi_est, w'*x_0(i,:)' - trapz(t,exp(-l*t).*g(x')',dim)];
    phi_est = [phi_est, trapz(t,exp(-l*t).*g(x')',dim)];

    % for flow reversed
    %phi_est = [phi_est, w'*x_0(i,:)' + trapz(t,exp(-l*t).*g(x')',dim)];
end

%% get true eig funs
Phi = zeros(size(grid));
for i = 1:length(grid)
    Phi(i) = phi(grid(i)) ;
end
% phi_est = reshape((phi_est),size(grid)-2);

%% Eigenfunctions
figure(1)
subplot(1,2,1)
plot(grid,Phi);
title('true $\phi(x) = \frac{x}{\sqrt{1-x^2}}$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$\phi(x)$','interpreter','latex');
set(gca,'fontsize',20)
xlim([-0.9,0.9])
axis square

subplot(1,2,2)
plot(grid,phi_est); hold on;
title('estimated $\phi(x)$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$\phi(x)$','interpreter','latex');
set(gca,'fontsize',20)
xlim([-0.9,0.9])
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

