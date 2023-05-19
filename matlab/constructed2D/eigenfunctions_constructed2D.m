% Path integral Eigenfunctions
% Anaytical Example
clear; clc; close all;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% set up the system for positive lambda2
Dom = [-0.9 0.9]; ds = 0.05;

%sys 2
lambda1 = -1;
lambda2 = -2;
f = @(t, x) [lambda1*(x(1,:)-x(1,:).^3); lambda2*(x(2,:)-x(2,:).^3)];

% get quiver
[X,Y] = meshgrid(Dom(1):0.05:Dom(2),Dom(1):0.05:Dom(2));
u = X-X.^2;
v = -Y+Y.^3;

% true eig functions
phi1 = @(x) x(1,:)./sqrt(1-x(1,:)^2);
phi2 = @(x) x(2,:)./sqrt(1-x(2,:)^2);
x = sym('x',[2;1], 'real');
Phi = @(x)[phi1(x);phi2(x)];

%% linearize the system
A = double(subs(jacobian(f(0,x),x),x,[0;0]));
fn = @(x) f(0,x)-A*x;
[V,D,W] = eig(A); 
l1 = D(1,1); l2 = D(2,2);
w1 = W(:,1); w2 = W(:,2);

g1 = @(x) w1'*fn(x);
g2 = @(x) w2'*fn(x);

%% Set up path Integral

grid = Dom(1):ds:Dom(2); %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);

x_0 = [q1(:),q2(:)]; phi1_est=[];phi2_est=[];
%options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
options = odeset('RelTol',1e-9,'AbsTol',1e-300);

tspan = [0 20]; % for negative time
parfor i = 1:length(x_0)
    [t,x] = ode45(@(t,x)f(t,x),tspan,x_0(i,:),options);
    phi1_est = [phi1_est, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x')')];
end

tspan = [0 20]; % for positive time
parfor i = 1:length(x_0)
    [t,x] = ode45(@(t,x)f(t,x),tspan,x_0(i,:),options);
    phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
end

%% get true eig funs
Phi11 = zeros(size(q1));
Phi22 = zeros(size(q2));
for i = 1:length(grid)
    for j = 1:length(grid)
        Phi11(i,j) = phi1([grid(i);grid(j)]) ;
        Phi22(i,j) = phi2([grid(i);grid(j)]) ;
    end
end
phi1_est = reshape((phi1_est),size(q2));
phi2_est = reshape((phi2_est),size(q2));

%% plots for
figure(2)
subplot(2,4,1)
p1 = pcolor(q1,q2,Phi11'); hold on;
set(p1,'Edgecolor','none')
colormap jet
%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlim(Dom)
ylim(Dom)
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(2,4,2)
p2 = pcolor(q1,q2,phi1_est); hold on;
set(p2,'Edgecolor','none')
colormap jet
%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlim(Dom)
ylim(Dom)
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar


subplot(2,4,3)
p1 = pcolor(q1,q2,Phi22'); hold on;
set(p1,'Edgecolor','none')
colormap jet
%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlim(Dom)
ylim(Dom)
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(2,4,4)
p2 = pcolor(q1,q2,phi2_est); hold on;
set(p2,'Edgecolor','none')
colormap jet
%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlim(Dom)
ylim(Dom)
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar


%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-3);
isterminal=1;
direction=0;
end

