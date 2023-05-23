% Path integral Eigenfunctions
% Anaytical Example
clear; clc; close all;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% set up the system for positive lambda2
x = sym('x',[2;1], 'real');

% sys 1 (simple analytical 2d example)
% lambda1=-4; lambda2=-3;
% f = @(t, x) [lambda1*x(1,:); (lambda2*x(2,:)-x(1,:).^2)];
% beta=1/(2*lambda1-lambda2); % for x2  
% 
% % true eig functions - sys 1 and 2
% phi1 = @(x) x(1,:);
% phi2 = @(x) x(2,:)+beta*x(1,:).^2;
% x = sym('x',[2;1], 'real');
% Phi = @(x)[phi1(x);phi2(x)];
% 
% % get quiver
% [X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
% u = (lambda1.*X);
% v = (lambda2.*Y-X.^2);

% sys 2
% https://epubs.siam.org/doi/pdf/10.1137/17M116207X
Dom = [-1 1]; ds = 0.05;
a1=-3; a2=-3;
f = @(t, x) [-2*a2*x(2,:).*(x(1,:).^2-x(2,:)-2*x(1,:).*x(2,:).^2+x(2,:).^4) +  ...
                   a1*(x(1,:)+4*x(1,:).^2.*x(2,:)-x(2,:).^2-8*x(1,:).*x(2,:).^3+4*x(2,:).^5);
           2*a1*(x(1,:)-x(2,:).^2).^2-a2*(x(1,:).^2-x(2,:)-2*x(1,:).*x(2,:).^2+x(2,:).^4)];

% % analytical Eigenfunctions
phi1 = @(x) x(1,:)-x(2,:).^2;
phi2 = @(x) -x(1,:).^2+x(2,:)+2*x(1,:).*x(2,:).^2-x(2,:).^4;

% get quiver
[X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
u = -2*a2*Y.*(X.^2-Y-2*X.*Y.^2+Y.^4)+  ...
                   a1*(X+4*X.^2.*Y-Y.^2-8*X.*Y.^3+4*Y.^5);
v = 2*a1*(X-Y.^2).^2-a2*(X.^2-Y-2*X.*Y.^2+Y.^4 );

% sys 3
% From shankar
% Dom = [-0.4 0.4]; ds = 0.05;
% lambda1 = -5;
% lambda2 = -10;
% f = @(t, x) [-lambda1*(-11*x(1,:).^3+11*x(1,:)+6*x(2,:).^3+6*x(2,:)) / (7*(3*x(1,:).^2-1));
%            -lambda2*(-x(1,:).^3+x(1,:)+5*x(2,:).^3+5*x(2,:)) / (7*(3*x(2,:).^2-1))];
% 
% % % analytical Eigenfunctions
% phi1 = @(x) 2*x(1,:).^3-2*x(1,:) -3*x(2,:).^3-3*x(2,:);
% phi2 = @(x) -x(1,:).^3+x(1,:) -2*x(2,:).^3-2*x(2,:);
% 
% % get quiver
% [X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
% u = -lambda1*(-11*X.^3+11*X+6*Y.^3+6*Y)/(7*(3*X.^2-1));
% v = -lambda2*(-X.^3+X+5*Y.^3+5*Y)/(7*(3*Y.^2-1));

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
% options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
options = odeset('RelTol',1e-9,'AbsTol',1e-300);

% for negative time
tspan = [0 10]; 
parfor i = 1:length(x_0)
    [t,x] = ode45(@(t,x)f(t,x),tspan,x_0(i,:),options);
    phi1_est = [phi1_est, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x')')];
    phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
end

% for positive time
% tspan = [0 10];
% parfor i = 1:length(x_0)
%     [t,x] = ode45(@(t,x)f(t,x),tspan,x_0(i,:),options);
%     phi1_est = [phi1_est, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x')')];
%     phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
% end

%% get true eig funs
Phi11 = zeros(size(q1));
Phi22 = zeros(size(q2));
for i = 1:length(grid)
    for j = 1:length(grid)
        Phi11(i,j) = phi1([grid(i);grid(j)]) ;
        Phi22(i,j) = phi2([grid(i);grid(j)]) ;
    end
end
phi2_est = reshape((phi2_est),size(q2));
phi1_est = reshape((phi1_est),size(q2));

%% Eigenfunctions
figure(1)
subplot(2,3,1)
surf(q1',q2',Phi11);
title('true $\phi_1(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(2,3,2)
surf(q1,q2,phi1_est); hold on;
title('estimated $\phi_1(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(2,3,3)
Error = abs((phi1_est - Phi11)./(Phi11));
error = norm(phi1_est - Phi11);
surf(q1,q2, Error)
title('Error = $||\phi_1(x)-\hat \phi_1(x)||$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(2,3,4)
surf(q1',q2',Phi22);
title('true $\phi_2(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(2,3,5)
surf(q1,q2,phi2_est); hold on;
title('estimated $\phi_2(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(2,3,6)
Error = abs((phi2_est - Phi22)./(Phi22));
error = norm(phi2_est - Phi22);
surf(q1,q2, Error)
title('Error = $||\phi_2(x)-\hat \phi_2(x)||$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square
%% plots for paper with positive lamda2
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
% clim([-2,2])

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
% clim([-2,2])

subplot(2,4,5)
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
% clim([-2,2])

subplot(2,4,6)
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
% clim([-2,2])


%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end
