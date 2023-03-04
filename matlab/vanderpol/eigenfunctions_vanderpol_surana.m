%% eigenfunctions for vanderpol system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
Dom = [-6 6]; ds = 0.05;
x = sym('x',[2;1]); 
mu = 1;
alpha = -1;
f = [alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))]; 

% get quiver
[X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
u = alpha.*Y;
v = alpha.*(mu.*Y - X - mu.*X.^2.*Y);

figure(1)
subplot(2,4,1)
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([-3,3])
ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

%% linearization at (0,0) unstable
eqb_point = [0 0];
A = eval(subs(jacobian(f),[x(1) x(2)],eqb_point));
[V,D,W] = eig(A);
D_diag = diag(D);

% split real and mag parts of eigval and eigvec
eig_val1_real = real(D_diag(1)); eig_val1_imag = imag(D_diag(1)); 
eig_val2_real = real(D_diag(2)); eig_val2_imag = imag(D_diag(2)); 
w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
w2_real = real(W(:,2)); w2_imag = imag(W(:,2));

% use complex eigval and real eigvec
[Wr,  Dr] = cdf2rdf(W ,D);
l1 = D(1,1); l2 = D(2,2);
w1_realified = Wr(:,1); w2_realified = Wr(:,2);

% define nonlinear part x_dot = Ax + fn(x)
fn = f - A*[x(1);x(2)];

% define matlab functions
g1_real = matlabFunction(w1_real'*fn,'vars',{x(1),x(2)});
g1_imag = matlabFunction(w1_imag'*fn,'vars',{x(1),x(2)});
g2_real = matlabFunction(w2_real'*fn,'vars',{x(1),x(2)});
g2_imag = matlabFunction(w2_imag'*fn,'vars',{x(1),x(2)});
g1_realified = matlabFunction(w1_realified'*fn,'vars',{x(1),x(2)});
g2_realified = matlabFunction(w2_realified'*fn,'vars',{x(1),x(2)});
f = matlabFunction(f);

%% setup path integral
dim = 1; % dimension for integraton (1 for scalar)

w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% %define grid where eigenfunction is well defined
% setup for eqb point at (0,0)
bounds = 3;
grid = -bounds:0.025:bounds; %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
U = f(q1(:)',q2(:)');

x_0 = [q1(:),q2(:)]; 
phi1_real=[]; phi2_real=[];
phi1_imag = []; phi2_imag = [];
phi1_est = []; phi2_est = [];
for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 10];
    %options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
    options = odeset('events',@(t, x)offFrame(t, x, Dom(2)));
    [t,x] = ode45(@(t,x)f(x(1),x(2)),tspan,x_0(i,:), options);

    % path integral for complex eigenfunctions
    % real part
    phi1_real = [phi1_real, w1_real'*x_0(i,:)' +...
        trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t)+sin(eig_val1_imag*t)).*g1_real(x(:,1),x(:,2)),dim)];

    %phi2_real = [phi2_real, w2_real'*x_0(i,:)' +...
    %    trapz(t,exp(-eig_val2_real*t).*(cos(eig_val2_imag*t)-sin(eig_val2_imag*t)).*g2_real(x(:,1),x(:,2)),dim)];

    % imaginary part
    phi1_imag = [phi1_imag, w1_real'*x_0(i,:)' +...
        -trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t)+sin(eig_val1_imag*t)).*g1_imag(x(:,1),x(:,2)),dim)];

    %phi2_imag = [phi2_imag, w2_real'*x_0(i,:)' +...
    %    trapz(t,exp(-eig_val2_real*t).*(cos(eig_val2_imag*t)+sin(eig_val2_imag*t)).*g2_imag(x(:,1),x(:,2)),dim)];
    
    % use full complex eig function and split real and imag using
    % Suarana's method
    phi1_est = [phi1_est, w1_realified'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1_realified(x(:,1),x(:,2)),dim)];
    %phi2_est = [phi2_est, w2_realified'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2_realified(x(:,1),x(:,2)),dim)];
end

        
% delete all waitbars
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% reshape
% phi for eqb point at (0,0) by splitting real and imag parts
phi1_real = reshape((phi1_real),size(q2));
phi1_imag = reshape((phi1_imag),size(q2));
phi1_mag  = sqrt(phi1_real.^2+phi1_imag.^2);
phi1_phase = angle(phi1_real + i*phi1_imag);
phi1_mag_log = log(phi1_mag);

% phi2_real = reshape((phi2_real),size(q2));
% phi2_imag = reshape((phi2_imag),size(q2));
% phi2_mag  = sqrt(phi2_real.^2+phi2_imag.^2);
% phi2_phase = angle(phi2_real + i*phi2_imag);
% phi2_mag_log = log(phi2_mag);

% use Surana to split real and imag
phi1_est_real = 2*real(phi1_est);
phi1_est_imag = -2*imag(phi1_est);
phi1_est_mag  = sqrt(phi1_est_real.^2+phi1_est_imag.^2);
phi1_est_phase = angle(phi1_est_real + i*phi1_est_imag);

% reshape
phi1_est_mag = reshape((phi1_est_mag),size(q2));
phi1_est_phase = reshape((phi1_est_phase),size(q2));

%% plot
figure(2)
% eigfun by splitting real and imag parts
% eigenfucntion 1 mag and phase
subplot(2,4,1)
p1 = pcolor(q1,q2,phi1_mag); hold on;
set(p1,'Edgecolor','none')
colormap jet

l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,100]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([-3,3])
ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(2,4,2)
p2 = pcolor(q1,q2,phi1_phase); hold on;
set(p2,'Edgecolor','none')
colormap jet

l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([-3,3])
ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

% eigfun mag by using Surana
subplot(2,4,3)
p2 = pcolor(q1,q2,phi1_est_mag); hold on;
set(p2,'Edgecolor','none')
colormap jet

l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([-3,3])
ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

% eigfun phase by using Surana
subplot(2,4,4)
p2 = pcolor(q1,q2,phi1_est_phase); hold on;
set(p2,'Edgecolor','none')
colormap jet

l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([-3,3])
ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

% eigenfunction 1 mag in log scale and phase
subplot(2,4,5)
p1 = pcolor(q1,q2,-phi1_mag_log); hold on;
set(p1,'Edgecolor','none')
colormap jet

l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([-3,3])
ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(2,4,6)
p2 = pcolor(q1,q2,phi1_phase); hold on;
set(p2,'Edgecolor','none')
colormap jet

l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([-3,3])
ylim([-3,3])
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
value = (max(abs(Y))>Dom) | (min(sum(abs(Y)))<1e-4);
isterminal=1;
direction=0;
end