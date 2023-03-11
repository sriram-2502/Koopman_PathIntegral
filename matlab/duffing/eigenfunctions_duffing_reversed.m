%% eigenfunctions for duffing system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
% linearization at (0,0) saddle
Dom = [-2 2];
x = sym('x',[2;1]); 
delta = 0.5;
f = -[x(2); + x(1) - delta*x(2) - x(1)^3 ]; 

% get quiver
[X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
u = Y;
v = (X - delta.*Y - X.^3);

figure(1)
subplot(2,4,1)
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[x(2); + x(1) - delta*x(2) - x(1)^3]; 
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

%% Liearnize system

eqb_point_saddle = [0 0];
A = eval(subs(jacobian(f),[x(1) x(2)],eqb_point_saddle));
[V,D,W] = eig(A);
D = diag(D);

% for real eig functions
eig_val1 = real(D(1)); w1 = W(:,1);
eig_val2 = real(D(2)); w2 = (W(:,2));

% define nonlinear part x_dot = Ax + fn(x)
fn = f - A*[x(1);x(2)];

% define matlab functions
g1 = matlabFunction(w1'*fn,'vars',{x(1),x(2)});
g2 = matlabFunction(w2'*fn,'vars',{x(1),x(2)});
f = matlabFunction(f);

%% setup path integral for saddle point at (0,0)
dim = 1; % dimension for integraton (1 for scalar)
bounds = Dom(2);
grid = -bounds:0.025:bounds; %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);

x_0 = [q1(:),q2(:)]; 
phi1=[];phi2=[];

% stable eigenfunction for saddle point
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
%options = odeset('events',@(t, x)offFrame(t, x, Dom(2)));
%options = odeset('RelTol',1e-9,'AbsTol',1e-300);

for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 10];
    [t,x] = ode45(@(t,x)f(x(1),x(2)),tspan,x_0(i,:), options);
    
    % terminate t and x when they enter linear region
    idx = find(norm(x(end,:)-eqb_point_saddle)<1e-6);
    if(~isempty(idx))
        idx = idx(1);
        t = t(1:idx); x = x(1:idx,:);
    end

    % stable eigenfunctions at saddle point (0,0)
    phi1 = [phi1, w1'*x_0(i,:)' + trapz(t,exp(-eig_val1*t).*g1(x(:,1),x(:,2)),dim)];

    % unstable eigenfunctions at saddle point (0,0)
    phi2 = [phi2, w2'*x_0(i,:)' + trapz(t,exp(-eig_val2*t).*g2(x(:,1),x(:,2)),dim)];
end
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

% w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% for i = 1:length(x_0)
%     waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
%     tspan = [0 20];
%     %options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
%     %options = odeset('events',@(t, x)offFrame(t, x, Dom(2)));
%     options = odeset('RelTol',1e-9,'AbsTol',1e-300);
%     [t,x] = ode45(@(t,x)f(x(1),x(2)),tspan,x_0(i,:),options);
%     
%     % terminate t and x when they enter linear region
%     idx = find(norm(x(end,:)-eqb_point_saddle)<1e-6);
%     if(~isempty(idx))
%         idx = idx(1);
%         t = t(1:idx); x = x(1:idx,:);
%     end
% 
%     % unstable eigenfunctions at saddle point (0,0)
%     phi2 = [phi2, w2'*x_0(i,:)' + trapz(t,exp(-eig_val2*t).*g2(x(:,1),x(:,2)),dim)];
% end
% F = findall(0,'type','figure','tag','TMWWaitbar');
% delete(F);

%% resize
% phi for saddle at (0,0)
phi1 = reshape((phi1),size(q2)); %phi1_saddle(phi1_saddle>10)=10;
phi2 = reshape((phi2),size(q2)); %phi2_saddle(phi2_saddle>10)=10;

%% plot only the saddle point
% saddle point at (0,0)
ic_pts = 1;
figure(2)

subplot(2,2,1)
p1 = pcolor(q1,q2,phi1); hold on;
set(p1,'Edgecolor','none')
colormap jet
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
f = @(t,x)[x(2); -delta*x(2) - x(1)^3 + x(1)]; 
xl = -bounds; xh = bounds;
yl = -bounds; yh = bounds;
for x0 = linspace(-bounds, bounds, ic_pts)
    for y0 = linspace(-bounds, bounds, ic_pts)
        [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end

axes1 = gca;
axis square
axis([-bounds bounds -bounds bounds])
set(axes1,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
%title ('Stable Eigenfunction $\psi_1(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
colorbar
%clim([-5e10, 5e10])
box on
axes.LineWidth=2;

subplot(2,2,2)
p2 = pcolor(q1,q2,phi2); hold on;
set(p2,'Edgecolor','none')
colormap jet
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
f = @(t,x)[x(2); -delta*x(2) - x(1)^3 + x(1)]; 
xl = -bounds; xh = bounds;
yl = -bounds; yh = bounds;
for x0 = linspace(-bounds, bounds, ic_pts)
    for y0 = linspace(-bounds, bounds, ic_pts)
        [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end

axes2 = gca;
axis square
axis([-bounds bounds -bounds bounds])
set(axes2,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
%title ('Unstable eigenfunction $\psi_2(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
colorbar
%clim([-5e-4, 5e-4])
box on
axes.LineWidth=2;

%% zero level set
phi2_level_set = phi2;
phi2_level_set(abs(phi2_level_set)<10e-3)=0;
phi2_level_set(abs(phi2_level_set)>10e-3)=-100;

subplot(2,2,4)
p2 = pcolor(q1,q2,phi2_level_set); hold on;
set(p2,'Edgecolor','none')
colormap jet
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
f = @(t,x)[x(2); -delta*x(2) - x(1)^3 + x(1)]; 
xl = -bounds; xh = bounds;
yl = -bounds; yh = bounds;
for x0 = linspace(-bounds, bounds, ic_pts)
    for y0 = linspace(-bounds, bounds, ic_pts)
        [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end

axes2 = gca;
axis square
axis([-bounds bounds -bounds bounds])
set(axes2,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
%title ('Unstable eigenfunction $\psi_2(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
colorbar
clim([-5e-4, 5e-4])
box on
axes.LineWidth=2;

%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end