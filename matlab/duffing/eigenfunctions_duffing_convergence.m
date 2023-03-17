% Path integral Eigenfunctions
%% eigenfunctions for duffing system
clc; clear; %close all;
%% system description
% nonlinear ode x_dot = f(x)
% linearization at (0,0) saddle
Dom = [-2 2];
x = sym('x',[2;1]); 
delta = 0.5;
f = -[x(2); + x(1) - delta*x(2) - x(1)^3]; 

% get quiver
grid = Dom(1):0.01:Dom(2);
[X,Y] = meshgrid(grid);
u = Y;
v = X - delta.*Y - X.^3;

% figure(1)
% subplot(2,4,1)
% ff = @(t,x)[x(2); + x(1) - delta*x(2) - x(1)^3]; 
% tspan = [0,20]; ic_pts = 2;
% xl = Dom(1); xh = Dom(2);
% yl = Dom(1); yh = Dom(2);
% Xs = [];
% start_idx = 1;
% for x0 = linspace(Dom(1), Dom(2), ic_pts)
%     for y0 = linspace(Dom(1), Dom(2), ic_pts)
%         [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
%         plot(xs(start_idx:end,1),xs(start_idx:end,2),'k','LineWidth',1); hold on;
%         Xs = [Xs;xs(start_idx:end,:)];
%     end
% end
% % xlim([-3,3])
% % ylim([-3,3])
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% zlabel('$x_3$','FontSize',20, 'Interpreter','latex')
% box on
% axes.LineWidth=2;

%% Linearize system
eqb_point_saddle = [0 0];
A = eval(subs(jacobian(f),[x(1) x(2)],eqb_point_saddle));
[V,D,W] = eig(A);
D = diag(D);

% for real eig functionsTh
eig_val1 = real(D(1)); w1 = W(:,1);
eig_val2 = real(D(2)); w2 = (W(:,2));

% define nonlinear part x_dot = Ax + fn(x)
fn = f - A*[x(1);x(2)];

% define matlab functions
g1 = matlabFunction(w1'*fn,'vars',{x(1),x(2)});
g2 = matlabFunction(w2'*fn,'vars',{x(1),x(2)});
f = matlabFunction(f);

%% Set up path Integral
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
dim = 1;
grid = Dom(1):0.025:Dom(2); %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
x_0 = [q1(:)';q2(:)'];

x_0 = [q1(:),q2(:)]; 
phi1=[];phi2=[];
phi1_time_ode45=[]; phi1_est_time_ode45 = [];
phi2_time_ode45=[]; phi2_est_time_ode45 = [];
options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));

for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))

    % iterate over time
    phi1_time_ode45 = []; phi2_time_ode45 = [];
    t_ode45 = linspace(0.001,10,5);
    for t_int = t_ode45
        [t,x] = ode45(@(t,x)f(x(1),x(2)),[0 t_int],x_0(i,:), options);
        %collect row vector of eig fun for each time
        phi1_time_ode45 = [phi1_time_ode45, w1'*x_0(i,:)'...
            + trapz(t,exp(-eig_val1*t).*g1(x(:,1),x(:,2)),dim)];
        phi2_time_ode45 = [phi2_time_ode45, w2'*x_0(i,:)'...
            + trapz(t,exp(-eig_val2*t).*g2(x(:,1),x(:,2)),dim)];
    end

    % collect columns of phi1_time for each initial condition
    phi1_est_time_ode45 = [phi1_est_time_ode45;phi1_time_ode45];
    phi2_est_time_ode45 = [phi2_est_time_ode45;phi2_time_ode45];

    % calculate phi1 phi2 for plot
    [t,x] = ode45(@(t,x)f(x(1),x(2)),[0 t_int],x_0(i,:), options);
    % stable eigenfunctions at saddle point (0,0)
    phi1 = [phi1, w1'*x_0(i,:)' + trapz(t,exp(-eig_val1*t).*g1(x(:,1),x(:,2)),dim)];
    %phi1 = [phi1, w1'*x_0(i,:)' + trapz(t,exp(-eig_val1*t).*g1(x(1),x(2))',dim)];

    % unstable eigenfunctions at saddle point (0,0)
    phi2 = [phi2, w2'*x_0(i,:)' + trapz(t,exp(-eig_val2*t).*g2(x(:,1),x(:,2)),dim)];
    %phi2 = [phi2, w2'*x_0(i,:)' + trapz(t,exp(-eig_val2*t).*g2(x(1),x(2))',dim)];
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% reshape
phi1 = reshape((phi1),size(q2));
phi2 = reshape((phi2),size(q2));
phi1_time_ode45 = phi1_est_time_ode45(:,end);
phi1_time_ode45 = reshape((phi1_time_ode45),size(q2));
phi2_time_ode45 = phi2_est_time_ode45(:,end);
phi2_time_ode45 = reshape((phi2_time_ode45),size(q2));

%% plot eigenfunctions
figure(6)

subplot(2,4,1)
p1 = pcolor(q1,q2,phi1); hold on;
set(p1,'Edgecolor','none')
colormap jet
%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
% xlim([-2,2])
% ylim([-2,2])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar
% clim([-5,5])


subplot(2,4,2)
[min_val, min_idx] = min(phi2,[],'all');
[r,c] = ind2sub(size(q1), min_idx);
phi2(r,c) = nan;
[max_val, max_idx] = max(phi2,[],'all');
[r,c] = ind2sub(size(q1), max_idx);
phi2(r,c) = nan;

p2 = pcolor(q1,q2,phi2); hold on;
set(p2,'Edgecolor','none')
colormap jet
%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
% xlim([-2,2])
% ylim([-2,2])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar
% clim([-5,5])

%% plot convergence rate
figure(2)
subplot(2,4,[5 6])
t_plot = t_ode45;
for i = 1:1:length(x_0)
    plot(t_plot,phi1_est_time_ode45(i,:)); hold on
end
axes = gca;
set(axes,'FontSize',15);
xlabel('integration time $(s)$','FontSize',20, 'Interpreter','latex')
ylabel('$\psi_s(x)$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

subplot(2,4,[7 8])
t_plot = t_ode45;
for i = 1:1:length(x_0)
    plot(t_plot,phi2_est_time_ode45(i,:)); hold on
end
axes = gca;
set(axes,'FontSize',15);
xlabel('integration time $(s)$','FontSize',20, 'Interpreter','latex')
ylabel('$\psi_u(x)$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end