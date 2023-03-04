% Path integral Eigenfunctions
% Anaytical Example
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% set up the system for positive lambda2
Dom = [-2 2]; ds = 0.05;

lambda1=-4; lambda2=-3;
% sys 1
f = @(t, x) [lambda1*x(1,:); lambda2*(x(2,:)-x(1,:).^2)];
beta=lambda2/(2*lambda1-lambda2); % for x2 
% sys 2
% f = @(t, x) [lambda1*x(1,:); lambda2*(x(2,:)+x(1,:).^2)];
% beta=-lambda2/(2*lambda1-lambda2);

% get quiver
[X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
u = (lambda1.*X);
v = (lambda2.*(Y-X.^2));

% true eig functions
phi1 = @(x) x(1,:);
phi2 = @(x) x(2,:)+beta*x(1,:).^2;
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
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

grid = Dom(1):0.1:Dom(2); %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
QQ = [q1(:)';q2(:)'];
U = f(0,QQ);
%quiver(q1(:),q2(:),U(1,:)',U(2,:)')

x_0 = [q1(:),q2(:)]; 
phi1_est=[];phi2_est=[];
phi2_time_euler=[]; phi2_est_time_euler = [];
phi2_time_ode45=[]; phi2_est_time_ode45 = [];
options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));

for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    
    % do euler itegration
    dt = 1e-4;
    x0 = x_0(i,:);
    x = x0'; dim=1;
    phi2_time_euler=[];N=1e4;
    
    for euler_steps = 1:N
        time_elapsed = euler_steps*dt;
        x_new = x + dt*f(0,x);
        % collect row vector of eig fun for each time
        phi2_time_euler =[phi2_time_euler, w2'*x0' + trapz(time_elapsed,exp(-l2*time_elapsed).*g2(x_new)',dim)];
        x = x_new;
    end
    % collect columns of phi1_time for each initial condition
    phi2_est_time_euler = [phi2_est_time_euler;phi2_time_euler];

    % usde ode45
    [t,x] = ode45(@(t,x)f(t,x),[0 100],x_0(i,:), options);
    phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];

end

for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    % iterate over time
    phi2_time_ode45 = [];
    t_ode45 = 0.001:0.1:2;
    for t = t_ode45
        % usde ode45
        [t,x] = ode45(@(t,x)f(t,x),[0 t],x_0(i,:), options);
        %collect row vector of eig fun for each time
        phi2_time_ode45 = [phi2_time_ode45, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
    end
    % collect columns of phi1_time for each initial condition
    phi2_est_time_ode45 = [phi2_est_time_ode45;phi2_time_ode45];
end


F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

Phi11 = zeros(size(q1));
Phi22 = zeros(size(q2));
for i = 1:length(grid)
    for j = 1:length(grid)
        Phi11(i,j) = phi1([grid(i);grid(j)]) ;
        Phi22(i,j) = phi2([grid(i);grid(j)]) ;
    end
end
phi2_est = reshape((phi2_est),size(q2));
phi2_est_euler = phi2_est_time_euler(:,end);
phi2_est_euler = reshape((phi2_est_euler),size(q2));
phi2_time_ode45 = phi2_est_time_ode45(:,end);
phi2_time_ode45 = reshape((phi2_time_ode45),size(q2));

%% Convergence plot
figure(1)

subplot(2,2,1)
p1 = pcolor(q1,q2,Phi22'); hold on;
set(p1,'Edgecolor','none')
colormap jet

%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlim([-2,2])
ylim([-2,2])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar
clim([-5,5])


subplot(2,2,2)
p1 = pcolor(q1,q2,phi2_est); hold on;
set(p1,'Edgecolor','none')
colormap jet

%plot quiver
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlim([-2,2])
ylim([-2,2])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar
clim([-5,5])


subplot(2,2,[3,4])
t_plot = t_ode45;
for i = 1:20:length(x_0)
    plot(t_plot,phi2_est_time_ode45(i,:)); hold on
end

axes = gca;
%axis square
set(axes,'FontSize',15);
xlabel('integration time $(s)$','FontSize',20, 'Interpreter','latex')
ylabel('$h_u(x)$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
xlim([0,1.9])

% subplot(2,2,4)
% t_plot = (1:N).*dt;
% for i = 1:length(x_0)
%     plot(t_plot,phi2_est_time_euler(i,:)); hold on
% end
% 
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% box on
% axes.LineWidth=2;

%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end

