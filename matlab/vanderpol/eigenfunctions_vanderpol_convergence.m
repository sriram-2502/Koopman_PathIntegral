%% eigenfunctions for vanderpol system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
Dom = [-3 3]; ds = 0.07;
x = sym('x',[2;1]); 
mu = 1;
alpha = 1;
f = [alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))]; 

% get quiver
[X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
u = alpha.*Y;
v = alpha.*(mu.*Y - X - mu.*X.^2.*Y);

%figure(1)
%subplot(2,4,1)
% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
% ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
% tspan = [0,10]; ic_pts = 8;
% xl = Dom(1); xh = Dom(2);
% yl = Dom(1); yh = Dom(2);
% for x0 = linspace(Dom(1), Dom(2), ic_pts)
%     for y0 = linspace(Dom(1), Dom(2), ic_pts)
%         [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
%         plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
%     end
% end
% xlim([-3,3])
% ylim([-3,3])
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% box on
% axes.LineWidth=2;

%% linearization at (0,0) unstable
eqb_point = [0 0];
A = eval(subs(jacobian(f),[x(1) x(2)],eqb_point));
[V,D,W] = eig(A);
l1 = D(1,1); l2 = D(2,2);
w1 = W(:,1); w2 = W(:,2);
D = diag(D);

eig_val1_real = real(D(1)); eig_val1_imag = imag(D(1)); 
eig_val2_real = real(D(2)); eig_val2_imag = imag(D(2)); 
w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
w2_real = real(W(:,2)); w2_imag = imag(W(:,2));

% define nonlinear part x_dot = Ax + fn(x)
fn = f - A*[x(1);x(2)];

% define matlab functions
g1 = matlabFunction(w1'*fn,'vars',{x(1),x(2)});
g2 = matlabFunction(w2'*fn,'vars',{x(1),x(2)});
g1_real = matlabFunction(w1_real'*fn,'vars',{x(1),x(2)});
g1_imag = matlabFunction(w1_imag'*fn,'vars',{x(1),x(2)});
g2_real = matlabFunction(w2_real'*fn,'vars',{x(1),x(2)});
g2_imag = matlabFunction(w2_imag'*fn,'vars',{x(1),x(2)});
f = matlabFunction(f);

%% Set up path Integral
% w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

dim = 1;
grid = Dom(1):0.025:Dom(2); %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
x_0 = [q1(:)';q2(:)'];

x_0 = [q1(:),q2(:)]; 
phi1_ode45_time_real = [];
phi1_ode45_time_imag = [];
phi1_real=[]; phi2_real=[];
phi1_imag = []; phi2_imag = [];
phi1_matlab = []; phi2 = [];

options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
time_to_simulate = 5;
t_ode45 = linspace(0.001,50,time_to_simulate);

parfor i = 1:length(x_0)
    %waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    % iterate over time
    phi1_ode45_real = []; phi1_ode45_imag = [];
    for t_int = t_ode45
        [t,x] = ode45(@(t,x)f(x(1),x(2)),[0 t_int],x_0(i,:), options);
        
        %collect row vector of eig fun for each time
        phi1_ode45_real = [phi1_ode45_real, w1_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_real(x(:,1),x(:,2) ...
            +sin(eig_val1_imag*t)).*g1_imag(x(:,1),x(:,2))),dim)];
        
        phi1_ode45_imag = [phi1_ode45_imag, w1_imag'*x_0(i,:)' +...
            trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_imag(x(:,1),x(:,2)) ...
            -sin(eig_val1_imag*t).*g1_real(x(:,1),x(:,2))),dim)];   
    end
    
    % collect columns of phi1_time for each initial condition
    phi1_ode45_time_real = [phi1_ode45_time_real;phi1_ode45_real];
    phi1_ode45_time_imag = [phi1_ode45_time_imag;phi1_ode45_imag];

    % calculate phi1 phi2 for plot
    [t,x] = ode45(@(t,x)f(x(1),x(2)),[0 t_int],x_0(i,:), options);

    phi1_real = [phi1_real, w1_real'*x_0(i,:)' +...
    trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_real(x(:,1),x(:,2) ...
        +sin(eig_val1_imag*t)).*g1_imag(x(:,1),x(:,2))),dim)];
    
    phi1_imag = [phi1_imag, w1_imag'*x_0(i,:)' +...
    trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_imag(x(:,1),x(:,2)) ...
        -sin(eig_val1_imag*t).*g1_real(x(:,1),x(:,2))),dim)];
end

% F = findall(0,'type','figure','tag','TMWWaitbar');
% delete(F);

%% reshape
phi1_real = reshape((phi1_real),size(q2));
phi1_imag = reshape((phi1_imag),size(q2));

phi1_mag  = sqrt(phi1_real.^2+phi1_imag.^2);
phi1_phase = angle(phi1_real + i*phi1_imag);

phi1_ode45_time_mag = sqrt(phi1_ode45_time_real.^2+phi1_ode45_time_imag.^2);

%% plot eigenfunctions
figure(3)
ic_pts = 1; Dom = [-0.1, 0.1]; tspan = [0,100];
% eigenfucntion 1 mag and phase
subplot(2,4,1)
p1 = pcolor(q1,q2,phi1_mag); hold on;
set(p1,'Edgecolor','none')
colormap jet

l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
ff = @(t,x)[alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
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

%% Convergence plot
figure(2)

subplot(2,2,[1,2])
t_plot = t_ode45;
phi1_convergence = reshape(phi1_ode45_time_mag,time_to_simulate,[]);
for i = 1:1:length(x_0)
    plot(t_plot,phi1_convergence(:,i),'LineWidth',2); hold on
end
axes = gca;
set(axes,'FontSize',15);
xlabel('integration time $(s)$','FontSize',20, 'Interpreter','latex')
ylabel('$\psi_s(x)$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

%% helper functions

function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end