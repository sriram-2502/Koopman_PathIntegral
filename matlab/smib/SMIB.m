%% eigenfunctions for vanderpol system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
ds = 0.05;
Dom = [-7,7];
Damp = 1.3; %damping
Pm = 5; % mechanical input
Pe = 10;
a = asin(Pm/Pe);
n= 2;
m = 1;
x = sym('x',[n;1]);
u = sym('u',[m;1]);
f = [x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm]; % syms form

% get quiver
% [X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
% u = alpha.*Y;
% v = alpha.*(mu.*Y - X - mu.*X.^2.*Y);

figure(1)
subplot(2,4,1)
% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
ff = @(t,x)[x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([Dom(1),Dom(2)])
ylim([Dom(1),Dom(2)])
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
D = diag(D);

eig_val1_real = real(D(1)); eig_val1_imag = imag(D(1)); 
eig_val2_real = real(D(2)); eig_val2_imag = imag(D(2)); 
w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
w2_real = real(W(:,2)); w2_imag = imag(W(:,2));

% define nonlinear part x_dot = Ax + fn(x)
fn = f - A*[x(1);x(2)];

% define matlab functions
g1_real = matlabFunction(w1_real'*fn,'vars',{x(1),x(2)});
g1_imag = matlabFunction(w1_imag'*fn,'vars',{x(1),x(2)});
g2_real = matlabFunction(w2_real'*fn,'vars',{x(1),x(2)});
g2_imag = matlabFunction(w2_imag'*fn,'vars',{x(1),x(2)});
f = matlabFunction(f);


% realified version
[V1,D1,W1] = eig(A);
[W, D] = cdf2rdf(W1, D1);
w1 = W(:,1); w2 = W(:,2);
g1 = matlabFunction(w1'*fn,'vars',{x(1),x(2)});
g2 = matlabFunction(w2'*fn,'vars',{x(1),x(2)});
lam_r = real(D1(1,1)); lam_i = imag(D1(1,1));
g_complex = matlabFunction(W1(:,1)'*fn,'vars',{x(1),x(2)});
%% setup path integral
dim = 1; % dimension for integraton (1 for scalar)

w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% %define grid where eigenfunction is well defined
% setup for eqb point at (0,0)
bounds = Dom(2);
grid = -bounds:0.05:bounds; %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
U = f(q1(:)',q2(:)');

x_0 = [q1(:),q2(:)]; 
phi1_real=[]; phi2_real=[];
phi1_imag = []; phi2_imag = [];
phi1_est = []; phi2_est = [];
Phi1_complex = [];
for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 60];
    options = odeset('events',@(t, x)offFrame(t, x, Dom(2)));
    [t,x] = ode45(@(t,x)f(x(1),x(2)),tspan,x_0(i,:), options);
    
    % check if sol stays within the region of interest
    % if sol leaves this circle, dont compute path integral
    if norm(x(end,:)-eqb_point)>0.5
        phi1_real = [phi1_real, -2];
        phi1_imag = [phi1_imag, -2];
        phi2_real = [phi2_real, -2];
        phi2_imag = [phi2_imag, -2];
        phi1_est = [phi1_est, -2];
        phi2_est = [phi2_est, -2];
        Phi1_complex = [Phi1_complex, -2];
    else
        % terminate t and x when they enter linear region
        idx = find(norm(x(end,:)-eqb_point)<1e-6);
        if(~isempty(idx))
            idx = idx(1);
            t = t(1:idx); x = x(1:idx,:);
        end

        % path integral for complex eigenfunctions
        % real part
        phi1_real = [phi1_real, w1_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_real(x(:,1),x(:,2) + ...
            sin(eig_val1_imag*t)).*g1_imag(x(:,1),x(:,2))),dim)];
    
        % imaginary part
        phi1_imag = [phi1_imag, w1_real'*x_0(i,:)' +...
            -trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_imag(x(:,1),x(:,2) - ...
            sin(eig_val1_imag*t)).*g1_real(x(:,1),x(:,2))),dim)];

        Phi1_complex = [Phi1_complex, W1(:,1)'*x_0(i,:)' + ...
            trapz(t,exp(-D1(1,1)*t).*g_complex(x(:,1),x(:,2)),dim)];
         
        % realified version
        if length(t)>1
            phi1_est = [phi1_est, w1'*x_0(i,:)' + ...
                trapz(t,exp(-lam_r*t).*cos(lam_i*t).*g1(x(:,1),x(:,2)),dim)];
            phi2_est = [phi2_est, w2'*x_0(i,:)' - ...
                trapz(t,exp(-lam_r*t).*sin(lam_i*t).*g2(x(:,1),x(:,2)),dim)];
        else
            phi1_est = [phi1_est, 0];
            phi2_est = [phi2_est, 0];
        end
    end
end

% normalize
% I_est1 = sqrt(trapz(grid,trapz(grid,reshape((phi1_est.^2),size(q2)),2)));
% I_est2 = sqrt(trapz(grid,trapz(grid,reshape((phi2_est.^2),size(q2)),2)));2;

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% reshape
% phi for eqb point at (0,0)
phi1_real = reshape((phi1_real),size(q2));
phi1_imag = reshape((phi1_imag),size(q2));
phi1_mag  = sqrt(phi1_real.^2+phi1_imag.^2);
phi1_phase = angle(phi1_real + i*phi1_imag);
phi1_mag_log = log(phi1_mag);
phi1_mag(phi1_mag>10)=10;

% realified
phi1_realified = reshape((phi1_est),size(q2));
phi2_realified = reshape((phi2_est),size(q2));
Phi_mag_realified  = sqrt(phi1_realified.^2+phi2_realified.^2);
Phi_phase_realified  = angle(phi1_realified+i*phi2_realified);
Phi_mag_realified(Phi_mag_realified>10)=10;

Phi1_complex = reshape(Phi1_complex,size(q2));
Phi1_complex_mag = abs(Phi1_complex);
Phi1_complex_mag(Phi1_complex_mag>10)=10;
%% plot
figure(1)

% eigenfucntion 1 mag and phase
subplot(2,4,3)
p1 = pcolor(q1,q2,phi1_mag); hold on;
set(p1,'Edgecolor','none')
colormap jet

% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
ff = @(t,x)[x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([Dom(1),Dom(2)])
ylim([Dom(1),Dom(2)])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(2,4,4)
p2 = pcolor(q1,q2,phi1_phase); hold on;
set(p2,'Edgecolor','none')
colormap jet

% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
ff = @(t,x)[x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([Dom(1),Dom(2)])
ylim([Dom(1),Dom(2)])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

% eigenfunction 1 mag in log scale and phase
subplot(2,4,7)
p1 = pcolor(q1,q2,Phi_mag_realified); hold on;
set(p1,'Edgecolor','none')
colormap jet

% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
ff = @(t,x)[x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([Dom(1),Dom(2)])
ylim([Dom(1),Dom(2)])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(2,4,8)
p2 = pcolor(q1,q2,Phi_phase_realified); hold on;
set(p2,'Edgecolor','none')
colormap jet

% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
ff = @(t,x)[x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([Dom(1),Dom(2)])
ylim([Dom(1),Dom(2)])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(2,4,5)
p2 = pcolor(q1,q2,Phi1_complex_mag); hold on;
set(p2,'Edgecolor','none')
colormap jet

% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
ff = @(t,x)[x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm];
tspan = [0,10]; ic_pts = 8;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
xlim([Dom(1),Dom(2)])
ylim([Dom(1),Dom(2)])
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
value = (max(abs(Y))>Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end