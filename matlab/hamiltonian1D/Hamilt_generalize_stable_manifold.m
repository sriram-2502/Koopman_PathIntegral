%% Stable Manifold for Hamiltonian system


clear all; close all; clc;
set(0,'DefaultLineLineWidth',2.5) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

f = @(t,x,u)(x - x.^3 + u); h = @(t,x,u)(x);

n=1; m=1;% dimension
x = sym('x',[n;1],'real');
p = sym('p',[n;1],'real'); u = sym('u',[m,1],'real');
f_sym = f(0,x,0); 
g = jacobian(f(0,x,u),u); 
fz1 = @(x)f(0,x,0);
A = double(subs(jacobian(f(0,x,0),x),x,0));
B = double(subs(jacobian(f(0,0,u),u),u,0));

%% 
qx = x^2 ;
C= 1;
Df  = jacobian(f_sym ,x);
R_x = g*inv(C)*g'
D_R_x = jacobian(p'*R_x*p, x);
Dq  = jacobian(qx,x);

H_sym = [f_sym - g*inv(C)*g'*p;
        -Df'*p+0.5*D_R_x-0.5*Dq']

z = [x;p];
eqp_u = zeros(n+m,1)
E = double(subs(jacobian(H_sym,z),z,eqp_u));
H_sys = matlabFunction(H_sym, 'Vars',{x,p});

%%
Dom = [-3 3];

[X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));

H_sym = H_sys(x,p)
H_u = H_sym(1)
H_v = H_sym(2)
H_u_fun =matlabFunction(H_u, 'Vars',{x,p});
H_v_fun =matlabFunction(H_v, 'Vars',{x,p});

figure(1); % Phase portrait
u =H_u_fun(X,Y);
v = H_v_fun(X,Y);
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');

eqb_saddle = [0 0];
[V,D,W] = eig(E);
D = diag(D);

eig_val1 = real(D(1)); w1 = W(:,1);
eig_val2 = real(D(2)); w2 = (W(:,2));

% fn = f - A*[x(1);x(2)];
fn = H_sys(x,p) - E*[x;p];

% define matlab functions
g1 = matlabFunction(w1'*fn,'vars',{x,p});
g2 = matlabFunction(w2'*fn,'vars',{x,p});
f = H_sys;

%%
num_points = 500;
delta_x =0.1;
x_min = eqb_saddle(1)-delta_x; % Minimum x coordinate of the rectangle
x_max = eqb_saddle(1)+delta_x; % Maximum x coordinate of the rectangle
y_min = eqb_saddle(2)-delta_x; % Minimum y coordinate of the rectangle
y_max = eqb_saddle(2)+delta_x; % Maximum y coordinate of the rectangle

% Generate random coordinates for the data points within the specified rectangle
x_coords = x_min + (x_max - x_min) * rand(num_points, 1);
y_coords = y_min + (y_max - y_min) * rand(num_points, 1);
data_points = [x_coords, y_coords];

%% Time Domain simulation Backward in time
options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
tspan = 0:0.1:20;
figure(11);
hold on;
x_all = [];
for i = 1: length(data_points)
    [t,x] = ode45(@(t,x)(-f(x(1),x(2))),tspan,data_points(i,:), options);
    x_all = [x_all ; x];
    plot(x(:,1),x(:,2))
end
%%
% Plot the data points

%%% Select the data points in the domain 
mag_lim =3;
selected_indices = abs(x_all(:, 1)) < mag_lim & abs(x_all(:, 2)) < mag_lim; %absolute value in the path integral domain
selected_data_points = x_all(selected_indices, :);
dim = 1; % dimension for integraton (1 for scalar)
x_0n = selected_data_points;
figure(1)
p2 = scatter(x_0n(:,1), x_0n(:,2));
p2.MarkerEdgeColor = [  0.3922    0.8314    0.0745];
p2.MarkerFaceColor = [0.4660 0.1740 0.5880];
p2.MarkerFaceAlpha = 0.4;
xlabel('$x_1$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
xlim(Dom);
ylim(Dom);
set(gca,'fontsize',16,'FontName','Times New Roman');
box on
axis square

%%
figure(1);
p1 = scatter(data_points(:, 1), data_points(:, 2)); hold on
p1.MarkerFaceColor =   [  0     0     1];
p1.MarkerEdgeColor = [ 0     1     1];
p1.MarkerFaceAlpha = 0.4;
sz = 100;
xlabel('$x_1$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
xlim(Dom)
ylim(Dom)
axis square;
set(gca,'fontsize',16,'FontName','Times New Roman');
box on

%% path integral for unstable eigenfunction 
% with initial points around stable manifold
phi2n=[];
options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
tspan = [0 2]; 

parfor k = 1:length(x_0n)
    [t,x] = ode45(@(t,x)f(x(1),x(2)),tspan,x_0n(k,:), options);
    phi2n = [phi2n, w1'*x_0n(k,:)' + trapz(t,exp(-eig_val1*t).*g1(x(:,1),x(:,2)),dim)];% unstable eigenfunctions at saddle point (0,0)
end


%% plot for path integral near boundary
dim = 1; % dimension for integraton (1 for scalar)

figure(3);hold on;
scatter3(x_0n(:,1), x_0n(:,2), phi2n, 50, phi2n, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
colormap jet;
colorbar('southoutside');
view(2)
xlabel('$x_1$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', 16, 'Interpreter', 'latex');
set(gca,'fontsize',20,'FontName','Times New Roman');
box on
axis square;
xlim(Dom);
ylim(Dom);
box on;
grid on;

%% ZERO LEVEL NEW TRIAL
% Define the range for phi2n values
phi2n_min = -0.001;
phi2n_max = 0.001;

% Find the indices of phi2n values within the specified range
indices = abs(phi2n) >= phi2n_min & abs(phi2n) <= phi2n_max;

% Filter the data points and corresponding phi2n values
filtered_x_0n = x_0n(indices, :);
filtered_phi2n = phi2n(indices);
%
% Plot the filtered data points with color based on phi2n values
figure(4);
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
p3 = scatter3(filtered_x_0n(:, 1), filtered_x_0n(:, 2), filtered_phi2n, 20, filtered_phi2n, 'filled');hold on;
p3.MarkerEdgeColor ='r';
p3.MarkerFaceColor = 'r';
view(2)
xlabel('$x_1$','FontSize',16, 'Interpreter','latex')
ylabel('$x_2$','FontSize',16, 'Interpreter','latex')
axes = gca;
axis square
set(gca,'fontsize',20,'FontName','Times New Roman');
grid on;

%% Least square
n=2;
n_x=n/2;
x = sym('x',[n;1]);

Phi =phi2n;
deg =13; %order basis (for monomial)
[Psi, ~] = monomial_basis(deg, n_x);
Psi_fun = matlabFunction(Psi,'Vars',{x});

Psi_grid = Psi_fun(x_0n');
Q= Phi*pinv(Psi_grid);
Phi_approx=  Q* Psi;
Phi_fun = matlabFunction(Phi_approx,'Vars',{x});

%Plots
Dom_n = Dom;
bounds_n = Dom_n(2); ds = 0.01;
grid_n = -bounds_n:ds:bounds_n; %define grid where eigenfunction is well defined
[q1_n,q2_n] = meshgrid(grid_n);
[X,Y] = meshgrid(Dom_n(1):0.25:Dom_n(2),Dom_n(1):0.25:Dom_n(2));
phi1_n = zeros(size(q1_n));
for i = 1:length(grid_n)
    for j = 1:length(grid_n)
        PHI = Phi_fun([grid_n(i);grid_n(j)]);
        phi1_n(i,j) = PHI(1);
    end
end

ic_pts = 1;
close(figure(11))
figure(11)
p1 = pcolor(q1_n,q2_n,phi1_n'.*1e4); hold on;
set(p1,'Edgecolor','none')
colormap jet
set(gca,'fontsize',20,'FontName','Times New Roman')
fimplicit(Phi_approx,'r',LineWidth=4)
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlabel('$x_1$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', 16, 'Interpreter', 'latex');
set(gca,'fontsize',16,'FontName','Times New Roman');
box on
axis square;
xlim(Dom);
ylim(Dom);
box on; 
grid on;

%% 
W1 = Q(:,n_x)
W2 = Q(:,n_x+1:end)
p = inv(W1)*W2*Psi(n_x+1:end);

x_dot = x(1)-x(1).^3- p;

dx = matlabFunction(x_dot,'Vars',{x(1)});

%% Euler integration for x_dot = f(x)
num_points = 5000;
x_euler = zeros(1, num_points);
x_euler(:, 1) = 3;
dt =0.1;
tspan = 0:dt:5000;
for i = 2:num_points
    x_euler(:, i) = x_euler(:, i-1) + dt * dx( x_euler(:, i-1));
end

figure; 
plot( x_euler,LineWidth=3 );
xlabel('$t$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$x(t)$', 'FontSize', 16, 'Interpreter', 'latex');
set(gca,'fontsize',16,'FontName','Times New Roman');
box on
axis square;
box on; 
grid on;
%% function value for scalar example
xlin  =-2:0.1:2;
dx_val = dx(xlin);
figure;
plot(xlin,dx_val,LineWidth=3 )
xlabel('$x$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$f(x)$', 'FontSize', 16, 'Interpreter', 'latex');
set(gca,'fontsize',16,'FontName','Times New Roman');
box on
axis square;
box on; 
grid on;
%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-3);
isterminal=1;
direction=0;
end