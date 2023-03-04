%% eigenfunctions for duffing system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
% linearization at (0,0) saddle

x = sym('x',[2;1]); 
delta = 0.5;
f = [x(2); + x(1) - delta*x(2) - x(1)^3 ]; 

eqb_point_saddle = [0 0];
A = eval(subs(jacobian(f),[x(1) x(2)],eqb_point_saddle));
[V,D,W] = eig(A);
D = diag(D);

% for real eig functions
eig_val1_saddle = real(D(1)); w1_saddle = W(:,1);
eig_val2_saddle = real(D(2)); w2_saddle = (W(:,2));

% define nonlinear part x_dot = Ax + fn(x)
fn_saddle = f - A*[x(1);x(2)];

% define matlab functions
g1_saddle = matlabFunction(w1_saddle'*fn_saddle,'vars',{x(1),x(2)});
g2_saddle = matlabFunction(w2_saddle'*fn_saddle,'vars',{x(1),x(2)});
f_saddle = matlabFunction(f);

%% setup path integral
dim = 1; % dimension for integraton (1 for scalar)

%define grid where eigenfunction is well defined
boundsx1 = -2; boundsx2 = 2;
boundsy1 = -2; boundsy2 = 2;

gridx = boundsx1:0.02:boundsx2; 
gridy = boundsy1:0.02:boundsy2; 
[q1_node,q2_node] = meshgrid(gridx,gridy);

% setup inital values 
x_0 = [q1_node(:),q2_node(:)]; 
phi1_node_real=[]; phi2_node_real=[];
phi1_node_imag = []; phi2_node_imag = [];

%% setup for saddle point at (0,0)
bounds = 2;
grid = -bounds:0.02:bounds; %define grid where eigenfunction is well defined
[q1_saddle,q2_saddle] = meshgrid(grid);

x_0 = [q1_saddle(:),q2_saddle(:)]; 
phi1_saddle=[];phi2_saddle=[];

% stable eigenfunction for saddle point
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 20];
    [t,x] = ode45(@(t,x)f_saddle(x(1),x(2)),tspan,x_0(i,:));
    
    % terminate t and x when they enter linear region
    idx = find(norm(x(end,:)-eqb_point_saddle)<1e-6);
    if(~isempty(idx))
        idx = idx(1);
        t = t(1:idx); x = x(1:idx,:);
    end
    % real valued eigenfunctions at saddle point (0,0)
    phi1_saddle = [phi1_saddle, w1_saddle'*x_0(i,:)' + trapz(t,exp(-eig_val1_saddle*t).*g1_saddle(x(:,1),x(:,2)),dim)];
end
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

% unstable eigenfunction for saddle point
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 5];
    [t,x] = ode45(@(t,x)f_saddle(x(1),x(2)),tspan,x_0(i,:));
    
    % terminate t and x when they enter linear region
    idx = find(norm(x(end,:)-eqb_point_saddle)<1e-6);
    if(~isempty(idx))
        idx = idx(1);
        t = t(1:idx); x = x(1:idx,:);
    end
    % real valued eigenfunctions at saddle point (0,0)
    phi2_saddle = [phi2_saddle, w2_saddle'*x_0(i,:)' + trapz(t,exp(-eig_val2_saddle*t).*g2_saddle(x(:,1),x(:,2)),dim)];
end
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);


% phi for saddle at (0,0)
phi1_saddle = reshape((phi1_saddle),size(q2_saddle)); %phi1_saddle(phi1_saddle>10)=10;
phi2_saddle = reshape((phi2_saddle),size(q2_saddle)); %phi2_saddle(phi2_saddle>10)=10;

%% plot only the saddle point
% saddle point at (0,0)
ic_pts = 5;
subplot(2,2,1)
p1 = pcolor(q1_saddle,q2_saddle,phi1_saddle); hold on;
set(p1,'Edgecolor','none')
colormap jet

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
title ('Stable Eigenfunction $\psi_1(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
colorbar

subplot(2,2,3)
p2 = pcolor(q1_saddle,q2_saddle,phi2_saddle); hold on;
set(p2,'Edgecolor','none')
colormap jet

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
title ('Unstable eigenfunction $\psi_2(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
colorbar