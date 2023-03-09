%% eigenfunctions for vanderpol system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
x = sym('x',[3;1]);  Dom = [-2 2];
% sigma = 10; rho= 28; beta = 8/3;
sigma = 1.0; rho= 2.8; beta = 8/30;
f = [sigma*(x(2)-x(1)); x(1)*(rho-x(3)) - x(2); x(1)*x(2) - beta*x(3)]; 

% get quiver
% [X,Y,Z] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
% u = sigma.*(Y-X);
% v = X.*(rho-Z) - Y;
% w = X.*Y - beta.*Z;

figure(1)
subplot(2,4,1)
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% l = streamslice(X,Y,Z,u,v,w); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');

ff = @(t,x)[sigma*(x(2)-x(1)); x(1)*(rho-x(3)) - x(2); x(1)*x(2) - beta*x(3)]; 
tspan = [0,10]; ic_pts = 5;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
zl = Dom(1); zh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    waitbar(i/length(x0),w_bar,sprintf(string(i)+'/'+string(length(x0))))
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        for z0 = linspace(Dom(1), Dom(2), ic_pts)
            [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0 z0]);
            plot3(xs(:,1),xs(:,2),xs(:,3),'k','LineWidth',1); hold on;
        end
    end
end
xlim([-3,3])
ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
zlabel('$x_3$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% linearization at (0,0,0) saddle
eqb_point = [0 0 0];
A = eval(subs(jacobian(f),[x(1) x(2) x(3)],eqb_point));
[V,D,W] = eig(A);
D = diag(D);

eig_val1 = D(1); 
eig_val2 = D(2); 
eig_val3 = D(3); 
w1_real = W(:,1);
w2_real = W(:,2);
w3_real = W(:,3);

% define nonlinear part x_dot = Ax + fn(x)
fn = f - A*[x(1);x(2);x(3)];

% define matlab functions
g1 = matlabFunction(w1_real'*fn,'vars',{x(1),x(2),x(3)});
g2 = matlabFunction(w2_real'*fn,'vars',{x(1),x(2),x(3)});
g3 = matlabFunction(w3_real'*fn,'vars',{x(1),x(2),x(3)});

f = matlabFunction(f);
%% setup path integral
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

dim = 1; % dimension for integraton (1 for scalar)

% %define grid where eigenfunction is well defined
% setup for eqb point at (0,0)
bounds = Dom(2);
gridx = -3:0.25:3; %define grid where eigenfunction is well defined
gridy = -3:0.25:3;
gridz = -2:0.25:6;
[q1,q2,q3] = meshgrid(gridx,gridy,gridz);
U = f(q1(:)',q2(:)', q3(:)');
x_0 = [q1(:),q2(:),q3(:)]; 
phi1=[]; phi2=[]; phi3=[];
for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 10];
    [t,x] = ode45(@(t,x)f(x(1),x(2),x(3)),tspan,x_0(i,:));
    
    % check if sol stays within the region of interest
    % if sol leaves this circle, dont compute path integral
    if norm(x(end,:)-eqb_point)>100
        phi1 = [phi1, -2];
        phi2 = [phi2, -2];
        phi3 = [phi3, -2];

    else
        % terminate t and x when they enter linear region
        idx = find(norm(x(end,:)-eqb_point)<1e-6);
        if(~isempty(idx))
            idx = idx(1);
            t = t(1:idx); x = x(1:idx,:);
        end

        % path integral for complex eigenfunctions
        % real part
        phi1 = [phi1, w1_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val1*t).*g1(x(:,1),x(:,2),x(:,3)),dim)];
    
        phi2 = [phi2, w2_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val2*t).*g2(x(:,1),x(:,2),x(:,3)),dim)];

        phi3 = [phi3, w3_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val3*t).*g3(x(:,1),x(:,2),x(:,3)),dim)];
    end
end
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% reshape to plot
% phi for eqb point at (0,0,0)
phi1 = reshape((phi1),size(q2));
phi2 = reshape((phi2),size(q2));
phi3 = reshape((phi3),size(q2));

%% plot
[t,x] = ode45(@(t,x)f(x(1),x(2),x(3)),tspan,[1,1,1]);
x_bin = discretize(x(:,1),gridx);
y_bin = discretize(x(:,2),gridy);
z_bin = discretize(x(:,3),gridz);
bin_idx = [x_bin,y_bin,z_bin];
bin = x(bin_idx);
figure(2)

subplot(1,3,1)
sz = 50;
plot3(x(:,1),x(:,2),x(:,3)); hold on
scatter3(bin(:,1),bin(:,2),bin(:,3),sz,phi1(bin_idx),'filled');

subplot(1,3,2)
plot3(x(1,:),x(2,:),x(3,:)); hold on
scatter3(bin(:,1),bin(:,2),bin(:,3),sz,phi2(bin_idx),'filled');

subplot(1,3,3)
plot3(x(1,:),x(2,:),x(3,:)); hold on
scatter3(bin(:,1),bin(:,2),bin(:,3),sz,phi3(bin_idx),'filled');
