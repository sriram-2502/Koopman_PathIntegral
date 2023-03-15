%% eigenfunctions for vanderpol system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
x = sym('x',[3;1]);  
Dom = [-0.02 0.02];
%Dom = [-2 2];
sigma = 10; rho= 28; beta = 8/3;
%sigma = 1.0; rho= 2.8; beta = 8/30;
f = [sigma*(x(2)-x(1)); x(1)*(rho-x(3)) - x(2); x(1)*x(2) - beta*x(3)]; 

% get quiver
% [X,Y,Z] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
% u = sigma.*(Y-X);
% v = X.*(rho-Z) - Y;
% w = X.*Y - beta.*Z;
% quiver3(X,Y,Z,u,v,w)

figure(1)
subplot(2,4,1)
ff = @(t,x)[sigma*(x(2)-x(1)); x(1)*(rho-x(3)) - x(2); x(1)*x(2) - beta*x(3)]; 
tspan = [0,20]; ic_pts = 2;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
zl = Dom(1); zh = Dom(2);
Xs = [];
start_idx = 1;
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        for z0 = linspace(Dom(1), Dom(2), ic_pts)
            [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0 z0]);
            plot3(xs(start_idx:end,1),xs(start_idx:end,2),xs(start_idx:end,3),'k','LineWidth',1); hold on;
            Xs = [Xs;xs(start_idx:end,:)];
        end
    end
end
% xlim([-3,3])
% ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
zlabel('$x_3$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

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
dim = 1; % dimension for integraton (1 for scalar)

% %define grid where eigenfunction is well defined
% setup for eqb point at (0,0)
bounds = Dom(2);
q1 = Xs(:,1); q2 = Xs(:,2); q3 = Xs(:,3);
x_0 = [q1(:),q2(:),q3(:)]; 
phi1=[]; phi2=[]; phi3=[];
tspan = [0,20];
parfor i = 1:length(x_0)
    [t,x] = ode45(@(t,x)f(x(1),x(2),x(3)),tspan,x_0(i,:));
    
    % path integral for complex eigenfunctions
    % real part
    phi2 = [phi2, w2_real'*x_0(i,:)' +...
        trapz(t,exp(-eig_val2*t).*g2(x(:,1),x(:,2),x(:,3)),dim)];
end

tspan = [-200 0];
parfor i = 1:length(x_0)    
    [t,x] = ode45(@(t,x)f(x(1),x(2),x(3)),tspan,x_0(i,:));
    % real part
    phi1 = [phi1, (w1_real'*x_0(i,:)' -...
        trapz(t,exp(-eig_val1*t).*g1(x(:,1),x(:,2),x(:,3)),dim))./1e-10];

    phi3 = [phi3, w3_real'*x_0(i,:)' -...
        trapz(t,exp(-eig_val3*t).*g3(x(:,1),x(:,2),x(:,3)),dim)];
end
% phi1(phi1>1e6)=1e6;
% phi1(phi1<-1e6)=-1e6;
phi3(phi3>1e6)=1e6;
phi3(phi3<-1e6)=-1e6;
%% reshape to plot
% phi for eqb point at (0,0,0)
phi1 = reshape((phi1),size(q2));
phi2 = reshape((phi2),size(q2));
phi3 = reshape((phi3),size(q2));

%% scatter plot
figure(1)
sz = 10; ic_pts = 5; alpha = 1;
subplot(2,2,2)
Xs = [];
scatter3(q1(:),q2(:),q3(:),sz,phi1(:),'filled','MarkerFaceAlpha',alpha); hold on;
% for x0 = linspace(Dom(1), Dom(2), ic_pts)
%     for y0 = linspace(Dom(1), Dom(2), ic_pts)
%         for z0 = linspace(Dom(1), Dom(2), ic_pts)
%             [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0 z0]);
%             plot3(xs(start_idx:end,1),xs(start_idx:end,2),xs(start_idx:end,3),'k','LineWidth',1); hold on;
%         end
%     end
% end
axes1 = gca;
axis square
% %axis([-bounds bounds -bounds bounds])
set(axes1,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
colorbar
box on
axes1.LineWidth=2;

subplot(2,2,3)
scatter3(q1(:),q2(:),q3(:),sz,phi2(:),'filled','MarkerFaceAlpha',alpha); hold on;
% for x0 = linspace(Dom(1), Dom(2), ic_pts)
%     for y0 = linspace(Dom(1), Dom(2), ic_pts)
%         for z0 = linspace(Dom(1), Dom(2), ic_pts)
%             [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0 z0]);
%             plot3(xs(start_idx:end,1),xs(start_idx:end,2),xs(start_idx:end,3),'k','LineWidth',1); hold on;
%         end
%     end
% end
axes2 = gca;
axis square
% %axis([-bounds bounds -bounds bounds])
set(axes2,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
colorbar
box on
axes2.LineWidth=2;

subplot(2,2,4)
scatter3(q1(:),q2(:),q3(:),sz,phi3(:),'filled','MarkerFaceAlpha',alpha); hold on;
% for x0 = linspace(Dom(1), Dom(2), ic_pts)
%     for y0 = linspace(Dom(1), Dom(2), ic_pts)
%         for z0 = linspace(Dom(1), Dom(2), ic_pts)
%             [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0 z0]);
%             plot3(xs(start_idx:end,1),xs(start_idx:end,2),xs(start_idx:end,3),'k','LineWidth',1); hold on;
%         end
%     end
% end
axes3 = gca;
axis square
% %axis([-bounds bounds -bounds bounds])
set(axes3,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
colorbar
box on
axes3.LineWidth=2;
