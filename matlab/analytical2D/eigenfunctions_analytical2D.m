% Path integral Eigenfunctions
% Anaytical Example
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% set up the system for positive lambda2
Dom = [-2 2]; ds = 0.05;

lambda1=-25; lambda2=-10;

% convergence check
if(lambda2>0)
    if(2*lambda1-lambda2 < 0)
        disp('Convergece Condition is satified')
    else
        disp('Convergece Condition failed!!')
    end
end

if(lambda2<0)
    if(2*lambda1+lambda2<0)
        disp('Convergece Condition is satified')
    else
        disp('Convergece Condition failed!!')
    end
end

% sys 1
f = @(t, x) [lambda1*x(1,:); lambda2*(x(2,:)-x(1,:).^2)];
beta=lambda2/(2*lambda1-lambda2); % for x2 

alpha=1;
% get quiver
[X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
u = alpha.*(lambda1.*X);
v = alpha.*(lambda2.*(Y-X.^2));

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

grid = Dom(1):ds:Dom(2); %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
QQ = [q1(:)';q2(:)'];
U = f(0,QQ);
%quiver(q1(:),q2(:),U(1,:)',U(2,:)')

x_0 = [q1(:),q2(:)]; phi1_est=[];phi2_est=[];
%options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
options = odeset('RelTol',1e-9,'AbsTol',1e-300);
for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    % method 1
    [t,x] = ode45(@(t,x)f(t,x),[0 3],x_0(i,:), options);
    phi1_est = [phi1_est, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x')')];
    phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
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

%% Eigenfunctions
figure(1)
subplot(1,3,1)
phi1_est = reshape((phi1_est),size(q2));
surf(q1',q2',Phi22);
title('true $\phi_2(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(1,3,2)
surf(q1,q2,phi2_est); hold on;
title('estimated $\phi_2(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

subplot(1,3,3)
Error = abs((phi2_est - Phi22)./(Phi22));
error = norm(phi2_est - Phi22);
surf(q1,q2, Error)
title('Error = $||\phi_2(x)-\hat \phi_2(x)||$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20)
axis square

%% plots for paper with positive lamda2
figure(2)
subplot(2,4,1)
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
title('True: $\lambda_1=$'+string(lambda1)+' $\lambda_2=$'+string(lambda2), 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar
clim([-5,5])

subplot(2,4,2)
p2 = pcolor(q1,q2,phi2_est); hold on;
set(p2,'Edgecolor','none')
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
title('Estimated: $\lambda_1=$'+string(lambda1)+' $\lambda_2=$'+string(lambda2), 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar
clim([-5,5])

% %% set up the system for negative lambda2
% Dom = [-2 2]; ds = 0.05;
% 
% lambda1=-2; lambda2=-3;
% % sys 1 flow reversed and alpha=-1
% alpha=1;
% f = @(t, x) [alpha*lambda1*x(1,:); alpha*lambda2*(x(2,:)-x(1,:).^2)];
% beta=lambda2/(2*lambda1-lambda2); % for x2 
% % sys 2
% % f = @(t, x) [lambda1*x(1,:); lambda2*(x(2,:)+x(1,:).^2)];
% % beta=-lambda2/(2*lambda1-lambda2);
% 
% % get quiver for original system
% alpha=1;
% [X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
% u = alpha.*(lambda1.*X);
% v = alpha.*(lambda2.*(Y-X.^2));
% 
% % true eig functions
% phi1 = @(x) x(1,:);
% phi2 = @(x) x(2,:)+beta*x(1,:).^2;
% x = sym('x',[2;1], 'real');
% Phi = @(x)[phi1(x);phi2(x)];
% 
% %% linearize the system
% A = double(subs(jacobian(f(0,x),x),x,[0;0]));
% fn = @(x) f(0,x)-A*x;
% [V,D,W] = eig(A); 
% l1 = A(1,1); l2 = A(2,2);
% w1 = W(:,1); w2 = W(:,2);
% 
% g1 = @(x) w1'*fn(x);
% g2 = @(x) w2'*fn(x);
% 
% %% Set up path Integral
% w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% 
% grid = Dom(1):ds:Dom(2); %define grid where eigenfunction is well defined
% [q1,q2] = meshgrid(grid);
% QQ = [q1(:)';q2(:)'];
% U = f(0,QQ);
% %quiver(q1(:),q2(:),U(1,:)',U(2,:)')
% 
% x_0 = [q1(:),q2(:)]; phi1_est=[];phi2_est=[];
% options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
% for i = 1:length(x_0)
%     waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
%     % method 1
%     [t,x] = ode45(@(t,x)f(t,x),[0 100],x_0(i,:), options);
%     phi1_est = [phi1_est, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x')')];
%     phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
% end
% 
% F = findall(0,'type','figure','tag','TMWWaitbar');
% delete(F);
% 
% Phi11 = zeros(size(q1));
% Phi22 = zeros(size(q2));
% for i = 1:length(grid)
%     for j = 1:length(grid)
%         Phi11(i,j) = phi1([grid(i);grid(j)]) ;
%         Phi22(i,j) = phi2([grid(i);grid(j)]) ;
%     end
% end
% phi2_est = reshape((phi2_est),size(q2));
% 
% %% plots for paper with negative lambda2
% figure(2)
% subplot(2,4,7)
% p1 = pcolor(q1,q2,Phi22'); hold on;
% set(p1,'Edgecolor','none')
% colormap jet
% 
% %plot quiver
% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
% 
% xlim([-2,2])
% ylim([-2,2])
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% title('True: $\lambda_1=$'+string(lambda1)+' $\lambda_2=$'+string(lambda2), 'FontSize',18,'Interpreter','latex')
% box on
% axes.LineWidth=2;
% colorbar
% clim([-5,5])
% 
% subplot(2,4,8)
% p2 = pcolor(q1,q2,phi2_est); hold on;
% set(p2,'Edgecolor','none')
% colormap jet
% 
% %plot quiver
% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
% 
% xlim([-2,2])
% ylim([-2,2])
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% title('Estimated (reversed): $\lambda_1=$'+string(lambda1)+' $\lambda_2=$'+string(lambda2),'FontSize',18, 'Interpreter','latex')
% box on
% axes.LineWidth=2;
% colorbar
% clim([-5,5])

%% helper functions
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end

