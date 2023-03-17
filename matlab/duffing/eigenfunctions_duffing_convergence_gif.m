% Path integral Eigenfunctions
%% eigenfunctions for duffing system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
% linearization at (0,0) saddle
Dom = [-2 2];
x = sym('x',[2;1]); 
delta = 0.5;
f = -[x(2); + x(1) - delta*x(2) - x(1)^3]; 

% get quiver
grid = Dom(1):0.5:Dom(2);
[X,Y] = meshgrid(grid);
u = Y;
v = X - delta.*Y - X.^3;

% % figure(1)
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

grid = Dom(1):0.1:Dom(2); %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
dim = 1;

x_0 = [q1(:),q2(:)]; 
phi1_est=[];phi2_est=[];
phi1=[]; phi1_est_time_ode45 = [];
phi2=[]; phi2_est_time_ode45 = [];
options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));

t_ode45 = linspace(0.001,20,5);
figure(2)
for t_span = t_ode45
    phi1 = []; phi2 = [];
    w_bar = waitbar(0,'1','Name','Path integral for t ='+string(t_span)+'s ...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    for i = 1:length(x_0)
        waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
        [t,x] = ode45(@(t,x)f(x(1),x(2)),[0 t_span],x_0(i,:), options);
        
        % terminate t and x when they enter linear region
        idx = find(norm(x(end,:)-eqb_point_saddle)<1e-6);
        if(~isempty(idx))
            idx = idx(1);
            t = t(1:idx); x = x(1:idx,:);
        end

        %collect row vector of eig fun
        % stable eigenfunctions at saddle point (0,0)
        phi1 = [phi1, w1'*x_0(i,:)' + trapz(t,exp(-eig_val1*t).*g1(x(:,1),x(:,2)),dim)];
    
        % unstable eigenfunctions at saddle point (0,0)
        phi2 = [phi2, w2'*x_0(i,:)' + trapz(t,exp(-eig_val2*t).*g2(x(:,1),x(:,2)),dim)];
    end

    % reshape for plot
    phi1 = reshape((phi1),size(q2));
    phi2 = reshape((phi2),size(q2));
    
    % plot eignfuns as a gif for each time
    ic_pts = 1;    
    subplot(2,2,1)
    p1 = pcolor(q1,q2,phi1); hold on;
    set(p1,'Edgecolor','none')
    colormap jet
    l = streamslice(X,Y,u,v); hold on;
    set(l,'LineWidth',1)
    set(l,'Color','k');
    axes1 = gca;
    axis square
    axis([Dom(1) Dom(2) Dom(1) Dom(2)])
    set(axes1,'FontSize',15);
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    colorbar
    phi1_min = min(phi1,[],'all'); phi1_max = max(phi1,[],'all');
    clim([phi1_min,phi1_max])
    box on
    axes.LineWidth=2;

    subplot(2,2,2)
    p2 = pcolor(q1,q2,phi2); hold on;
    set(p2,'Edgecolor','none')
    colormap jet
    l = streamslice(X,Y,u,v); hold on;
    set(l,'LineWidth',1)
    set(l,'Color','k');
    axes2 = gca;
    axis square
    axis([Dom(1) Dom(2) Dom(1) Dom(2)])
    set(axes2,'FontSize',15);
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    colorbar
    phi2_min = min(phi2,[],'all'); phi2_max = max(phi2,[],'all');
    clim([phi2_min,phi2_max])
    box on
    axes.LineWidth=2;
    
    % collect columns of phi1_time for each time
    phi1_est_time_ode45 = [phi1_est_time_ode45;phi1];
    phi2_est_time_ode45 = [phi2_est_time_ode45;phi2];

    F = findall(0,'type','figure','tag','TMWWaitbar');
    delete(F);
end

%% Convergence plot
figure(2)

subplot(2,2,3)
t_plot = t_ode45;
phi1_convergence = reshape(phi1_est_time_ode45,5,[]);
for i = 1:1:length(x_0)
    plot(t_plot,phi1_convergence(:,i)); hold on
end
axes = gca;
set(axes,'FontSize',15);
xlabel('integration time $(s)$','FontSize',20, 'Interpreter','latex')
ylabel('$\psi_s(x)$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

subplot(2,2,4)
phi2_convergence = reshape(phi2_est_time_ode45,5,[]);
t_plot = t_ode45;
for i = 1:1:length(x_0)
    plot(t_plot,phi2_convergence(:,i)); hold on
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

