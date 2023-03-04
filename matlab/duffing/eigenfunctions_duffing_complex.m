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

%% linearization at (-1,0) stable and complex
eqb_point_node = [-1 0];

% change coordinates to set eqb point at 0,0
x = sym('x',[2;1]); 
delta = 0.5;
f1 = x(2)-eqb_point_node(2);
f2 = x(1)-eqb_point_node(1) -delta*(x(2)-eqb_point_node(2)) - (x(1)-eqb_point_node(1))^3;
f = [f1;f2]; 

A = eval(subs(jacobian(f),[x(1) x(2)],[0,0]));
[V,D,W] = eig(A);
D = diag(D);

% for imag eig functions
eig_val1_real = real(D(1)); eig_val1_imag = imag(D(1)); 
eig_val2_real = real(D(2)); eig_val2_imag = imag(D(2)); 
w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
w2_real = real(W(:,2)); w2_imag = imag(W(:,2));

% define nonlinear part x_dot = Ax + fn(x)
fn_node = f - A*[x(1);x(2)];

% define matlab functions
g1_real = matlabFunction(w1_real'*fn_node,'vars',{x(1),x(2)});
g1_imag = matlabFunction(w1_imag'*fn_node,'vars',{x(1),x(2)});
g2_real = matlabFunction(w2_real'*fn_node,'vars',{x(1),x(2)});
g2_imag = matlabFunction(w2_imag'*fn_node,'vars',{x(1),x(2)});
f_node = matlabFunction(f);

%% linearization at (1,0) stable and complex
% eqb_point_node = [1 0];
% A = eval(subs(jacobian(f),[x(1) x(2)],eqb_point_node));
% [V,D,W] = eig(A);
% D = diag(D);
% 
% % for imag eig functions
% eig_val1_real = real(D(1)); eig_val1_imag = imag(D(1)); 
% eig_val2_real = real(D(2)); eig_val2_imag = imag(D(2)); 
% w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
% w2_real = real(W(:,2)); w2_imag = imag(W(:,2));
% 
% % define nonlinear part x_dot = Ax + fn(x)
% fn_node = f - A*[x(1);x(2)];
% 
% % define matlab functions
% g1_real = matlabFunction(w1_real'*fn_node,'vars',{x(1),x(2)});
% g1_imag = matlabFunction(w1_imag'*fn_node,'vars',{x(1),x(2)});
% g2_real = matlabFunction(w2_real'*fn_node,'vars',{x(1),x(2)});
% g2_imag = matlabFunction(w2_imag'*fn_node,'vars',{x(1),x(2)});

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

%% setup for stable node at (-1,0)
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 20];
    [t,x] = ode45(@(t,x)f_node(x(1),x(2)),tspan,x_0(i,:));
    
    % check if sol stays within the region of interest
    % if sol leaves this circle, dont compute path integral
    if norm(x(end,:)-eqb_point_node)>0.7
        phi1_node_real = [phi1_node_real, -2];
        phi1_node_imag = [phi1_node_imag, -2];
        phi2_node_real = [phi2_node_real, -2];
        phi2_node_imag = [phi2_node_imag, -2];
    else
        idx = find(norm(x(end,:)-eqb_point_node)<1e-6);
        if(~isempty(idx))
            idx = idx(1);
            t = t(1:idx); x = x(1:idx,:);
        end

        % complex eigenfunctions at stable node (-1,0)
        % real part
        phi1_node_real = [phi1_node_real, w1_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t)+sin(eig_val1_imag*t)).*g1_real(x(:,1),x(:,2)),dim)];
    
        phi2_node_real = [phi2_node_real, w2_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val2_real*t).*(cos(eig_val2_imag*t)-sin(eig_val2_imag*t)).*g2_real(x(:,1),x(:,2)),dim)];
    
        % imaginary part
        phi1_node_imag = [phi1_node_imag, w1_real'*x_0(i,:)' +...
            -trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t)+sin(eig_val1_imag*t)).*g1_imag(x(:,1),x(:,2)),dim)];
    
        phi2_node_imag = [phi2_node_imag, w2_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val2_real*t).*(cos(eig_val2_imag*t)+sin(eig_val2_imag*t)).*g2_imag(x(:,1),x(:,2)),dim)];
    end
end
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

% reshape for plots
% phi for node at (-1,0)
phi1_node_real = reshape((phi1_node_real),size(q2_node));
phi1_node_imag = reshape((phi1_node_imag),size(q2_node));

% find magitude and phase
phi1_node_mag  = sqrt(phi1_node_real.^2+phi1_node_imag.^2);
phi1_node_phase  = angle(phi1_node_real + i*phi1_node_imag);

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

% reshape for plots
% normalize
% I_est1 = sqrt(trapz(grid,trapz(grid,reshape((phi1_est.^2),size(q2)),2)));
% I_est2 = sqrt(trapz(grid,trapz(grid,reshape((phi2_est.^2),size(q2)),2)));2;

% phi for saddle at (0,0)
phi1_saddle = reshape((phi1_saddle),size(q2_saddle)); %phi1_saddle(phi1_saddle>10)=10;
phi2_saddle = reshape((phi2_saddle),size(q2_saddle)); %phi2_saddle(phi2_saddle>10)=10;

%% plot only the saddle point
% % saddle point at (0,0)
% ic_pts = 5;
% subplot(2,2,1)
% p1 = pcolor(q1_saddle,q2_saddle,phi1_saddle); hold on;
% set(p1,'Edgecolor','none')
% colormap jet
% 
% f = @(t,x)[x(2); -delta*x(2) - x(1)^3 + x(1)]; 
% xl = -bounds; xh = bounds;
% yl = -bounds; yh = bounds;
% for x0 = linspace(-bounds, bounds, ic_pts)
%     for y0 = linspace(-bounds, bounds, ic_pts)
%         [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
%         plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
%     end
% end
% 
% axes1 = gca;
% axis square
% axis([-bounds bounds -bounds bounds])
% set(axes1,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% title ('Stable Eigenfunction $\psi_1(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
% colorbar
% 
% subplot(2,2,3)
% p2 = pcolor(q1_saddle,q2_saddle,phi2_saddle); hold on;
% set(p2,'Edgecolor','none')
% colormap jet
% 
% f = @(t,x)[x(2); -delta*x(2) - x(1)^3 + x(1)]; 
% xl = -bounds; xh = bounds;
% yl = -bounds; yh = bounds;
% for x0 = linspace(-bounds, bounds, ic_pts)
%     for y0 = linspace(-bounds, bounds, ic_pts)
%         [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
%         plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
%     end
% end
% 
% axes2 = gca;
% axis square
% axis([-bounds bounds -bounds bounds])
% set(axes2,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% title ('Unstable eigenfunction $\psi_2(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
% colorbar

%% plot all
% Node point at (-1,0)
set(gcf,'color','w');
ic_pts = 6;  % number of initial condictions

% convert to log
%phi1_node_mag = log(phi1_node_mag);

tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact');
for i=1:4  
    %subplot(4,8,[1,2,9,10])
    nexttile(1)
    p1 = pcolor(q1_node,q2_node,phi1_node_mag); hold on;
    set(p1,'Edgecolor','none')
    colormap jet
    
    % plot traj on top
    f = @(t,x)[x(2); -delta*x(2) - x(1)^3 + x(1)]; 
    xl = boundsx1; xh = boundsx2;
    yl = boundsy1; yh = boundsy2;
    for x0 = linspace(boundsx1, boundsx2, ic_pts)
        for y0 = linspace(boundsy1, boundsy2, ic_pts)
            [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
            plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
        end
    end
    
    axes1 = gca;
    axis([boundsx1 boundsx2 boundsy1 boundsy2])
    set(axes1,'FontSize',15);
    axis square
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    % title ('Stable Eigenfunction $|\psi_1(x)|$ at (-1,0)','FontSize',20, 'Interpreter','latex')
    colorbar
    
    nexttile(2)
    p2 = pcolor(q1_node,q2_node,phi1_node_phase); hold on;
    set(p2,'Edgecolor','none')
    colormap jet
    
    f = @(t,x)[x(2); -delta*x(2) - x(1)^3 + x(1)]; 
    xl = boundsx1; xh = boundsx2;
    yl = boundsy1; yh = boundsy2;
    for x0 = linspace(boundsx1, boundsx2, ic_pts)
        for y0 = linspace(boundsy1, boundsy2, ic_pts)
            [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
            plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
        end
    end
    
    axes2 = gca;
    axis([boundsx1 boundsx2 boundsy1 boundsy2])
    set(axes2,'FontSize',15);
    axis square
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    % title ('Stable Eigenfunction Phase $\psi_1(x)$ at (-1,0)','FontSize',20, 'Interpreter','latex')
    colorbar
    
    % saddle point at (0,0)
    nexttile(5)
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
    
    axes3 = gca;
    axis([-bounds bounds -bounds bounds])
    set(axes3,'FontSize',15);
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    axis square
    % title ('stable Eigenfunction $\psi_1(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
    colorbar
    
    %subplot(4,8,[19,20,27,28])
    nexttile(6)
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
    
    axes4 = gca;
    axis([-bounds bounds -bounds bounds])
    set(axes4,'FontSize',15);
    axis square
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    % title ('Unstable eigenfunction $\psi_2(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
    colorbar
end