%% eigenfunctions for duffing system
clc; clear; %close all;
%% system description
% nonlinear ode x_dot = f(x)
%eqb_point_node = [-1 0];
eqb_point = -1;
shift = -eqb_point;
Dom = [-3 2];

x = sym('x',[2;1]); 
delta = 0.5;
f = [x(2); + x(1) - delta*x(2) - x(1)^3 ]; 
% shift eqb point to -1,0
f_shifted = [x(2); + x(1) + shift - delta*x(2) - (x(1)+shift).^3 ]; 

% get quiver
[X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
u = Y;
v = (X + shift - delta.*Y - (X+shift).^3);

figure(1)
subplot(2,4,1)
ff = @(t,x)[x(2); + x(1) + shift - delta*x(2) - (x(1)+shift).^3 ]; 
tspan = [0,20]; ic_pts = 10;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
Xs = [];
start_idx = 1;
for x0 = linspace(-1, -3, ic_pts)
    for y0 = -2
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        %plot(xs(start_idx:end,1),xs(start_idx:end,2),'k','LineWidth',1); hold on;
        Xs = [Xs;xs(start_idx:end,:)];
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


%% linearization at (-1,0) stable and complex
eqb_point = [0 0];
A = eval(subs(jacobian(f_shifted),[x(1) x(2)],eqb_point));
[V,D,W] = eig(A);
l1 = D(1,1); l2 = D(2,2);
w1 = W(:,1); w2 = W(:,2);
D = diag(D);

eig_val1_real = real(D(1)); eig_val1_imag = imag(D(1)); 
eig_val2_real = real(D(2)); eig_val2_imag = imag(D(2)); 
w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
w2_real = real(W(:,2)); w2_imag = imag(W(:,2));

% define nonlinear part x_dot = Ax + fn(x)
fn = f_shifted - A*[x(1);x(2)];

% define matlab functions
g1 = matlabFunction(w1'*fn,'vars',{x(1),x(2)});
g2 = matlabFunction(w2'*fn,'vars',{x(1),x(2)});
g1_real = matlabFunction(w1_real'*fn,'vars',{x(1),x(2)});
g1_imag = matlabFunction(w1_imag'*fn,'vars',{x(1),x(2)});
g2_real = matlabFunction(w2_real'*fn,'vars',{x(1),x(2)});
g2_imag = matlabFunction(w2_imag'*fn,'vars',{x(1),x(2)});
f_shifted = matlabFunction(f_shifted);
f = matlabFunction(f);

%% setup path integral
dim = 1; % dimension for integraton (1 for scalar)

%define grid where eigenfunction is well defined
boundsx1 = Dom(1); boundsx2 = Dom(2);
boundsy1 = Dom(1); boundsy2 = Dom(2);

gridx = boundsx1:0.01:boundsx2; 
gridy = boundsy1:0.01:boundsy2; 
[q1,q2] = meshgrid(gridx,gridy);

% setup inital values 
x_0 = [q1(:),q2(:)]; 
phi1_real=[]; phi2_real=[];
phi1_imag = []; phi2_imag = [];
phi1_matlab = []; phi2 = [];

%% setup for stable node at (-1,0)
% w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% pick points along trajectory
q1 = Xs(:,1); q2 = Xs(:,2);
x_0 = [q1(:),q2(:)];

parfor i = 1:length(x_0)
    %waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    %t_span = [0 100];
    tspan = [-200 0]; % negative time
    [t,x] = ode45(@(t,x)f_shifted(x(1),x(2)),tspan,x_0(i,:));
    
    % complex eigenfunctions at stable node (-1,0)
    % real part
    phi1_real = [phi1_real, w1_real'*x_0(i,:)' -...
        trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_real(x(:,1),x(:,2) ...
            +sin(eig_val1_imag*t)).*g1_imag(x(:,1),x(:,2))),dim)];

    % imaginary part
    phi1_imag = [phi1_imag, w1_imag'*x_0(i,:)' -...
        trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_imag(x(:,1),x(:,2)) ...
        -sin(eig_val1_imag*t).*g1_real(x(:,1),x(:,2))),dim)];

    % compute directly
    phi1_matlab = [phi1_matlab, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x(:,1),x(:,2)), dim)];

end

% F = findall(0,'type','figure','tag','TMWWaitbar');
% delete(F);

%% compute lyapunov lambda^T * P + P * Lambda < 0  
[V,D,W] = eig(A);
[V_realified,D_realified] = cdf2rdf(V,D);
Q = eye(2);
P = lyap(D_realified,Q);
V_lyap = (phi1_real.^2-phi1_imag.^2)*P(1,1) + (phi1_real.^2-phi1_imag.^2)*P(2,2);

%% reshape for plots
% phi for node at (-1,0)
phi1_real = reshape((phi1_real),size(q2));
phi1_imag = reshape((phi1_imag),size(q2));
V_lyap = reshape((V_lyap),size(q2));
phi1_matlab = reshape((phi1_matlab),size(q2));

phi1_mag  = sqrt(phi1_real.^2+phi1_imag.^2);
phi1_mag_matlab = sqrt(real(phi1_matlab).^2+imag(phi1_matlab).^2);

phi1_phase = angle(phi1_real + i*phi1_imag);
phi1_phase2 = angle(real(phi1_matlab) + i*imag(phi1_matlab));
phi1_mag_log = log(phi1_mag);

%% scatter plot
figure(1)
sz = 25; ic_pts = 5; alpha = 1;
subplot(2,4,5)
Xs = [];
scatter(q1(:),q2(:),sz,phi1_mag(:),'filled','MarkerFaceAlpha',alpha); hold on;
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
axes1 = gca;
xlim([-2-shift,2-shift])
ylim([-2,2])
axis square
% %axis([-bounds bounds -bounds bounds])
set(axes1,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
colorbar
box on
axes1.LineWidth=2;

subplot(2,4,6)
scatter(q1(:),q2(:),sz,V_lyap(:),'filled','MarkerFaceAlpha',alpha); hold on;
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color','k');
xlim([-2-shift,2-shift])
ylim([-2,2])
axes2 = gca;
axis square
% %axis([-bounds bounds -bounds bounds])
set(axes2,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
colorbar
box on
axes2.LineWidth=2;

%% plot pcolor
% figure(1)
% 
% % eigenfucntion 1 mag and phase
% subplot(2,4,2)
% [max_val, max_idx] = max(phi1_mag,[],'all');
% [r,c] = ind2sub(size(q1), max_idx);
% phi1_mag(r,c) = nan;
% p1 = pcolor(q1,q2,phi1_mag); hold on;
% set(p1,'Edgecolor','none')
% colormap jet
% 
% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
% xlim([-2-shift,2-shift])
% ylim([-2,2])
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% box on
% axes.LineWidth=2;
% colorbar
% 
% subplot(2,4,3)
% p2 = pcolor(q1,q2,phi1_phase); hold on;
% set(p2,'Edgecolor','none')
% colormap jet
% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
% xlim([-2-shift,2-shift])
% ylim([-2,2])
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% box on
% axes.LineWidth=2;
% colorbar
% 
% subplot(2,4,4)
% [min_val, max_idx] = min(V_lyap,[],'all');
% [r,c] = ind2sub(size(q1), max_idx);
% V_lyap(r,c) = nan;
% 
% [min_val, max_idx] = min(V_lyap,[],'all');
% [r,c] = ind2sub(size(q1), max_idx);
% V_lyap(r,c) = nan;
% 
% p2 = pcolor(q1,q2,V_lyap); hold on;
% set(p2,'Edgecolor','none')
% colormap jet
% l = streamslice(X,Y,u,v); hold on;
% set(l,'LineWidth',1)
% set(l,'Color','k');
% xlim([-2-shift,2-shift])
% ylim([-2,2])
% axes = gca;
% axis square
% set(axes,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% box on
% axes.LineWidth=2;
% colorbar

%% helper functions

function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-3);
isterminal=1;
direction=0;
end