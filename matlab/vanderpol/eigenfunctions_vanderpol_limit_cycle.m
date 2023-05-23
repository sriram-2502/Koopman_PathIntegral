%% eigenfunctions for vanderpol system with unstable origin
clc; clear; close all;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])
%% system description
% nonlinear ode x_dot = f(x)
Dom = [-1 1];
x = sym('x',[2;1]); 
mu = 1;
alpha = 1;
f = [alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))]; 

% get quiver
[X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
u = alpha.*Y;
v = alpha.*(mu.*Y - X - mu.*X.^2.*Y);


ff = @(t,x)[-alpha*x(2); -alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))];
tspan = [0,10]; ic_pts = 10;
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
Xs = [];
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)ff(t,x),tspan,[x0 y0]);
        Xs = [Xs;xs];
    end
end

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

%% setup path integral
dim = 1; % dimension for integraton (1 for scalar)

% w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% %define grid where eigenfunction is well defined
% setup for eqb point at (0,0)
% bounds = 3;
% grid = -bounds:0.025:bounds; %define grid where eigenfunction is well defined
% [q1,q2] = meshgrid(grid);

q1 = Xs(:,1); q2 = Xs(:,2);
x_0 = [q1(:),q2(:)]; 
phi1_real=[]; phi2_real=[];
phi1_imag = []; phi2_imag = [];
phi1_matlab = []; phi2 = [];
options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
%options = odeset('events',@(t, x)offFrame(t, x, Dom(2)));
%options = odeset('RelTol',1e-9,'AbsTol',1e-300);
% tspan = [-20 0];
tspan = [0 20];
parfor i = 1:length(x_0)
%     waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    
    [t,x] = ode45(@(t,x)f(x(1),x(2)),tspan,x_0(i,:),options);
    % path integral for complex eigenfunctions
    % real part
    phi1_real = [phi1_real, w1_real'*x_0(i,:)' +...
        trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_real(x(:,1),x(:,2) ...
            +sin(eig_val1_imag*t)).*g1_imag(x(:,1),x(:,2))),dim)];

    % imaginary part
    phi1_imag = [phi1_imag, w1_imag'*x_0(i,:)' +...
        trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t).*g1_imag(x(:,1),x(:,2)) ...
        -sin(eig_val1_imag*t).*g1_real(x(:,1),x(:,2))),dim)];

    % compute directly
    phi1_matlab = [phi1_matlab, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x(:,1),x(:,2)), dim)];
    %phi2 = [phi2, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x(:,1),x(:,2)), dim)];
end

% F = findall(0,'type','figure','tag','TMWWaitbar');
% delete(F);

%% reshape
% phi for eqb point at (0,0)
phi1_real = reshape((phi1_real),size(q2));
phi1_imag = reshape((phi1_imag),size(q2));
phi1_matlab = reshape((phi1_matlab),size(q2));

phi1_mag  = sqrt(phi1_real.^2+phi1_imag.^2);
phi1_mag_matlab = sqrt(real(phi1_matlab).^2+imag(phi1_matlab).^2);

phi1_phase = angle(phi1_real + i*phi1_imag);
phi1_phase2 = angle(real(phi1_matlab) + i*imag(phi1_matlab));
phi1_mag_log = log(phi1_mag);

%% scatter plots
figure(1)
sz = 10; alpha = 1;
subplot(2,4,1)
Xs = [];
scatter(q1(:),q2(:),sz,phi1_mag(:),'filled','MarkerFaceAlpha',alpha); hold on;
axes1 = gca;
axis square
set(axes1,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
colorbar('southoutside')
box on
axes1.LineWidth=2;

%% helper functions

function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end