% Path integral Eigenfunctions
% Van der pol Example
clear; close all; clc;
set(0,'DefaultLineLineWidth',2.5) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

mu=1; alpha = 1;
%f = @(t, x)0.5*[x(2,:);  mu*x(2,:)-x(1,:)-mu*x(1,:).^2.*x(2,:) ];
%f = @(t, x) [x(2,:); -2*x(1,:)+(1/3)*x(1,:).^3-x(2,:)  ];
f = @(t, x) [alpha*x(2); alpha*(mu*x(2) - x(1) - mu*x(1)^2*x(2))]; 

%linearization
x = sym('x',[2;1], 'real');
A = double(subs(jacobian(f(0,x),x),x,[0;0]));
fn = @(x) f(0,x)-A*x;

%spectrum
[~,D,W] =eig(A);
[Wr,  Dr] = cdf2rdf(W ,D);
Dtemp = diag(D);
eig_val1_real = real(Dtemp(1)); eig_val1_imag = imag(Dtemp(1)); 
eig_val2_real = real(Dtemp(2)); eig_val2_imag = imag(Dtemp(2)); 
w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
w2_real = real(W(:,2)); w2_imag = imag(W(:,2));

l1 = D(1,1); l2 = D(2,2);
w1 = W(:,1); w2 = W(:,2);
g1 = @(x) w1'*fn(x);
g2 = @(x) w2'*fn(x);

%% Set up path Integral
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

Dom = [-4 4]; ds = 0.05;
grid = Dom(1):ds:Dom(2); %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
QQ = [q1(:)';q2(:)'];
U = f(0,QQ);

%quiver(q1(:),q2(:),U(1,:)',U(2,:)')
x_0 = [q1(:),q2(:)]; phi1_est=[];phi2_est=[];
options = odeset('RelTol',1e-9,'AbsTol',1e-300,'events',@(t, x)offFrame(t, x, Dom(2)));
%options = odeset('events',@(t, x)offFrame(t, x, Dom(2)));
for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))

    [t,x] = ode45(@(t,x)f(t,x),[0 3],x_0(i,:), options);
    phi1_est = [phi1_est, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x')')];
    phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% Plots
phi1_est_real = 2*real(phi1_est);
phi1_est_imagin = -2*imag(phi1_est);

% phi1_phase = atan(phi1_est_imagin./phi1_est_real);
% phi1_mag = sqrt(phi1_est_real.^2+phi1_est_imagin.^2);

phi1_phase = atan(imag(phi1_est)./real(phi1_est));
phi1_mag = sqrt(real(phi1_est).^2+imag(phi1_est).^2);


figure(1)
subplot(2,2,1)
phi1_estp = reshape((phi1_phase),size(q2));
surf(q1,q2,phi1_estp);  hold on;
title('$\phi_{phase}$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$p$','interpreter','latex');
set(gca,'fontsize',20)

subplot(2,2,2)
phi2_estp = reshape((phi1_mag),size(q2));
surf(q1,q2,phi2_estp); hold on;
title('$\phi_{mag}$','Interpreter','latex')
xlabel('$x$','interpreter','latex');
ylabel('$p$','interpreter','latex');
set(gca,'fontsize',20)

subplot(2,2,3)
zls = phi1_estp;
%zls(abs(phi1_estp)>5e-2)=-1; %zero-level set
p1=pcolor(q1,q2,zls); 
colorbar
colormap parula
set(p1,'Edgecolor','none')

subplot(2,2,4)
zls1 = phi2_estp;
%zls1(abs(phi2_estp)>5e-2)=-1; %zero-level set
p2=pcolor(q1,q2,zls1); 
colorbar
colormap parula
set(p2,'Edgecolor','none')

% figure
% p2=pcolor(q1,q2,zls1+zls); 
%     colorbar
% colormap parula
% set(p2,'Edgecolor','none')

%% plots for paper
ic_pts = 5;
tspan = [0 10];
figure(2)
subplot(1,2,1)
phi1_mag = reshape((phi1_mag),size(q2));
p1 = pcolor(q1,q2,phi1_mag); hold on;
set(p1,'Edgecolor','none')
colormap jet

f = @(t,x)[x(2); mu*x(2) - x(1) - mu*x(1)^2*x(2)]; 
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end

axes1 = gca;
axis square
axis([Dom(1) Dom(2) Dom(1) Dom(2)])
set(axes1,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
%title ('Unstable Eigenfunction $|\psi_1(x)|$ at (0,0)','FontSize',20, 'Interpreter','latex')
colorbar

subplot(1,2,2)
phi1_phase = reshape((phi1_phase),size(q2));
p2 = pcolor(q1,q2,phi1_phase); hold on;
set(p2,'Edgecolor','none')
colormap jet

f = @(t,x)[x(2); mu*x(2) - x(1) - mu*x(1)^2*x(2)]; 
xl = Dom(1); xh = Dom(2);
yl = Dom(1); yh = Dom(2);
for x0 = linspace(Dom(1), Dom(2), ic_pts)
    for y0 = linspace(Dom(1), Dom(2), ic_pts)
        [ts,xs] = ode45(@(t,x)f(t,x),tspan,[x0 y0]);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end

axes2 = gca;
axis square
axis([Dom(1) Dom(2) Dom(1) Dom(2)])
set(axes2,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
%title ('Unstable Eigenfunction Phase $\psi_1(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
colorbar

%% helper functions

function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (max(abs(Y))>4.*Dom) | (min(sum(abs(Y)))<1e-2);
isterminal=1;
direction=0;
end
