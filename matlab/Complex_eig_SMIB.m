%% Path integral method for eigenfunction computation for SMIB complex eigenvalue case
%% Bhagyashree version 1

clc; clear; close all;

Damp = 1.3; %damping
Pm = 5; %mechanical input
Pe = 10;
a = asin(Pm/Pe);
x = sym('x',[2;1]);
f = [x(2); -Pe*sin(x(1)+a)-Damp*x(2)+Pm ]; % syms form
A = eval(subs(jacobian(f),[x(1) x(2)],[0 0]));
fn = f - A*[x(1);x(2)];
%% Realification of eigenvalues will be needed
[V,D1,W1] = eig(A);
[W, D] = cdf2rdf(W1, D1);
w1 = W(:,1); w2 = W(:,2);
f = matlabFunction(f);
g1 = matlabFunction(w1'*fn,'vars',{x(1),x(2)});
g2 = matlabFunction(w2'*fn,'vars',{x(1),x(2)});
lam_r = real(D1(1,1)); lam_i = imag(D1(1,1));

%% set up path integral

pplim =7;
grid = -pplim:0.2:pplim; %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
U = f(q1(:)',q2(:)');
dtt= 0.01;
tspan = 0:dtt:60;
x_0 = [q1(:),q2(:)]; phi1_est=[];phi2_est=[];
figure(1)
for i = 1:length(x_0)
    [t,x] = ode45(@(t,x)f(x(1),x(2)),tspan,x_0(i,:));
%         idx = find(abs(x(:,1))>pplim);
%         if(~isempty(idx))
%             idx = idx(1);
%             t = t(1:idx); x = x(1:idx,:);
%         end
    if norm(x(end,:),2)>0.5
        phi1_est = [phi1_est, -2];
        phi2_est = [phi2_est, -2];
    else
        idx = find(sqrt(sum(x.*x,2))<1e-6);
        if(~isempty(idx))
            idx = idx(1);
            t = t(1:idx); x = x(1:idx,:);
        end
        if length(t)>1
            phi1_est = [phi1_est, w1'*x_0(i,:)' + ...
                trapz(t,exp(-lam_r*t).*cos(lam_i*t).*g1(x(:,1),x(:,2)))];
            phi2_est = [phi2_est, w2'*x_0(i,:)' - ...
                trapz(t,exp(-lam_r*t).*sin(lam_i*t).*g2(x(:,1),x(:,2)))];
        else
            phi1_est = [phi1_est, 0];
            phi2_est = [phi2_est, 0];
        end
    end
%     plot(x(:,1),x(:,2),0,0,'*')
%     title(num2str(phi1_est(end)))
%     phi1_est(end)
end

% normalize the phi and phi_est
% I_est1 = sqrt(trapz(grid,trapz(grid,reshape((phi1_est.^2),size(q2)),2)));
% % I1 = sqrt(trapz(grid,trapz(grid,Phi1(q1,q2).^2,2)));
%
% I_est2 = sqrt(trapz(grid,trapz(grid,reshape((phi2_est.^2),size(q2)),2)));
% I2 = sqrt(trapz(grid,trapz(grid,Phi2(q1,q2).^2,2)));
I_est1=1; I_est2=1; I1=1; I2=1;

%% plot

phi1_normalized = reshape((phi1_est),size(q2))/I_est1;
phi2_normalized = reshape((phi2_est),size(q2))/I_est2;
Phi_mag  = sqrt(phi1_normalized.^2+phi2_normalized.^2);
Phi_mag(Phi_mag>10)=10;
% Phi_mag(Phi_mag<1e-1) = 0;

%% Curve to compare and get the minimum value of Phi_mag
a = 1; b = 0.6; c=0;
W_t = W';
C_fit = @(z1,z2) a*(W_t(1,1)*z1+W_t(1,2)*z2).^2+b*(W_t(2,1)*z1+W_t(2,2)*z2).^2 +...
    2*c*(W_t(1,1)*z1+W_t(1,2)*z2).*(W_t(2,1)*z1+W_t(2,2)*z2);
C_val = C_fit(q1,q2);

h_sq =  @(z1,z2) z1.^2
h_sq_val =  h_sq(q1,q2);
figure (10)
p4 = surf(q1,q2,C_val,'FaceAlpha',0.5); hold on;
p5 = surf(q1,q2,h_sq_val,'FaceAlpha',0.5); hold on;
p4.EdgeColor = 'none';
p5.EdgeColor = 'none';
p5.FaceColor = 'g';
colorbar
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel(' value','interpreter','latex');
set(gca,'fontsize',20,'FontName','Times New Roman')
title('h(x) and linear part for bound','Interpreter','latex')
zlim([-3 12])
minDiff_Lin = min(min(C_val-h_sq_val))
%%
figure(11)
phiQphi = @(pr,pi)a*pr.^2+b*pi.^2;
phiQphi_val = phiQphi(phi1_normalized,phi2_normalized);
p6 = surf(q1,q2,h_sq_val,'FaceAlpha',0.5); hold on;
p7 = surf(q1,q2,phiQphi_val,'FaceAlpha',0.5); hold on;
p6.EdgeColor = 'none';
p7.EdgeColor = 'none';
p7.FaceColor = 'g';
colorbar
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel(' value','interpreter','latex');
set(gca,'fontsize',20,'FontName','Times New Roman')
title('h(x) and $\Phi$ for bound','Interpreter','latex')
zlim([-3 12])
minDiff_phi = min(min(phiQphi_val-h_sq_val))
%%

figure(8)
p4 = surf(q1,q2,Phi_mag,'FaceAlpha',0.5); hold on;
p5 = surf(q1,q2,C_val,'FaceAlpha',0.5); hold on;
p4.EdgeColor = 'none';
p5.EdgeColor = 'none';
p5.FaceColor = 'g';
colorbar
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('Phi magnitude value','interpreter','latex');
set(gca,'fontsize',20,'FontName','Times New Roman')
title('Value Magnitude of $\Phi(x)$, t= 30','Interpreter','latex')
zlim([-3 12])
%%
figure(9)
Phi_val = 50*phi1_normalized.^2- q1.^2;
p9 = surf(q1,q2,Phi_val,'FaceAlpha',0.5); hold on;
p10 = surf(q1,q2,zeros(size(Phi_val)),'FaceAlpha',0.5); hold on;
p10.EdgeColor = 'none';
p10.FaceColor = 'g';
p9.EdgeColor = 'none';
 zlim([-50 200])
%%
Phi_zero = ones(size(Phi_mag));
Phi_zero(abs(Phi_mag)<1e-1) = 0;
figure(6)
p6 = surf(q1,q2,phi1_normalized,'FaceAlpha',0.5); hold on;
p6.EdgeColor = 'none';
figure(7)
p7 = surf(q1,q2,phi2_normalized,'FaceAlpha',0.5); hold on;
p7.EdgeColor = 'none';

figure(5)
p1 = pcolor(q1,q2,Phi_mag); hold on;
colormap jet
colorbar
set(p1,'Edgecolor','none')
title('Magnitude of $\Phi(x)$','Interpreter','latex')
% title(sprintf(' V(x),  $T_{cr}= %g$ sec', Tc1),'Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20,'FontName','Times New Roman')
% phase portrait and vectror field


options = odeset('RelTol',1e-6,'AbsTol',1e-30);

xl = -pplim;
xh = pplim;
yl = -pplim;
yh = pplim;
pts1=1.5; % number of initial condictions
axX = max(abs(xl),abs(xh));
axY = max(abs(xl),abs(xh));
pts =1; % number of arrows in pahse portrait
figure(5);
hold on
f = @(t, x) [x(2,:); -Pe*sin(x(1,:)+a)-Damp*x(2,:)+Pm ]; % Stable

for x0 = xl:pts1:xh
    for y0 = yl:pts1:yh
        [ts,xs] = ode45(@(t,x)f(t, x),tspan,[x0 y0], options);
        plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    end
end
axis([-pplim pplim -pplim pplim])
