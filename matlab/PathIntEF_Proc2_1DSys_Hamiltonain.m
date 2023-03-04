clear; close all; clc;
% set(0,'DefaultLineLineWidth',2.5)
% set(0,'defaultfigurecolor',[1 1 1])

%% Construct the Hamiltonian System 
lam=1; beta=2;
f = @(t,x,u)(-lam.*x + beta*x.^3 + u);
h = @(t,x,u)(x);

n=1; % dimension
x = sym('x',[n;1]);
u = sym('u',[1;1]);
x0 = 0;

A = double(subs(jacobian(f(0,x,0),x),x,x0)); 
B = double(subs(jacobian(f(0,x0,u),u),u,0));
%linear eigenfunctions realization from complex eigenfunctions
% [~,D1,W1] = eig(A);
% [W, D] = cdf2rdf(W1, D1);

x = sym('x',[n;1],'real');
u = sym('u',[1;1],'real');
x0 = 0;

A0 = double(subs(jacobian(f(0,x,0),x),x,x0)); 
B0 = double(subs(jacobian(f(0,x0,u),u),u,0));

z = sym('z',[2*n;1],'real'); 

f_sym = f(0,z(1:n),0);
g = jacobian(f(0,z(1:n),u),u);
Df  = jacobian(f_sym, z(1:n));
R_x = g*g'; %R(x)
R_0 = B0*B0'; %R(0)

% q_x = 1/2*(xp(1)^2+xp(2)^2+xp(3)^2+xp(4)^2); %cost function with Q=1
q_x = h(0,z(1:n),0).^2;
Dq  = jacobian(q_x, z(1:n));
DDq = hessian(q_x, z(1:n));

p = z(n+1:end);
p1 = jacobian(p'*R_x*p, z(1:n));

H_sym =  [f_sym-R_x*z(n+1:end);  -Df'*z(n+1:end)+0.5.*p1'-Dq']; 
H_sys = matlabFunction(H_sym, 'Vars',{z(1:n), z(n+1:2*n)}); %equation 75

H = [A0, -R_0; -double(DDq), -A0']; %equation 76

%% Path Integral setup
Hn = H_sym - H*z;
[V,D1,W1] = eig(H);
[W, D] = cdf2rdf(W1, D1);
w1 = W(:,1); w2 = W(:,2);

g1 = matlabFunction(w1'*Hn,'vars',{z(1),z(2)});
g2 = matlabFunction(w2'*Hn,'vars',{z(1),z(2)});
l1 = D(1,1); l2 = D(2,2);

pplim = 1;
grid = -pplim:0.01:pplim; %define grid where eigenfunction is well defined
[q1,q2] = meshgrid(grid);
U = H_sys(q1(:)',q2(:)');
dtt= 0.01;
tspan = 0:dtt:10;
x_0 = [q1(:),q2(:)]; 
% ll=-0.0012:0.0002:0.0012;
% x_1=(V(:,1)*ll)';
phi1_est=[];phi2_est=[];phi2_nl=[];
% figure(1)
for i = 1:length(x_0)
    [t,x] = ode45(@(t,x)H_sys(x(1),x(2)),tspan,x_0(i,:));
%         idx = find(abs(x(:,1))>pplim);
%         if(~isempty(idx))
%             idx = idx(1);
%             t = t(1:idx); x = x(1:idx,:);
%         end
    if 0%norm(x(end,:),2)>10
%         phi1_est = [phi1_est, -2];
        phi2_est = [phi2_est, -2];
    else
%         idx = find(sqrt(sum(x.*x,2))>1e-6);
        idx = find((abs(x(:,1))>1.75*pplim)|(abs(x(:,2))>1.75*pplim));
        if(~isempty(idx))
            idx = idx(1);
            t = t(1:idx); x = x(1:idx,:);
        end
        if length(t)>1
%             phi1_est = [phi1_est, w1'*x_0(i,:)' + ...
%                 trapz(t,exp(-l1*t).*g1(x(:,1),x(:,2)))];
            ef2_nlp = trapz(t,exp(-l2*t).*g2(x(:,1),x(:,2)));
            phi2_est = [phi2_est, w2'*x_0(i,:)' + ef2_nlp];
            phi2_nl = [phi2_nl, ef2_nlp];
        else
%             phi1_est = [phi1_est, 0];
            phi2_est = [phi2_est, 0];
            phi2_nl = [phi2_nl, 0];
        end
    end
%     hold on
%     plot(x(:,1),x(:,2),'--',x(1,1),x(1,2),'k+',x(end,1),x(end,2),'bx',0,0,'r*')
%     title(num2str(phi1_est(end)))
%     phi1_est(end)
end
% plotv(V); plotv(-V);
I_est1=1; I_est2=1; I1=1; I2=1;

%% plot
% hold off
% figure(2);
% phi1_normalized = reshape((phi1_est),size(q2))/I_est1;
% phi1_mod = phi1_normalized;
% phi1_mod((phi1_mod)>1000)=1000;
% phi1_mod((phi1_mod)<-1000)=-1000;
% surf(q1,q2,phi1_mod); hold on;
% surf(q1,q2,Phi1(q1,q2)/I2) % phi2 corrsponds to eig -5 which is phi1_norm
% axes1 = gca;
% set(axes1,'FontSize',15);
% xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
% ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
% title('Stable Eigenfunction'); colorbar

figure(3);
phi2_normalized = reshape((phi2_est),size(q2))/I_est2;
phi2_mod = phi2_normalized;
phi2nl_mod = reshape((phi2_nl),size(q2))/I_est2;
% phi2_mod((phi2_mod)>1)=1;
% phi2_mod((phi2_mod)<-1)=-1;
p3 = surf(q1,q2,phi2_mod); hold on;
p3.EdgeColor = 'none';
% surf(q1,q2, Phi2(q1,q2)/I1)
axes3 = gca;
set(axes3,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title('Unstable Eigenfunction'); colorbar

figure(4)
p4 = surf(q1,q2,phi2nl_mod); hold on;
p4.EdgeColor = 'none';
axes4 = gca;
set(axes4,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title('Unstable Eigenfunction: Nonlinear Part'); colorbar

figure(5)
zls = phi2_mod;
zls(abs(phi2_mod)>4e-3)=-1;
p5 = pcolor(q1,q2,zls); %colorbar;
p5.EdgeColor = 'none';
title("Zero-level set: Unstable eigenfunction :: "+ ...
    "$\dot{x}=-x+$"+num2str(beta)+"$x^3+u$", 'FontSize',11, ...
    "Interpreter","latex")
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
hold on

%% Stable Lagrangian submanifold
pres = 0.01; x0 = -pplim:pres:pplim;
f_an = @(x)(f(0,x,0)); g_an = @(x)1;
R_an = @(x)g_an(x)*g_an(x)'; q_an = @(x)x.^2;
hje = @(x,p)p*f_an(x) - 0.5*R_an(x).*p.^2 + q_an(x);
Psol_an = solve(hje(z(1),z(2))==0,z(2)); % x==z(1); p==z(2)
Pfun_an = matlabFunction(Psol_an);
pVals_an = Pfun_an(x0);
Plt1 = plot(x0,pVals_an(1,:),'r--','LineWidth',2.5);
legend(Plt1,'Analytical','Location','northwest')
