clc; clear all; close all;

lam = -1;
f = @(t,x)(lam*x + x.^3);
z = sym('z','real');
f_sym = f(0,z);
A = subs(jacobian(f_sym,z),z,0);
fn = @(t,x)(f(t,x)-A*x);
fn_sym = fn(0,z);

[~,d,w] = eig(A);
g = @(t,x)(w*fn(t,x));
g_sym = g(0,z);

phi_an = @(x)(x./sqrt(1-x.^2));

% c1_exp = @(x)(log((1/x.^2)-1));
% st_1 = @(t,x)(1/(sqrt(exp(c1_exp(x)+2*t)+1)));
% st_2 = @(t,x)(-1/(sqrt(exp(c1_exp(x)+2*t)+1)));

dTs = [0.001, 0.01, 0.1, 1];
Tf=10;
x_inits = -0.9:0.05:0.9;
hx_all = zeros(1,size(x_inits,2));
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
figure; hold on;
for k=1:length(dTs)
    dT=dTs(k);
    tspan = 0:dT:Tf; %[0 Tf]; %
    for i=1:length(x_inits)
        y0 = x_inits(i);
        [t,y] = ode45(f, tspan, y0, opts);
        %     plot(t,y,'-o')
        %     hx=0; t=0; x=y0;
        %     for j=1:length(tspan)
        %         hx = hx + exp(t)*(y(j)^3)*dT;
        %         t = t + dT;
        %     end
        %     hx_all(i) = hx;
        tDiff = dT; %tDiff = [t(2:end); t(end)] - t;
        hx_temp = exp(lam*t).*(y.^3).*tDiff;
        hx_all(i) = sum(hx_temp);
    end
    phiPI_vals = w*x_inits + hx_all;
    plot(x_inits,phiPI_vals,'LineWidth',1.5);
end
phiAN_vals = phi_an(x_inits);
plot(x_inits,phiAN_vals,'*','LineWidth',1.5); 
legend_cell = {};
for k=1:length(dTs)
    legend_cell{k}="dT: "+string(dTs(k));
end
legend_cell{end+1} = "Analytical";
legend(legend_cell,'Location','NorthWest')