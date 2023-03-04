%% eigenfunctions for vanderpol system
clc; clear; close all;
%% system description
% nonlinear ode x_dot = f(x)
x = sym('x',[3;1]); 
% sigma = 10; rho= 28; beta = 8/3;
sigma = 10; rho= 28; beta = 8/3;
f = [sigma*(x(2)-x(1)); x(1)*(rho-x(3)) - x(2); x(1)*x(2) - beta*x(3)]; 

%% linearization at (0,0,0) saddle
eqb_point = [0 0 0];
A = eval(subs(jacobian(f),[x(1) x(2) x(3)],eqb_point));
[V,D,W] = eig(A);
D = diag(D);

eig_val1_real = real(D(1)); eig_val1_imag = imag(D(1)); 
eig_val2_real = real(D(2)); eig_val2_imag = imag(D(2));
eig_val3_real = real(D(3)); eig_val3_imag = imag(D(3)); 
w1_real = real(W(:,1)); w1_imag = imag(W(:,1));
w2_real = real(W(:,2)); w2_imag = imag(W(:,2));
w3_real = real(W(:,3)); w3_imag = imag(W(:,3));

% define nonlinear part x_dot = Ax + fn(x)
fn = f - A*[x(1);x(2);x(3)];

% define matlab functions
g1_real = matlabFunction(w1_real'*fn,'vars',{x(1),x(2),x(3)});
g1_imag = matlabFunction(w1_imag'*fn,'vars',{x(1),x(2),x(3)});
g2_real = matlabFunction(w2_real'*fn,'vars',{x(1),x(2),x(3)});
g2_imag = matlabFunction(w2_imag'*fn,'vars',{x(1),x(2),x(3)});
g3_real = matlabFunction(w3_real'*fn,'vars',{x(1),x(2),x(3)});
g3_imag = matlabFunction(w3_imag'*fn,'vars',{x(1),x(2),x(3)});

f = matlabFunction(f);
%% setup path integral
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

dim = 1; % dimension for integraton (1 for scalar)

% %define grid where eigenfunction is well defined
% setup for eqb point at (0,0)
bounds = 2;
grid = -bounds:0.25:bounds; %define grid where eigenfunction is well defined
[q1,q2,q3] = meshgrid(grid);
U = f(q1(:)',q2(:)', q3(:)');
x_0 = [q1(:),q2(:),q3(:)]; 
phi1_real=[]; phi2_real=[]; phi3_real=[];

%
Xs = [];
for i = 1:length(x_0)
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    tspan = [0 10];
    [t,x] = ode45(@(t,x)f(x(1),x(2),x(3)),tspan,x_0(i,:));
    Xs = [Xs;x];
    
    % check if sol stays within the region of interest
    % if sol leaves this circle, dont compute path integral
    if norm(x(end,:)-eqb_point)>100
        phi1_real = [phi1_real, -2];
        phi2_real = [phi2_real, -2];
        phi3_real = [phi3_real, -2];

    else
        % terminate t and x when they enter linear region
        idx = find(norm(x(end,:)-eqb_point)<1e-6);
        if(~isempty(idx))
            idx = idx(1);
            t = t(1:idx); x = x(1:idx,:);
        end

        % path integral for complex eigenfunctions
        % real part
        phi1_real = [phi1_real, w1_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val1_real*t).*(cos(eig_val1_imag*t)+sin(eig_val1_imag*t)).*g1_real(x(:,1),x(:,2),x(:,3)),dim)];
    
        phi2_real = [phi2_real, w2_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val2_real*t).*(cos(eig_val2_imag*t)-sin(eig_val2_imag*t)).*g2_real(x(:,1),x(:,2),x(:,3)),dim)];

        phi3_real = [phi3_real, w3_real'*x_0(i,:)' +...
            trapz(t,exp(-eig_val3_real*t).*(cos(eig_val3_imag*t)-sin(eig_val3_imag*t)).*g3_real(x(:,1),x(:,2),x(:,3)),dim)];
    end
end
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% reshape to plot
% phi for eqb point at (0,0,0)
phi1_real = reshape((phi1_real),size(q2));
phi2_real = reshape((phi2_real),size(q2));
phi3_real = reshape((phi3_real),size(q2));

%% plot
ic_pts = 10;

tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for i=1:4 
    nexttile(1)
     for x0 = grid
        for y0 = grid
            for z0 = grid
                ic = [x_0(x0,1);x_0(y0,2);x_0(z0,3)]; 
                if(norm(ic-[x0;y0;z0])<=1e-2)
                    scatter3(x0,y0,z0,10,phi1_real(x0,y0,z0))
                end
            end
        end
    end

    axes1 = gca;
    axis square
    axis([-bounds bounds -bounds bounds])
    set(axes1,'FontSize',15);
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    %title ('Unstable Eigenfunction $|\psi_1(x)|$ at (0,0)','FontSize',20, 'Interpreter','latex')
    colorbar
    
    nexttile(2) 
    for x0 = grid
        for y0 = grid
            for z0 = grid
                ic = [x_0(x0,1);x_0(y0,2);x_0(z0,3)]; 
                if(norm(ic-[x0;y0;z0])<=1e-2)
                    scatter3(x0,y0,z0,10,phi2_real(x0,y0,z0))
                end
            end
        end
    end
    
    axes2 = gca;
    axis square
    axis([-bounds bounds -bounds bounds])
    set(axes2,'FontSize',15);
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    %title ('Unstable Eigenfunction Phase $\psi_1(x)$ at (0,0)','FontSize',20, 'Interpreter','latex')
    colorbar
end