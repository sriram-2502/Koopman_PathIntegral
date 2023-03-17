clc; clear; close all;
load('matlab_matrix.mat')

%% get mag and phase
phi_mag = sqrt(Eigfunc(:,:,1).^2 + Eigfunc(:,:,2).^2);
phi_phase = angle(Eigfunc(:,:,1) + j*Eigfunc(:,:,2));

%% plot along trajectories
% plot magnitude
figure(1)
subplot(2,4,1)
Xs = [];
start_idx = 1;
alpha = 1; sz = 35;
for idx = 1:1000:1000
    plot(X(:,1,idx),X(:,2,idx),'k','LineWidth',1); hold on;
    scatter(X(:,1,idx),X(:,2,idx),sz,phi_mag(:,idx),'filled','MarkerFaceAlpha',alpha); hold on;
    %plot(X(:,1,idx),X(:,3,idx),'k','LineWidth',1); hold on;
end
colorbar
colormap jet
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$q_1$','FontSize',20, 'Interpreter','latex')
ylabel('$q_2$','FontSize',20, 'Interpreter','latex')
title('$|\psi|$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

subplot(2,4,5)
Xs = [];
start_idx = 1;
for idx = 1:1000:1000
    plot(X(:,3,idx),X(:,4,idx),'k','LineWidth',1); hold on;
    scatter(X(:,3,idx),X(:,4,idx),sz,phi_mag(:,idx),'filled','MarkerFaceAlpha',alpha); hold on;
    %plot(X(:,1,idx),X(:,3,idx),'k','LineWidth',1); hold on;
end
colorbar
colormap jet
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$\dot{q}_1$','FontSize',20, 'Interpreter','latex')
ylabel('$\dot{q}_2$','FontSize',20, 'Interpreter','latex')
title('$|\psi|$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

%% plot phase

subplot(2,4,2)
Xs = [];
start_idx = 1;
alpha = 1; sz = 35;
for idx = 1:1000:1000
    plot(X(:,1,idx),X(:,2,idx),'k','LineWidth',1); hold on;
    scatter(X(:,1,idx),X(:,2,idx),sz,phi_phase(:,idx),'filled','MarkerFaceAlpha',alpha); hold on;
    %plot(X(:,1,idx),X(:,3,idx),'k','LineWidth',1); hold on;
end
colorbar
colormap jet
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$q_1$','FontSize',20, 'Interpreter','latex')
ylabel('$q_2$','FontSize',20, 'Interpreter','latex')
title('$\psi_{phase}$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

subplot(2,4,6)
Xs = [];
start_idx = 1;
for idx = 1:1000:1000
    plot(X(:,3,idx),X(:,4,idx),'k','LineWidth',1); hold on;
    scatter(X(:,3,idx),X(:,4,idx),sz,phi_phase(:,idx),'filled','MarkerFaceAlpha',alpha); hold on;
    %plot(X(:,1,idx),X(:,3,idx),'k','LineWidth',1); hold on;
end
colorbar
colormap jet
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$\dot{q}_1$','FontSize',20, 'Interpreter','latex')
ylabel('$\dot{q}_2$','FontSize',20, 'Interpreter','latex')
title('$\psi_{phase}$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;