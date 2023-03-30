clc; clear; %close all;
load('robotic_arm_unstable.mat')

%% get real value
phi1 = Eigfunc;
q = X;

%% wrap angles
for dt = 1:1000
    for i = 1:200
        if(mod(q(dt,1,i),2*pi) > 0)
            %disp('mod 1')
            q(dt,1,i) = mod(q(dt,1,i),2*pi);
        end
        if(mod(q(dt,2,i),2*pi) > 0)
            %disp('mod 2')
            q(dt,2,i) = mod(q(dt,2,i),2*pi); 
        end
        if(q(dt,1,i) < 0)
            disp('mod 3')
            q(dt,1,i) = mod(-q(dt,1,i),2*pi);
        end
        if(q(dt,2,i) < 0)
            disp('mod 4')
            q(dt,2,i) = mod(-q(dt,2,i),2*pi);
        end
    end
end

%% plot along trajectories
% plot on q1 q2
figure(1)
subplot(2,4,5)
Xs = [];
start_idx = 1;
alpha = 1; sz = 35;
for idx = 1:200
    %plot(q(:,1,idx),q(:,2,idx),'k','LineWidth',1); hold on;
    scatter(q(:,1,idx),q(:,2,idx),sz,phi1(:,idx),'filled','MarkerFaceAlpha',alpha); hold on;
    %plot(X(:,1,idx),X(:,3,idx),'k','LineWidth',1); hold on;

    % plot start
    %scatter(q(1,1,idx),q(1,2,idx),100,'x','MarkerEdgeColor','green'); hold on;
    % plot end
    %scatter(q(end,1,idx),q(end,2,idx),100,'square','MarkerEdgeColor','black'); hold on;
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

subplot(2,4,6)
Xs = [];
start_idx = 1;
for idx = 1:200
    %plot(q(:,3,idx),q(:,4,idx),'k','LineWidth',1); hold on;
    scatter(q(:,3,idx),q(:,4,idx),sz,phi1(:,idx),'filled','MarkerFaceAlpha',alpha); hold on;
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