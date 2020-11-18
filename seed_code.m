%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was used to generate Figures 4-9 in "Red light and the
% dormancy-germination decision in Arabidopsis seeds" (in revision @ The
% Bulletin of Mathematical Biology)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter values
p1=1; p2=50; p3=1; p4=1; p5=0.1; p6=10; p7=1; p8=0.1; p9=10; p10=1; p11=1;
ep=0.01; p12=1000;
% Figure 4
% Slow DoG model (equations 9-12)
% x(1) = b, x(2) = p, x(3) = r, x(4) = c
bprc_slow = @(t,x) [p1*(p2*sech(ep*(t-p12))/(1+p3*x(2)*x(2))-x(1));...
    p4*(p5+p6*x(3)/(1+x(3))-(x(2)+x(1)*x(2))-p7*(x(2)*x(3)-x(4)));...
    p8+p9*x(2)/(1+x(2))-x(3)-p10*(x(2)*x(3)-x(4));...
    p11*(x(2)*x(3)-x(4));];
% Solve system
[t2,x2] = ode45(bprc_slow,[0 2000],[1 1 1 1]);
% Plotting
figure(4)
plot(t2,p2*sech(ep*(t2-p12)),'LineWidth',3,'Color','b')
hold on
plot(t2,x2(:,4),'LineWidth',3,'Color','r')
xlabel('Time')
axis([0 2000 0 100])
legend('s(\epsilon,p_2,p_{12},\tau)','[PIF1-RVE1]')
set(gca,'Fontsize',32);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5
% Parameter values
p1=1; p2=25; p3=1; p4=1; p5=0.1; p6=10; p7=1; p8=0.1; p9=10; p10=1; p11=1;
% Vary p13 values
j=[0.005,0.01,0.015, 0.02, 0.05, 0.1,0.15];
for k=1:length(j)
    % Wigwag DoG model (equations 13-16)
    % x(1) = b, x(2) = p, x(3) = r, x(4) = c
    bprc_var = @(t,x) [p1*((p2*sin(j(k)*t)^2)/(1+p3*x(2)*x(2))-x(1));...
        p4*(p5+p6*x(3)/(1+x(3))-(x(2)+x(1)*x(2))-p7*(x(2)*x(3)-x(4)));...
        p8+p9*x(2)/(1+x(2))-x(3)-p10*(x(2)*x(3)-x(4));...
        p11*(x(2)*x(3)-x(4));];
    % Solve system
    [t3,x3] = ode45(bprc_var,[0 10000],[1 1 1 1]);
    % Plotting
    txt=['p_{13} =' num2str(j(k))];
    figure(5)
    % Avoid initial transients
    X=x3(5000:end,1);
    Y=x3(5000:end,4);
    plot(X,Y,'LineWidth',3,'DisplayName',txt)
    hold on
    set(gca,'Fontsize',32);
    xlabel('[PhyB*]');
    ylabel('[PIF1-RVE1]');
    legend show
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 6-8
% Parameter Values
p1=1; p3=1; p4=1; p5=0.1; p6=10; p7=1; p8=0.1; p9=10; p10=1; p11=1;
% Vary p2
jj=1:0.2:50;
% Vary p13
j=0.005:0.005:0.5;
% Arrays for storing minima and maxima
minY=zeros(length(jj),length(j));
maxY=zeros(length(jj),length(j));
for b=1:length(jj)
    for k=1:length(j)
        % Wigwag DoG model (equations 13-16)
        % x(1) = b, x(2) = p, x(3) = r, x(4) = c
        bprc_var = @(t,x) [p1*((jj(b)*sin(j(k)*t)^2)/(1+p3*x(2)*x(2))-x(1));...
            p4*(p5+p6*x(3)/(1+x(3))-(x(2)+x(1)*x(2))-p7*(x(2)*x(3)-x(4)));...
            p8+p9*x(2)/(1+x(2))-x(3)-p10*(x(2)*x(3)-x(4));...
            p11*(x(2)*x(3)-x(4));];
        % Solve the system
        [t3,x3] = ode45(bprc_var,[0 10000],[1 1 1 1]);
        % Avoid initial transients
        Y=x3(5000:end,4);
        % Store minimum
        minY(b,k)=min(Y);
        % Store maximum
        maxY(b,k)=max(Y);
    end
end
% Plotting
[X,Y]=meshgrid(j,jj);
figure(6)
contourf(X,Y,minY)
set(gca,'Fontsize',32);
xlabel('p_{13}')
ylabel('p_{2}')
colorbar

figure(7)
contourf(X,Y,maxY)
set(gca,'Fontsize',32);
xlabel('p_{13}')
ylabel('p_{2}')
colorbar

figure(8)
contourf(X,Y,maxY-minY)
set(gca,'Fontsize',32);
xlabel('p_{13}')
ylabel('p_{2}')
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 9
% Parameter Values
p1=1; p3=1; p4=1; p5=0.1; p6=10; p7=1; p8=0.1; p9=10; p10=1; p11=1; p13=0.01;
% p2
v=5:5:25;
for b=1:length(v)
    % Wigwag DoG model (equations 13-16)
    % x(1) = b, x(2) = p, x(3) = r, x(4) = c
    bprc_var = @(t,x) [p1*((v(b)*sin(p13*t)^2)/(1+p3*x(2)*x(2))-x(1));...
        p4*(p5+p6*x(3)/(1+x(3))-(x(2)+x(1)*x(2))-p7*(x(2)*x(3)-x(4)));...
        p8+p9*x(2)/(1+x(2))-x(3)-p10*(x(2)*x(3)-x(4));...
        p11*(x(2)*x(3)-x(4));];
    % Solve the system
    [t3,x3] = ode45(bprc_var,[0 10000],[1 1 1 1]);
    % Plotting
    txt=['p_{2} =' num2str(v(b))];
    % Avoid initial transients
    Y=x3(5000:end,4);
    figure(9)
    plot(t3,x3(:,4),'LineWidth',2,'DisplayName',txt)
    legend('Location','northeastoutside')
    set(gca,'Fontsize',32);
    xlabel('time');
    ylabel('[PIF1-RVE1]');
    hold on
end

