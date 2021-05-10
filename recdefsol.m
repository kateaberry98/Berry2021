%% Figure params k_r,k_d,m,n,pb0_1,pb0_2
% Code for Section 2 plots
% Data:
% 1) logisticgrowth1: 0.8,0.75,1,1,[0.4;0.4],[0.2;0.6],tfinal=200
% 2) logisticgrowth2: 0.75,0.8, ""
% 3) strongallee1: 0.6,0.2,3,4,"",tfinal=50
% 4) strongallee2: "",m=4,n=3,""
% 5) weakandstrongallee1: 0.1,0.9,3,1,"",tfinal=80
% 6) weakandstrongallee2: 0.9,0.1,""
%%
k_r = 0.9;
k_d = 0.1;
m = 3;
n = 1;

t0 = 0;
tfinal = 80;
pb0_1 = [0.4; 0.4];
pb0_2 = [0.2; 0.6];
K = sum(pb0_1);


dpbdt = @(t,pb) recdefODEs(t,pb,k_r,k_d,m,n);
[t1,pb1] = ode45(dpbdt, [t0 tfinal], pb0_1);
[t2,pb2] = ode45(dpbdt, [t0 tfinal], pb0_2);
figure(1)
plot(t1,pb1(:,1),'r','LineWidth',2)
hold on
plot(t1,pb1(:,2),'b','LineWidth',2)
hold on
plot(t2,pb2(:,1),'g','LineWidth',2)
hold on
plot(t2,pb2(:,2),'m','LineWidth',2)
xlabel('Time, t')
ylabel('p(t),b(t)')
legend('p_0 = 0.4','b_0 = 0.4','p_0 = 0.2','b_0 = 0.6','Location','Best')
hold off
ylim([0 K])
% print -dpng 

p = linspace(0,K,1000);
dpdt = @(p,K) p.*(K-p).*(k_r.*(p.^(m-1)) - k_d.*((K-p).^(n-1)));
figure(2)
plot(p,dpdt(p,K),'k','LineWidth',2)
hold on
p1=plot(pb0_1(1),dpdt(pb0_1(1),K),'r*','LineWidth',8);
hold on
p2=plot(pb0_2(1),dpdt(pb0_2(1),K),'g*','LineWidth',8);
xlabel('Panic-buyer density, p')
ylabel('dp/dt')
hold on
line(xlim,[0,0],'Color','k','LineStyle','--','LineWidth',2)
hold off
legend([p1,p2],'p_0 = 0.4','p_0 = 0.2','Location','Best')
% print -dpng 