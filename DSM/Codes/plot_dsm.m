%plot
pol_iter_7 = pol_iter;
pol_tariff_7 = pol_tariff;
pol_total_tariff_7 = pol_total_tariff;
policy_evol_7 = policy_evol;

pol_total_tariff_alp0p5_eps0p2=pol_total_tariff_1;
pol_total_tariff_alp0p5_eps0p4=pol_total_tariff_5;
pol_total_tariff_alp0p5_eps0p0=pol_total_tariff_7;

%q vs sarsa
x=1:30;
y1=pol_total_tariff_alp0p5_eps0p0(:,1:30);
y2=pol_total_tariff_alp0p5_eps0p2(:,1:30);
y3=pol_total_tariff_alp0p5_eps0p4(:,1:30);

figure
plot(x,y1/100,'k',x,y2/100,'--k',x,y3/100,':k','LineWidth',2);
title('Cost Function Gradient - SARSA vs Q-learning')
xlabel('Iteration Number - Policy Improvement') % x-axis label
ylabel('Estimated Tariff') % y-axis label
legend('Q-learning','SARSA Epsilon = 0.2', 'SARSA Epsilon = 0.4');

%q schedule- plot double y
figure
stackedbar = @(x, A) bar(x, A, 0.75, 'stack');
prettyline = @(x, y) plot(x, y, 'k', 'LineWidth', 1);
[ax, h1, h2] = plotyy(1:24, [q_app1' q_app2' q_app3' q_app4' q_app5' q_app6'], 1:24, price_per_hour, stackedbar, prettyline);
set(ax, 'XLim', [0 24]);
set(ax, 'XTick', 1:24);
set(ax(2), 'YLim', [0 70]);
set(ax(2), 'YTick', 0:10:70);
colormap gray;
title('Q-learning - Appliance Schedule at 30 iterations');
xlabel('Hour of Day');
ylabel(ax(1),'Estimated Tariff (in Cents)');
ylabel(ax(2),'Price per Hour (in Cents)');
legend('Appliance 1', 'Appliance 2', 'Appliance 3', 'Appliance 4', 'Appliance 5', 'Appliance 6', 'Tariff per KWh: 1st Jan 2016', 'Location', 'southeast');


%sarsa schedule- plot double y
figure
stackedbar = @(x, A) bar(x, A, 0.75, 'stack');
prettyline = @(x, y) plot(x, y, 'k', 'LineWidth', 1);
[ax, h1, h2] = plotyy(1:24, [sarsa_app1' sarsa_app2' sarsa_app3' sarsa_app4' sarsa_app5' sarsa_app6'], 1:24, price_per_hour, stackedbar, prettyline);
set(ax, 'XLim', [0 24]);
set(ax, 'XTick', 1:24);
set(ax(2), 'YLim', [0 70]);
set(ax(2), 'YTick', 0:10:70);
colormap gray;
title('SARSA - Appliance Schedule at 30 iterations');
xlabel('Hour of Day');
ylabel(ax(1),'Estimated Tariff (in Cents)');
ylabel(ax(2),'Price per Hour (in Cents)');
legend('Appliance 1', 'Appliance 2', 'Appliance 3', 'Appliance 4', 'Appliance 5', 'Appliance 6', 'Tariff per KWh: 1st Jan 2016', 'Location', 'southeast');

