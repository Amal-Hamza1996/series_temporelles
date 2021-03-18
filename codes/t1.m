D = 100:100:1000;
SDR_5000 = [9.7773, 9.7952, 9.7817, 9.7984, 9.7818, 9.7476, 9.7642, 9.7450, 9.7618, 9.7721];
SDR_10000 = [9.1408, 9.1858, 9.1811, 9.1834, 9.1657, 9.1822, 9.1844, 9.1791, 9.1755, 9.1772];
%% Les figures
figure(1)
plot(D,SDR_5000,'k-','Linewidth',2);
hold on;
plot(D,SDR_10000,'r-','Linewidth',2);
title('SDR en fonction de horizon de prédiction pour 5000, et 10000 utilisateurs.');
legend('5000 utilisateurs','10000 utilisateurs');
ylabel('SDR');
xlabel('Horizon de prédiction en (ms)');
grid on;