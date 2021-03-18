%% Les valeurs de SDR
%D = 1:10; % distance de prédiction
D = 100:100:1000;
%D_ar_nn = 100:100:600;
D_ar_nn = 10:10:60;
SDR_5000 = [0.0861, 0.0857, 0.0857, 0.0861, 0.0868, 0.0900, 0.0718, 0.0674, 0.0545, 0.0575]; % 5000 utilisateurs.
SDR_10000 = [-0.0047, -0.0050, -0.0054, -0.0053, -0.0042, -0.0042, -0.0078, -0.0072, 0.0023, -0.0012]; % 10000 utilisateurs.
SDR_10000_ar_nn = [-0.0058, -0.1547, -0.1472, -0.0991, -0.1049, -0.2214];

%% Agrégats de flux égale à 10s.
D_ar_nn_10 = 10:10:60;
SDR_500_10_s = [1.4490, -1.3300, -0.6346, -0.9639, -0.2251, -0.1863];
SDR_750_10_s = [-0.3399, -1.1526, -0.4016, -0.6246, -0.2836, -0.4905 ];
SDR_5000_10_s = [2.8091,-0.4105, -0.3568, -0.7780, -0.4448, -1.1945 ];
SDR_10000_10_s = [1.6298, -0.3204, -0.5629, -1.3501, -1.7768, -0.8659  ];

%% nombre d'utilisateur variable
D_ar_1s = 1:1:10;
D_ar_10s = 10:10:120;
D_ar_1s = [D_ar_1s,20,30,40,50,60];
D_ar_10s = [D_ar_10s,180,240,300,360,420,480,540,600,660,720,780,840,900,960];
SDR_ar_10s = [3.1009,2.2555,1.9923,1.9673,2.1978,2.0177,2.0596,2.0693,2.0750,2.2355, 2.1320, 1.9824, 2.0994,1.9948,1.7412,1.7004,1.5412,1.3735,1.3568,1.1337,1.1822,0.9727,0.7826,0.8326,0.5605,0.1909,];

%% nombre variable d'uitilisateur pour un flux égale à 2 min.
D_ar_variable = 2:2:14;
SDR_scenario_1_entrant = [8.0392,7.0262,5.4526,3.9520,2.4751,1.7484,0.8819]; 
SDR_scenario_1_sortant = [8.0919,7.0546,5.4621,3.9618,2.4798,1.7492,0.8836];
SDR_scenario_2_entrant = [];
SDR_scenario_2_entrant = [];
SDR_scenario_3_entrant = [5.7948,5.3450,3.8000,3.6009,3.5102,1.6184,0.7091];
SDR_scenario_3_sortant = [5.8099,5.3444,3.9091,3.6130,3.5285,2.5298,0.7904];
SDR_scenario_4_entrant = [1.6125, 1.5903,1.6239,0.7172,0.4964,0.6067,0.4591];
SDR_scenario_4_sortant = [1.5965,1.5962,1.6032,0.7100,0.4945,0.5974,0.4527];
SDR_scenario_5_entrant = [];
SDR_scenario_5_sortant = [];
SDR_scenario_6_entrant = [];
SDR_scenario_6_sortant = [];

%% les figures.
% figure(1)
% plot(D,SDR_5000,'k-','Linewidth',2);
% hold on;
% plot(D,SDR_10000,'r-','Linewidth',2);
% title('SDR en fonction de horizon de prédiction pour 5000, et 10000 utilisateurs.');
% legend('5000 utilisateurs','10000 utilisateurs');
% ylabel('SDR');
% xlabel('Horizon de prédiction en (ms)');
% grid on;

% figure(2)
% plot(D_ar_nn,SDR_500_10_s,'g-','Linewidth',2);
% hold on;
% plot(D_ar_nn,SDR_750_10_s,'r-','Linewidth',2);
% hold on;
% plot(D_ar_nn,SDR_5000_10_s,'k-','Linewidth',2);
% hold on;
% plot(D_ar_nn,SDR_10000_10_s,'y-','Linewidth',2);
% title('SDR en fonction de horizon de prédiction');
% legend('500 utilisateurs','750 utilisateurs','5000 utilisateurs','10000 utilisateurs');
% ylabel('SDR');
% xlabel('Horizon de prédiction en seconde (s)');
% grid on;
% 
% figure(3)
% plot(D_ar_nn,SDR_10000_ar_nn,'g-','Linewidth',2);
% title('SDR en fonction de horizon de prédiction pour une granularité de 10s');
% ylabel('SDR');
% xlabel('Horizon de prédiction en seconde (s)');
% grid on;
% 
% figure(4)
% plot(D_ar_1s,SDR_ar_1s,'c-','Linewidth',2);
% title('SDR en fonction de horizon de prédiction pour une granularité de 1s');
% ylabel('SDR');
% xlabel('Horizon de prédiction en seconde (s)');
% grid on;

% figure(5)
% plot(D_ar_10s,SDR_ar_10s,'c-','Linewidth',2);
% title('SDR en fonction de horizon de prédiction pour une granularité de 10s');
% ylabel('SDR');
% xlabel('Horizon de prédiction en seconde (s)');
% grid on;

figure(6)
plot(D_ar_variable,SDR_scenario_4_entrant,'g-','Linewidth',2);
hold on;
plot(D_ar_variable,SDR_scenario_4_sortant,'r-','Linewidth',2);
title('SDR en fonction de horizon de prédiction pour une granularité de 2 min');
legend('Débit entrant','Débit sortant');
ylabel('SDR');
xlabel('Horizon de prédiction en minute (min)');
grid on;