% test AR pr�diction.

%% Le chemin de donn�es. 
PPATH = 'D:\prediction_trafic_reseau\data-07-2018\';
addpath(genpath([PPATH,'Codes_LOT2_cloud']))

%% Lire les donn�es.
[data]=textread([PPATH,'retour.txt'],['%n']); 
y=data(:)';
W = y;
%% Les param�tres..
Ntrain = 10000;
Ntest = length(y);
D1 = 1; % embedding dimension.
D2 = 1; % horizon de pr�diction.
P = 10; % ordre.
AGGR = 1;
DSA = 0;
dt = 100;
%% test agr�agats de flux.
%[ya,ta] = agregate_data(y,AGGR,DSA,dt);
y=y(:);
yt = y;
y=y(1:floor(length(y)/AGGR)*AGGR);
yi=cumsum(y);
for k=1:AGGR;
    yy(k,:)=diff(yi(k+AGGR:AGGR:end));
end
ya=yy(:);
ta=(1/2+(0:length(ya)-1))*AGGR*dt + dt;

%% Les param�tres..
Ntrain = 10000;
Ntest = length(y);
D1 = 1; % embedding dimension.
D2 = 1; % horizon de pr�diction.
P = 10;
%% test AR pr�diction.
[py,SDR,CX,ID] = ARprediction(y,Ntrain,Ntest,D1,D2,P);

