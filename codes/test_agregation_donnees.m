%% test agrégation des données.

AGGR = 10;
DSA=0; 
DDIFF=0;
dt = 100;
PPATH_aller_500 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-500\default\';
addpath(genpath([PPATH_aller_500,'Codes_LOT2_cloud']))
[data]=textread([PPATH_aller_500,'rx_throughput.txt'],['%n']);
%[data]=textread([PPATH_aller_500,'tt.txt'],['%n']);
y=data(:)'; %dt=100;  clear data A
y1 = y;
y=y(:);y=y(1:floor(length(y)/AGGR)*AGGR);
y2 = y;
yi=cumsum(y);
yy1=diff(yi(1+AGGR:AGGR:end))';
yy2=diff(yi(2+AGGR:AGGR:end))';
yy3=diff(yi(3+AGGR:AGGR:end))';
yy4=diff(yi(4+AGGR:AGGR:end))';
yy5=diff(yi(5+AGGR:AGGR:end))';
yy6=diff(yi(6+AGGR:AGGR:end))';
yy7=diff(yi(7+AGGR:AGGR:end))';
yy8=diff(yi(8+AGGR:AGGR:end))';
yy9=diff(yi(9+AGGR:AGGR:end))';
yy10=diff(yi(10+AGGR:AGGR:end))';
YY = [yy1;yy2;yy3;yy4;yy5;yy6;yy7;yy8;yy9;yy10];
ya = YY(:);
ta=(1/2+(0:length(ya)-1))*AGGR*dt + dt;
% for k=1:AGGR;
%         yy(k,:)=diff(yi(k+AGGR:AGGR:end));
% end
% ya=yy(:);

%%if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);