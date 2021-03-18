% script_SVM_Pred.m



try BATCHMODE=BATCHMODE; catch BATCHMODE=0; end

if ~BATCHMODE
    clear all;
    BATCHMODE=0;
    
    %PPATH='C:\Users\uayesta\Dropbox\RT-CNES-HW-JYT\data-06-2018\';
    %PPATH='/Users/hwendt/Dropbox/RT-CNES-HW-JYT/data-06-2018/';
    %PPATH = 'D:\prediction_trafic_reseau\data-07-2018\';
    
    %addpath(genpath('C:\Users\uayesta\Dropbox\RT-CNES-HW-JYT\data-06-2018\Codes_LOT2_cloud'))
    PPATH = 'D:\prediction_trafic_reseau\data-07-2018\';
    PPATH_aller_500 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-500\default\';
    PPATH_aller_750 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-750\default\';
    PPATH_aller_5000 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-5000\default\';
    PPATH_aller_10000 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-10000\default\';
    
    
    addpath(genpath([PPATH,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_500,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_750,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_5000,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_10000,'Codes_LOT2_cloud']))
    
    CROSS_PRED=1; %% PREDICT y{2} from y{1}
    
    %--- signals
    DATASEL=1;  % select data set
    AGGR=0;     % aggregate over windows of size AGGR (0 or 1 = no aggretation)
    DSA=0;      % downsample aggregate by factor DSA  (0 = no downsampling)
    DDIFF=0;    % work with difference signal
    switch DATASEL
        case 0;
            N=27998;dt=100; % signal length
            NF=20; % NF=20;
            rng(1234,'twister'); b = fir1(NF-1, .5); [d,p0] = lpc(b,NF); u = sqrt(p0)*randn(1,N);
            y = filter(1,d,u);
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 1
            [data]=textread([PPATH,'aller.txt'],['%n']); SigStr='flow_500_3000s'; SVMP=1;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 2
            [data]=textread([PPATH,'retour.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=2;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 3
Cshift=6;  % shift 2nd time series by Cshift to compensate that measurements are taken at SAT, not GW
            [data]=textread([PPATH,'retour.txt'],['%n']); SigStr='flow_500_3000s'; SVMP=1;
            y{1}=data(:)'; dt=100; clear data A
            if Cshift;y{1}=y{1}(Cshift+1:end);end
            if DDIFF; y{1}=diff(y{1}); end;[y{1},ta] = agregate_data(y{1},AGGR,DSA,dt);
            [data]=textread([PPATH,'aller.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=2;
            y{2}=data(:)'; dt=100; clear data A
            if Cshift;y{2}=y{2}(1:end-Cshift);end
            if DDIFF; y{2}=diff(y{2}); end;[y{2},ta] = agregate_data(y{2},AGGR,DSA,dt);
        case 4
            [data]=textread([PPATH,'aller.txt'],['%n']); SigStr='flow_500_3000s'; SVMP=1;
            y{1}=data(:)'; dt=100; clear data A
            if DDIFF; y{1}=diff(y{1}); end;[y{1},ta] = agregate_data(y{1},AGGR,DSA,dt);
            [data]=textread([PPATH,'retour.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=2;
            y{2}=data(:)'; dt=100; clear data A
            if DDIFF; y{2}=diff(y{2}); end;[y{2},ta] = agregate_data(y{2},AGGR,DSA,dt);
        case 5
            [data]=textread([PPATH_aller_500,'rx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=3;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 6
            [data]=textread([PPATH_aller_500,'tx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=4;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 7
            [data]=textread([PPATH_aller_750,'rx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=5;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 8
            [data]=textread([PPATH_aller_750,'tx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=6;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 9
            [data]=textread([PPATH_aller_5000,'rx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=7;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 10
            [data]=textread([PPATH_aller_5000,'tx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=8;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 11
            [data]=textread([PPATH_aller_10000,'rx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=9;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
        case 12
            [data]=textread([PPATH_aller_10000,'tx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=10;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
    end
    y0=y;
    yy=y; 
    if ~iscell(y)
        %y=y-min(y);y=y/max(y); y=y-mean(y); y0=y(:);
        ymin=min(y); y=y-min(y); ymax=max(y);y=y/max(y);ymean=mean(y); y=y-mean(y); y0=y(:);
        N=length(y);
    else
        for i=1:length(y)
            %y{i}=y{i}-min(y{i});y{i}=y{i}/max(y{i}); y{i}=y{i}-mean(y{i}); y0{i}=y{i}(:);
            ymin(i)=min(y{i});y{i}=y{i}-min(y{i});ymax(i)=max(y{i});y{i}=y{i}/max(y{i}); ymean(i)=mean(y{i});y{i}=y{i}-mean(y{i}); y0{i}=y{i}(:);
            N=length(y{i});
        end
    end
          
    
    % --- set training and test parameters
    Ntrain=10000; % training examples
    Ntest=N;
    NtrainCV=5000;
    NtestCV=5000;
    D10=1;       % embedding distance
    %D20=12;       % D2 - step ahead prediction
    %D20=100;
    D20=1;
    D1=D10;D2=D20;
    %if DSA==0; D1=AGGR*D10; D2=AGGR*D20; end
    
%     %% PREDICT y{2} from y{1}
%     if CROSS_PRED
%         y0=y;clear y; y{1}=y0;y{2}=y0.^2;
%     end
end

%--- set SVM parameters
 %DDV=[4 6 8 10 12]; % embedding dimension
 %MFACT = [0.5 1 1.5 2 3 5 10]; % data scaling
 %CV = [10 100 1000 10000];

 DDV=[4]; % embedding dimension
 MFACT = [0.5]; % data scaling
 CV = [10];

Llambda=1/2000;
Eepsilon = .001;
Kernel='poly';      Kerneloption = [1 2 3 4 5];
%Kernel='gaussian';  Kerneloption = [0.1 0.25 0.5 1 2 3 5 10];


 %MFACT = [0.5 2 5 10]; % data scaling
 %CV = [100 1000 10000];

 MFACT = [0.5]; % data scaling
 CV = [100];

%Kernel='gaussian';  Kerneloption = [0.1 0.5 2 10];



% --- cross validation for SVM parameter selection
clear SD* CC* arSD* arCC* maSD* maCC* i*
ii=1;
params.kernel=Kernel;
for did=1:length(DDV);                              DD=DDV(did);
    for mid=1:length(MFACT);                        Mfact=MFACT(mid);
        for cid=1:length(CV);                       params.C = CV(cid);
            for lid=1:length(Llambda);              params.lambda=Llambda(lid);
                for eid=1:length(Eepsilon);         params.epsilon = Eepsilon(eid);
                    for kid=1:length(Kerneloption); params.kerneloption = Kerneloption(kid);
                        
                        [~, SDR,CX] = SVMprediction(y,NtrainCV,NtestCV,D1,D2,DD,Mfact,params, 1 );
                        SDtr(ii)=SDR.train;
                        SD(ii)  =SDR.test;
                        CCtr(ii)=CX.train;
                        CC(ii)  =CX.test;
                        
                        iDD(ii)=DD;
                        iMM(ii)=Mfact;
                        iCC(ii)=params.C;
                        iLL(ii)=params.lambda;
                        iEE(ii)=params.epsilon;
                        iKK(ii)=params.kerneloption;
                        
                        ii=ii+1;
                    end;end;end;end;end;end

% select best parameter combination
[a,optid]=max(SD);
%[a,Did]=max(CC);
params0=params;clear params;
Dstar=iDD(optid);
Mfstar=iMM(optid);
params.kernel=Kernel;
params.C=iCC(optid);
params.lambda=iLL(optid);
params.epsilon=iEE(optid);
params.kerneloption=iKK(optid);

if ~BATCHMODE
    Mfstar
    params
end


% --- predict training and test data
[py,SDR,CX,ID] = SVMprediction(y,Ntrain,Ntest,D1,D2,Dstar,Mfstar,params);

% --- baseline prediction
[inputs,targets,inputsT,targetsT,id,targetsall] = vectorize_data(y,Dstar,Ntrain,Ntest,D1,D2);
baseline = inputsT(:,1);
tmp=corrcoef(targetsT,baseline); CXbaseline=tmp(2);
tmp=corrcoef(sign(targetsT),sign(baseline)); CXSbaseline=tmp(2);
SDRbase=10*log10(mean(targetsT.^2)/mean((targetsT-baseline).^2));

% --- compute average absolute errors
if ~iscell(y)
    targetsTr=(y(id{2})+ymean)*ymax+ymin; xtmp=(py.test+ymean)*ymax+ymin; AbsError=mean(abs(targetsTr-xtmp));
    %MAD_prediction = 100*(mean(abs(targetsT-xtmp))/mean(abs(targetsT))); 
    xerr=targetsTr-xtmp;
    MAD_prediction = mad(xtmp);
    targetsTr=(targetsT+ymean)*ymax+ymin; xtmp=(baseline+ymean)*ymax+ymin; AbsErrorbase=mean(abs(targetsTr-xtmp)); 
    %MAD_baseline = 100*(mean(abs(targetsT-xtmp))/mean(abs(targetsT))); 
    xerrbase=targetsTr-xtmp;
    MAD_baseline = mad(xtmp);
else
    targetsTr=(y{1}(id{2})+ymean(1))*ymax(1)+ymin(1); xtmp=(py.test+ymean(1))*ymax(1)+ymin(1); AbsError=mean(abs(targetsTr-xtmp)); 
    %MAD_prediction = 100*(mean(abs(targetsT-xtmp))/mean(abs(targetsT)));
    MAD_prediction = mad(xtmp);
    xerr=targetsTr-xtmp;
    targetsTr=(targetsT+ymean(1))*ymax(1)+ymin(1); xtmp=(baseline+ymean(1))*ymax(1)+ymin(1); AbsErrorbase=mean(abs(targetsTr-xtmp));
    %MAD_baseline = 100*(mean(abs(targetsT-xtmp))/mean(abs(targetsT))) ;
    MAD_baseline = mad(xtmp);
    xerrbase=targetsTr-xtmp;
end


% for plots
iZOOM=2000; nZOOM=100; % for plots
if ~BATCHMODE
    disp('Prediction:');
    SDRtest=SDR.test,AbsError,MAD_prediction
    disp('Baseline:');
    SDRbase,AbsErrorbase,MAD_baseline
    
    figure(111);clf;
    subplot(221); plot(ta(id{2}),xerr);axis tight;ylabel('absolute error');grid on;
    subplot(222); hist(xerr,50);grid on;
    subplot(223); plot(ta(id{2}),xerrbase);axis tight;ylabel('absolute error baseling');grid on;
    subplot(224); hist(xerrbase,50);grid on;
    
    % --- plots
    if ~iscell(y)
        figure(1); clf;
        subplot(211);plot(ta,y,'k-');  grid on;  title('Signal');axis tight;
        [R2] = ACF_fft(y);
        %subplot(212);loglog(ta(1:floor(N/20)),abs(R2(1:floor(N/20))),'k-'); grid on; hold off;  title('Signal auto-correlation');axis tight;
        subplot(212);plot(ta(1:50),(R2(1:50))/R2(1),'k-'); grid on; hold off;  title('Signal auto-correlation');axis tight;
    else
        figure(1); clf;
        subplot(211);plot(ta,y{1},'k-');  grid on;  title('Input Signal');axis tight;
        subplot(212);plot(ta,y{2},'k-');  grid on;  title('Target Signal');axis tight;
    end
    %
    figure(2);clf;
    subplot(211);
    hold on;
    plot(ta(id{3}),targetsall,'k');
    plot(ta(id{2}),py.test,'b');
    plot(ta(id{1}),py.train,'r');
    xlim(ta([1 N]));
    hold off; grid on;
    title(['SVM predict- Embedding=',num2str(Dstar),', Train=',num2str(Ntrain),', D1=',num2str(D1),', D2=',num2str(D2),'  - SDR=',num2str(SDR.test),'dB - Rho=',num2str(CX.test),'  - SDR base=',num2str(SDRbase),'dB']);
    subplot(212);
    hold on;
    plot(ta(id{3}),targetsall,'k');
    plot(ta(id{2}),py.test,'b');
    plot(ta(id{1}),py.train,'r');
    xlim(ta(iZOOM+[0 nZOOM-1]));
    hold off; grid on;
else
    % --- save
    save([SPTH,BASESTR,'_SVM_',Kernel,'.mat'],'py', 'SDR', 'CX', 'ID', 'baseline', 'CXbaseline', 'SDRbase', 'id', 'targetsall', ...
        'SDtr', 'SD', 'CCtr', 'CC', 'iDD','iMM','Mfact','iCC','iLL','iEE','iKK', 'Dstar', 'Mfstar', 'params');
end

