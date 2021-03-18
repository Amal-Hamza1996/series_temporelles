% script_SVM_Pred.m


try BATCHMODE=BATCHMODE; catch BATCHMODE=0; end

if ~BATCHMODE
    clear all;
    BATCHMODE=0;
    
    %PPATH='D:\Dropbox\Dropbox\RT-CNES-HW-JYT\data-07-2018\';
    PPATH = 'D:\prediction_trafic_reseau\data-07-2018\';
    PPATH_aller_500 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-500\default\';
    PPATH_aller_750 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-750\default\';
    PPATH_aller_5000 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-5000\default\';
    PPATH_aller_10000 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-10000\default\';
    PPATH_Users_variable = 'D:\prediction_trafic_reseau\data-07-2018\nb_variable_utilisateurs\scenario1\';
    
    %PPATH='C:\Users\uayesta\Dropbox\RT-CNES-HW-JYT\data-06-2018\';
    %PPATH='/Users/hwendt/Dropbox/RT-CNES-HW-JYT/data-06-2018/';
    
    %addpath(genpath('C:\Users\uayesta\Dropbox\RT-CNES-HW-JYT\data-06-2018\Codes_LOT2_cloud'))
    addpath(genpath([PPATH,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_500,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_750,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_5000,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_aller_10000,'Codes_LOT2_cloud']))
    addpath(genpath([PPATH_Users_variable,'Codes_LOT2_cloud']))
    
    CROSS_PRED=1; %% PREDICT y{2} from y{1}
    
    %--- signals
    DATASEL=13;  % select data set
    %AGGR=3;     % aggregate over windows of size AGGR (0 or 1 = no aggretation)
    AGGR=1200;
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
Cshift=6; % 7  % shift 2nd time series by Cshift to compensate that measurements are taken at SAT, not GW
% --- Cshift = 6 --> cross-correlation has maximum at 0.
            [data]=textread([PPATH,'retour.txt'],['%n']); SigStr='flow_500_3000s'; SVMP=1;
            y{1}=data(:)'; dt=100; clear data A
            %if Cshift;y{1}=y{1}(Cshift+1:end);end
            if Cshift;y{1}=circshift(y{1},-Cshift);y{1}=y{1}(abs(Cshift)+1:end-abs(Cshift));end
            if DDIFF; y{1}=diff(y{1}); end;[y{1},ta] = agregate_data(y{1},AGGR,DSA,dt);
            [data]=textread([PPATH,'aller.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=2;
            %[data]=textread([PPATH,'queueretourterminal.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=2;
            y{2}=data(:)'; dt=100; clear data A
            %if Cshift;y{2}=y{2}(1:end-Cshift);end
            if Cshift;y{2}=y{2}(abs(Cshift)+1:end-abs(Cshift));end
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
        case 13
            [data]=textread([PPATH_Users_variable,'rx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=11;
            y=data(:)'; dt=100; clear data A
            if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);
      
    end
    y0=y;
    yy=y; 
    if ~iscell(y)
        %y=y-min(y);y=y/max(y); y=y-mean(y); y0=y(:);
        ymin=min(y); y=y-min(y); ymax=max(y);y=y/max(y);ymean=mean(y); y=y-mean(y); 
        y = y + min(y);
        y0=y(:);
        N=length(y);
    else
        for i=1:length(y)
            %y{i}=y{i}-min(y{i});y{i}=y{i}/max(y{i}); y{i}=y{i}-mean(y{i}); y0{i}=y{i}(:);
            ymin(i)=min(y{i});y{i}=y{i}-min(y{i});ymax(i)=max(y{i});y{i}=y{i}/max(y{i}); ymean(i)=mean(y{i});y{i}=y{i}-mean(y{i}); y0{i}=y{i}(:);
            N=length(y{i});
        end
    end
    
    
%     x1=y{1};x2=y{2};
% cc12=xcorr(x1,x2,'coeff');cc120=cc12;cc12=cc12(ceil(length(cc12)/2):end);
% figure(7); plot((0:Mlag-1)-25,cc120(1:Mlag),'bo-');grid on;title('autocorr?lation croiss?e');ylim([0 1]);
% return
    
    % --- set training and test parameters
    Ntrain=10000; % training examples
    Ntest=N;
    NtrainCV=5000;
    NtestCV=5000;
    D10=1;       % embedding distance
    D20=1200;
    %D20=12;       % D2 - step ahead prediction
   % D20 = 6;
    D1=D10;D2=D20;
    %if DSA==0; D1=AGGR*D10; D2=AGGR*D20; end
    
%     %% PREDICT y{2} from y{1}
%     if CROSS_PRED
%         y0=y;clear y; y{1}=y0;y{2}=y0.^2;
%     end
    
end


%--- set NN parameters
DDV=[4 8 12 16]; % embedding dimension
NNHLV(:,1) = [2 3 4 5 7 10 15]; % Network size

DDV=[12]; % embedding dimension
%DDV=[24]; % embedding dimension

NNHLV(:,1) = [3]; % Network size
%NNHLV(:,1) = [5]; % Network size
%NNHLV(:,1) = [7]; % Network size

NNHLV(:,2)=NNHLV(:,1);

NNHLV(:,3)=NNHLV(:,1);
%NNHLV(:,4)=NNHLV(:,1);
%NNHLV(:,5)=NNHLV(:,1);











% --- cross validation for parameter selection
clear SD* CC* arSD* arCC* maSD* maCC* i*
ii=1;
for did=1:length(DDV)
    DD=DDV(did);
    for nid=1:length(NNHLV)
        NNHL=NNHLV(nid,:);
        [~, SDR,CX] = NNprediction(y,NtrainCV,NtestCV,D1,D2,DD,NNHL);
        SDtr(ii)=SDR.train;
        SD(ii)  =SDR.test;
        CCtr(ii)=CX.train;
        CC(ii)  =CX.test;
        
        iDD(ii)=DDV(did);
        iNN(ii,:)=NNHLV(nid,:);
        ii=ii+1;
    end
end

% select best parameter combination
[a,optid]=max(SD);
%[a,Did]=max(CC);
Dstar=iDD(optid);
NNHLstar=iNN(optid,:);
if ~BATCHMODE
    NNHLstar
end

% --- predict training and test data
[py,SDR,CX,ID] = NNprediction(y,Ntrain,Ntest,D1,D2,Dstar,NNHLstar);

% --- baseline prediction
[inputs,targets,inputsT,targetsT,id,targetsall] = vectorize_data(y,Dstar,Ntrain,Ntest,D1,D2);
baseline = inputsT(:,1);
tmp=corrcoef(targetsT,baseline); CXbaseline=tmp(2);
tmp=corrcoef(sign(targetsT),sign(baseline)); CXSbaseline=tmp(2);
SDRbase=10*log10(mean(targetsT.^2)/mean((targetsT-baseline).^2));

% --- compute average absolute errors
if ~iscell(y)
    targetsTr=(y(id{2})+ymean)*ymax+ymin; xtmp=(py.test+ymean)*ymax+ymin; 
    MAD_prediction=100*mean(abs(targetsTr-xtmp))/mean(abs(targetsTr)); 
    %MAD_prediction = mad(xtmp);
    AbsError=mean(abs(targetsTr-xtmp));
    xerr=targetsTr-xtmp;
    targetsTr=(targetsT+ymean)*ymax+ymin; xtmp=(baseline+ymean)*ymax+ymin; 
    MAD_baseline=100*mean(abs(targetsTr-xtmp))/mean(abs(targetsTr)); 
    %MAD_baseline = mad(xtmp);
    AbsErrorbase=mean(abs(targetsTr-xtmp));
    xerrbase=targetsTr-xtmp;
else
    targetsTr=(y{1}(id{2})+ymean(1))*ymax(1)+ymin(1); xtmp=(py.test+ymean(1))*ymax(1)+ymin(1); 
    MAD_prediction=100*mean(abs(targetsTr-xtmp))/mean(abs(targetsTr)); 
    %MAD_prediction = mad(xtmp);
    AbsError=mean(abs(targetsTr-xtmp));
    xerr=targetsTr-xtmp;
    targetsTr=(targetsT+ymean(1))*ymax(1)+ymin(1); xtmp=(baseline+ymean(1))*ymax(1)+ymin(1); 
    MAD_baseline=100*mean(abs(targetsTr-xtmp))/mean(abs(targetsTr)); 
    %MAD_baseline = mad(xtmp);
    AbsErrorbase=mean(abs(targetsTr-xtmp));
    xerrbase=targetsTr-xtmp;
end

% compute MAD
%MAD_baseline = 100*(mean(abs(targetsT-baseline))/mean(abs(targetsT)));
%MAD_prediction = 100*(mean(abs(targetsTr-xtmp))/mean(abs(targetsTr)));

%MAD_baseline
%MAD_prediction
% for plots
iZOOM=2000; nZOOM=100; % for plots

if ~BATCHMODE
    disp('Prediction:');
    SDRtest=(SDR.test),AbsError, MAD_prediction
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
        subplot(212);plot(ta(1:50),(R2(1:50))/R2(1),'ko-'); grid on; hold off;  title('Signal auto-correlation');axis tight;
    else
        figure(1); clf;
        subplot(211);plot(ta,y{1},'k-');  grid on;  title('Input Signal');axis tight;
        subplot(212);plot(ta,y{1},'k-');  grid on;  title('Target Signal');axis tight;
    end
    %
    figure(2);clf;
%    subplot(211);
    hold on;
    plot(ta(id{3}),targetsall,'k');
    plot(ta(id{2}),py.test,'b');
    plot(ta(id{1}),py.train,'r');
    legend('input data','test data','training data');
    xlim(ta([1 N]));
    xlabel('time (ms)');
    ylabel('normalized and mean centering data');
    hold off; grid on;
    title(['FFNN predict- Embedding=',num2str(Dstar),', Train=',num2str(Ntrain),', D1=',num2str(D1),', D2=',num2str(D2),'  - SDR=',num2str(SDRtest),'dB - Rho=',num2str(CX.test),'  - SDR base=',num2str(SDRbase),'dB']);
%     subplot(212);
%     hold on;
%     plot(ta(id{3}),targetsall,'k');
%     plot(ta(id{2}),py.test,'b');
%     plot(ta(id{1}),py.train,'r');
%     xlim(ta(iZOOM+[0 nZOOM-1]));
%     hold off; grid on;
else
    % --- save
    save([SPTH,BASESTR,'_NN.mat'],'py', 'SDR', 'CX', 'ID', 'baseline', 'CXbaseline', 'SDRbase', 'id', 'targetsall', ...
        'SDtr', 'SD', 'CCtr', 'CC', 'iDD', 'iNN','Dstar','NNHLstar');
end
