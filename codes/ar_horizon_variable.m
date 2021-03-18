% script_SVM_Pred.m

try BATCHMODE=BATCHMODE; catch BATCHMODE=0; end

if ~BATCHMODE
    clear all;
    BATCHMODE=0;
    PPATH_aller_10000 = 'D:\prediction_trafic_reseau\data-07-2018\simulation_aller\multi-app-simulator-10000\default\';
    addpath(genpath([PPATH_aller_10000,'Codes_LOT2_cloud']))
    
    
    CROSS_PRED=1; %% PREDICT y{2} from y{1}
    
    %--- signals
   
  %AGGR=3;     % aggregate over windows of size AGGR (0 or 1 = no aggretation)
   % AGGR = 5;
    AGGR = 1;
     
    DSA=0;      % downsample aggregate by factor DSA  (0 = no downsampling)
    
    DDIFF=0;    % work with difference signal
   
   [data]=textread([PPATH_aller_10000,'tx_throughput.txt'],['%n']); SigStr='flow_750_3000s'; SVMP=10;
    y=data(:)'; dt=100; clear data A
    if DDIFF; y=diff(y); end;[y,ta] = agregate_data(y,AGGR,DSA,dt);

        
    
    y0=y;
    yy=y; 
    if ~iscell(y)
        ymin=min(y); y=y-min(y); ymax=max(y);y=y/max(y);ymean=mean(y); y=y-mean(y); y0=y(:);
        N=length(y);
        W=y;
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
   % D20 = 6;
    D20 = 1:10;
    D1=D10;D2=D20;
   
    
end


%--- AR model order
PP=10:20;

% --- cross validation for model order selection
clear SD* CC*
SDR_value = []
for k = 1:length(D20)
    for pid=1:length(PP)
    [~, SDR,CX] = ARprediction(y,NtrainCV,NtestCV,D1,D2(k),PP(pid));
    SDtr(pid)=SDR.train;
    SD(pid)  =SDR.test;
    CCtr(pid)=CX.train;
    CC(pid)  =CX.test;
    end




    % select best model order
    [a,optid]=max(SD);
    Pstar=PP(optid);
    
    % --- predict training and test data
    [py,SDR,CX,ID] = ARprediction(y,Ntrain,Ntest,D1,D2(k),Pstar);
    n = SDR
    sdr_value = [SDR_value n.test]
    % --- baseline prediction
    [inputs,targets,inputsT,targetsT,id,targetsall] = vectorize_data(y,Pstar,Ntrain,Ntest,D1,D2(k));
    baseline = inputsT(:,1);
    tmp=corrcoef(targetsT,baseline); CXbaseline=tmp(2);
    tmp=corrcoef(sign(targetsT),sign(baseline)); CXSbaseline=tmp(2);
    SDRbase=10*log10(mean(targetsT.^2)/mean((targetsT-baseline).^2));
 

    % --- compute average absolute errors and MAD.
    targetsTr=(y(id{2})+ymean)*ymax+ymin; 
    xtmp=(py.test+ymean)*ymax+ymin; 
    AbsError=mean(abs(targetsTr-xtmp)); xerr=targetsTr-xtmp;

    MAD_prediction = mad(xtmp);
    targetsTr=(targetsT+ymean)*ymax+ymin; xtmp=(baseline+ymean)*ymax+ymin; AbsErrorbase=mean(abs(targetsTr-xtmp)); xerrbase=targetsTr-xtmp;
    MAD_baseline = mad(xtmp);
end
% if ~BATCHMODE
% 
%         
%         
%     
%     disp('Prediction:');
%     SDRtest=SDR.test,AbsError,MAD_prediction%coef_corr_prediction
%     disp('Baseline:');
%     SDRbase,AbsErrorbase,MAD_baseline %coef_corr_baseline
%     
%      
% 
%     
% 
%     
%     
%     
% else
%     % --- save
%     save([SPTH,BASESTR,'_AR.mat'],'py', 'SDR', 'CX', 'ID', 'baseline', 'CXbaseline', 'SDRbase', 'id', 'targetsall', ...
%         'SDtr', 'SD', 'CCtr', 'CC', 'PP', 'Pstar');
% end
% 
% 
