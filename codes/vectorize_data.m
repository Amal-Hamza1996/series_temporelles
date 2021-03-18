function [xtrain,ytrain,xtest,ytest,id,yall] = vectorize_data(y,DD,Ntrain,Ntest,D1,D2,iRD)
%% function [xtrain,ytrain,xtest,ytest,id,yall] = vectorize_data(y,DD,Ntrain,Ntest,D1,D2,iRD)

CUMULATE=1; % CUMULATE=0: if two time series are given, predict first one from second one
            % CUMULATE=1: if two time series are given, use both to predict first one

if ~iscell(y)
    y=y(:)';
    N=length(y); if nargin<3; Ntrain=floor(N/10); end; if nargin<4; Ntest=N-Ntrain; end; if nargin<5; D1=1; end; if nargin<6; D2=1; end; if nargin<7; iRD=0; end
    
    X0=y;
    for i=D1:D1:DD*D1-1; X0=[X0;circshift(y',-i)']; end
    Y0=circshift(y',-i-D2)';
    % remove last circular elements
    X0=flipud(X0(:,1:end-DD -D2+1 -(DD-1)*(D1-1)));
    Y0=Y0(1:end-DD -D2+1 -(DD-1)*(D1-1));
    Ntest=min(Ntest,length(Y0)-Ntrain);
    Ntrain=min(Ntrain,length(Y0));
    
    if iRD
        rID=randperm(Ntrain+Ntest);
    else
        rID=1:Ntrain+Ntest;
    end
    idTrain=sort(rID(1:Ntrain));
    idTest=sort(rID(Ntrain+1:Ntrain+Ntest));
    
    % training data
    xtrain=X0(:,idTrain)';
    ytrain=Y0(idTrain)';
    % test data
    xtest=X0(:,idTest)';
    ytest=Y0(idTest)';
    
    yall=[ytrain;ytest];
    
    id{1}=i+D2+(idTrain);
    id{2}=i+D2+(idTest);
    id{3}=i+D2+sort(rID);
    
    % id{1}=i+D2+(1:Ntrain);
    % id{2}=i+D2+(Ntrain+1:length(yall));
    % id{3}=i+D2+(1:length(yall));
else
    
    N=length(y{1}); if nargin<3; Ntrain=floor(N/10); end; if nargin<4; Ntest=N-Ntrain; end; if nargin<5; D1=1; end; if nargin<6; D2=1; end; if nargin<7; iRD=0; end
    
    for ii=1:2
        y{ii}=y{ii}(:)';
        X0{ii}=y{ii};
        for i=D1:D1:DD*D1-1; X0{ii}=[X0{ii};circshift(y{ii}',-i)']; end
        Y0{ii}=circshift(y{ii}',-i-D2)';
        % remove last circular elements
        X0{ii}=flipud(X0{ii}(:,1:end-DD -D2+1 -(DD-1)*(D1-1)));
        Y0{ii}=Y0{ii}(1:end-DD -D2+1 -(DD-1)*(D1-1));
        Ntest=min(Ntest,length(Y0{ii})-Ntrain);
        Ntrain=min(Ntrain,length(Y0{ii}));
    end
    
    if iRD
        rID=randperm(Ntrain+Ntest);
    else
        rID=1:Ntrain+Ntest;
    end
    idTrain=sort(rID(1:Ntrain));
    idTest=sort(rID(Ntrain+1:Ntrain+Ntest));
    
    % training data
    xtrain=X0{1}(:,idTrain)';
    ytrain=Y0{2}(idTrain)';
    % test data
    xtest=X0{1}(:,idTest)';
    ytest=Y0{2}(idTest)';
    
    % predict 1 from 1 and 2
    if CUMULATE
        % training data
        xtrain=[X0{1}(:,idTrain)',X0{2}(:,idTrain)'];
        ytrain=Y0{1}(idTrain)';
        % test data
        xtest=[X0{1}(:,idTest)',X0{2}(:,idTest)'];
        ytest=Y0{1}(idTest)';
        
        %disp('CUMUL');
    end
    
    yall=[ytrain;ytest];
    
    id{1}=i+D2+(idTrain);
    id{2}=i+D2+(idTest);
    id{3}=i+D2+sort(rID);
    
    % id{1}=i+D2+(1:Ntrain);
    % id{2}=i+D2+(Ntrain+1:length(yall));
    % id{3}=i+D2+(1:length(yall));
    
%        keyboard
        
end
