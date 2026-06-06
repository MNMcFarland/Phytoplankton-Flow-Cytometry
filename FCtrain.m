function [mdl,cfg] = FCtrain(fin,varargin)
%%
% Train machine learning model from Flow Cytometry data files
%
% USAGE: [mdl,cfg] = FCtrain(fin,'property',value,...)
%
% INPUT:
%   fin - output from FCwtrshd.m, .fcs or file name as string, cell array of strings, or directory structure
%   Optional 'property', value pairs:
%       'mthd' - string 'nn' for neural network or 'rf' for random forest, default = 'nn' 
%       'nnet' - neural network created by FCnet.m, default = FCnet([size(X,2) 100 90 80 70 60 length(unique(Y))]);
%       'prp' - proportion of data to use for training, default = 0.2
%       'order' - height parameters and order to use for clustering, >= 4 element array, default = [7 11 10 12 8]
%       'epochs' - number of training epochs, default = 50
%       'mbs' - mini batch size, default = 50
%       'eta' - learning rate, default = 0.002
%       'radius' - neighborhood radius, default = 0.45
%       'minsz' - minimum density (number of events) in a cluster, default = 5
%       'meta' - optional logical flag for metaclustering, true = do metacluster analysis (default = false)
%       'metaradius' - metacluster neighborhood radius, default = 0.5
%       'cnf' - if not empty, plot confusion matrix up to the specified label, default = [] 
%
% OUTPUT:
%   mdl - trained model
%   cfg - structure of parameters used

%   RF - trained random forest
%   RFacc - random forest accuracy

% set defaults
cfg = struct;
cfg.mthd = 'nn';
mdl = [];
cfg.prp = .2; % proportion of data to select for training
cfg.r = 0.45; % density radius
cfg.mr = 0.6; % density meta radius
cfg.minsz = 5; % minimum cluster size
cfg.ord = [7 11 10 12 8];
cfg.cnf = [];

cfg.epochs = 50; % epochs - the number of epochs to train for
cfg.mbs = 50; % mbs - size of the mini-batches to use when sub-sampling
cfg.eta = 0.002; % eta - learning rate
% cfg.lmbda = 0; % lmbda - L2 regularization parameter

% parse optional input
if rem(length(varargin),2)~=0
    error('unmatched property, value inputs')
end
for m = 1:2:length(varargin)
    if ischar(varargin{m})
        switch lower(varargin{m})
            case 'mthd'
                cfg.mthd = varargin{m+1};
            case 'nnet'
                mdl = varargin{m+1};
            case 'prp'
                cfg.prp = varargin{m+1};
            case 'order'
                cfg.ord = varargin{m+1};
            case 'epochs'
                cfg.epochs = varargin{m+1};
            case 'mbs'
                cfg.mbs = varargin{m+1};
            case 'eta'
                cfg.eta = varargin{m+1};
            case 'radius'
                cfg.r = varargin{m+1};
            case 'minsz'
                cfg.minsz = varargin{m+1};
            case 'meta'
                cfg.meta = varargin{m+1};
            case 'metaradius'
                cfg.mr = varargin{m+1};
            case 'cnf'
                cfg.cnf = varargin{m+1};
            otherwise 
                error('FCtrain.m input error: Unrecognized input')
        end
    else
        error('FCclsfy.m input error');
    end
end

if ~isempty(fin) && ischar(fin)
    fnm = dir(fin);
elseif ~isempty(fin) && isstruct(fin) && isfield(fin,'name') && isfield(fin,'folder')
    fnm = fin;
elseif ~isempty(fin) && iscell(fin)
    for m = 1:length(fin)
        [pth, nm, ext] = fileparts(fin{m});
        fnm(m).name = [nm ext];
        if isempty(pth)
            fnm(m).folder = pwd;
        else
            fnm(m).folder = pth;
        end
    end
else % select files for training set manually
    [fname,fpath] = uigetfile('*.fcs','multiselect','on');
    if ~isempty(fname)
        for m = 1:length(fname)
            fnm(m).name = fname{m};
            fnm(m).folder = fpath;
%             dat.path{m,1} = [fpath fname{m}];
        end
    end
end

%% compile training data
% if prp<1, randomly select proportion of data, weighted according to probability density
X = [];
Y = [];
if ~isfield(fnm,'lbl')
    fnm = FCwtrshd(fnm,'radius',cfg.r,'order',cfg.ord,'printplots',false,...
        'meta',true,'metaradius',cfg.mr,'minsz',cfg.minsz,'metaminsz',1);
end

for m = 1:length(fnm)
    raw = readAccuri([fnm(m).folder filesep fnm(m).name]);
    raw = raw(:,cfg.ord);
    z = any(raw<=0,2); % identify events with zeros
    X1 = log10(1+raw(~z,:)); % log convert and remove zeros
    if cfg.prp<1
        n = size(X1,1);
        D = fnm(m).D(~z);
        W = 1-log(D)/max(log(D)); % weights for random sampling
        idx = randsample(n,round(cfg.prp*n),true,W); % samples with replacement
%         [~,idx] = datasample(X1,round(dat.prp*n),'Replace',false,'Weights',W); % samples without replacement
    else
        idx = true(size(X1,1),1);
    end
    X = [X; X1(idx,:)];
    Y1 = fnm(m).lbl(~z);
    Y = [Y; Y1(idx)];
end

clear X1 Y1 W D idx raw

X = X/7.2247; % normalize training data

% classify clusters using neural net if provided
% if strcmp(cfg.mthd,'nn') && ~isempty(mdl)
    % do some stuff like this except for clusters
    % z = any(X<=0,2); % identify events with zeros
    % lbl1 = predict(mdl,X(~z,:)); % classify
    % lbl = zeros(size(X,1),1); % initialize label vector for this file
    % lbl(~z) = lbl1; % avoid events with zeros but don't remove
% end

% split into training and validation data sets
n = size(X,1);
p = randperm(n,round(0.2*n));
Xv = X(p,:);
X(p,:) = [];
Yv = Y(p);
Y(p) = [];

%% train the model
if strcmp(cfg.mthd,'nn')
    if isempty(mdl)
        % NN = FCnet([size(X,2) 100 90 80 70 60 length(unique(Y))],[0 .2 .2 .2 .2 .2 0]); % create neural net
        mdl = FCnet([size(X,2) 100 90 80 70 60 length(unique(Y))]); % create neural net
    end
    mdl = train(mdl,X,Y,Xv,Yv,cfg.epochs,cfg.mbs,cfg.eta,0); 
elseif strcmp(cfg.mthd,'rf')
    mdl = TreeBagger(100,X,Y,'Method','classification','MinLeafSize',50,'OOBPrediction','on'); % random forest
else
    error('FCtrain error: cfg.mthd must be "nn" or "rf"')
end

if strcmp(cfg.mthd,'rf')
    pred = predict(mdl, Xv);
    pred = str2double(pred);
    cfg.acc = sum(pred==Yv)/length(Yv);
end    

% plot confusioin matrix
if ~isempty(cfg.cnf)
    pred = predict(mdl, Xv);
    figure;
    % plotconfusion(categorical(Yv(Yv<10)),categorical(pred(Yv<10)));
    plotconfusion(categorical(Yv(Yv<cfg.cnf)),categorical(pred(Yv<cfg.cnf)));
end

%% random forest
% if nargout>1
%     RF = TreeBagger(100,X,Y,'Method','classification','MinLeafSize',50,'OOBPrediction','on'); % acc = 85.1%
%     
%     % plot(oobError(RF));
%     pred = RF.predict(Xv);
%     pred = str2double(pred);
%     RFacc = sum(pred==Yv)/length(Yv);
% end

% confusioin matrix
% figure;
% plotconfusion(categorical(Yv(Yv<10)),categorical(pred(Yv<10)));


