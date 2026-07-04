function [dat,cpm,apm,cfg] = FCwtrshd(fin, varargin)

% Cluster flow cytometer data using density based watershed algorithm.
% Does optional 2 pass metaclustering of initial cluster 
% means and relabels clusters accordingly.
% Requires: wtrshd.m, clustdplot4.m, & readAccuri.m
%
% USAGE: [dat,cpm,apm,cfg] = FCwtrshd(fin, 'property', value,...)
%
% EXAMPLES:
%   dat = FCwtrshd('project\FCdata\*.fcs');
%
%   d = dir('*.fcs');
%   dat = FCwtrshd(d, 'radius', .3, 'meta', false); % don't metacluster
%
% INPUT:
%   fin - .fcs file names: path as string, cell array of strings, or directory structure
%   Optional 'property', value pairs:
%       'radius' - neighborhood radius, default = 0.3
%       'minsz' - minimum density (number of events) in a cluster, default = 5
%       'meta' - optional logical flag for metaclustering, true = do metacluster analysis (default = false)
%       'metaradius' - metacluster neighborhood radius, default = 0.5
%       'metaminsz' - metacluster minimum density, default = 1
%       'printplots' - optional logical flag to print plots to disk, true (default) = write to disk,
%           or directory as character string to print plots to
%       'order' - height parameters and order to use for clustering, >=4 element vector, default = [7 11 10 12 8]
%
% OUTPUT: plots of clusters or metaclusters in folder \plots
%   dat - structure with fields for each input file:
%       lbl - cluster labels
%       hdr - column headers (parameter names)
%       mu - initial cluster means for each file
%       D - density for each event
%       name - file name
%       folder - folder name
%       mu - cluster means
%       sigma - cluster covariance matrices
%       cnt - # events for each cluster
%   cpm - cells per mL for every file and cluster
%   apm - area per mL for every file, cluster, and channel
%   cfg - structure with input options used

% Malcolm McFarland
% Harbor Branch Oceanographic Institute
% Florida Atlantic University
% mmcfarland@fau.edu

% set defaults
r = 0.3;
minsz = 5;
meta = false;
mr = 0.5;
mminsz = 1;
plt = true;
ord = [7 11 10 12 8];

% parse optional input
if rem(length(varargin),2)~=0
    error('unmatched property, value inputs')
end
for m = 1:2:length(varargin)
    if ischar(varargin{m})
        switch lower(varargin{m})
            case 'radius'
                r = varargin{m+1};
            case 'minsz'
                minsz = varargin{m+1};
            case 'meta'
                meta = varargin{m+1};
            case 'metaradius'
                mr = varargin{m+1};
            case 'metaminsz'
                mminsz = varargin{m+1};
            case 'printplots'
                plt = varargin{m+1};
            case 'order'
                ord = varargin{m+1};
            otherwise 
                error('FCwtrshd.m input error: Unrecognized input')
        end
    else
        error('FCwtrshd.m input error');
    end
end

if ischar(fin)
    fnm = dir(fin);
elseif isstruct(fin) && isfield(fin,'name') && isfield(fin,'folder')
    fnm = fin;
elseif iscell(fin)
    for m = 1:length(fin)
        [pth, nm, ext] = fileparts(fin{m});
        fnm(m).name = [nm ext];
        if isempty(pth)
            fnm(m).folder = pwd;
        else
            fnm(m).folder = pth;
        end
    end
end

% include configuration options in output
cfg.r = r;
cfg.minsz = minsz;
cfg.meta = meta;
cfg.mr = mr;
cfg.mminsz = mminsz;
cfg.order = ord;
cfg.plt = plt;
fnum = length(fnm);

%% set up a progress bar with cancel button
wb = waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','classifying data');
wb.Position = [985 374 270 62.75];
hnd = findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
cntr = 0; % counter for progress bar
if meta && fnum>1
    tot = 2*fnum; % total # of files to process in 2 passes
else
    tot = fnum;
end

%% initial cluster
for m = 1:fnum
    if getappdata(wb,'canceling') % check if cancelled
        delete(wb)
        error('canceled')
    end
        
    [X,vol,hdr] = FC_read([fnm(m).folder filesep fnm(m).name], ord);
    z = any(X<=0,2); % identify events with zeros
%     Xs = zscore(X(~z,:)); % standardize variables
%     [lbl1,D,mu,sigma,cnt] = wtrshd(Xs,r,minsz); % <------------------- watershed on standardized variables
%     [lbl1,D,mu,sigma,cnt] = wtrshd(X(~z,:),r,minsz); % <------------------- watershed
    [lbl1,D1,mu] = wtrshd(X(~z,:),r,minsz); % <------------------- watershed
    lbl = zeros(size(X,1),1); % initialize lbl vector for this file
    lbl(~z) = lbl1; % avoid events with zeros but don't remove
    D = zeros(size(X,1),1);
    D(~z) = D1;
    
    if ~meta || fnum==1
%         clustplot4(X, lbl{m}, 'lim', [2 2.9 0 0; 7.23 7.23 7.23 7.23], 'hdr', hdr);
        if plt
            clustdplot4(X, lbl, 'lim', 7.23, 'hdr', hdr);
            if ischar(plt)
                if ~isfolder(plt)
                    mkdir(plt);
                end
                j = max(strfind(fnm(m).folder,filesep))+1;
                print([plt '\' fnm(m).folder(j:end) '_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300'); 
            else
                if ~isfolder([fnm(m).folder '\plots'])
                    mkdir([fnm(m).folder '\plots']);
                end
%                 clustdplot4(X, lbl, 'lim', 7.23, 'hdr', hdr);
    %             print(['plots\' fnm(m).folder(j:end) '_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300');
                print([fnm(m).folder '\plots\' fnm(m).name(1:end-4) '.png'],'-dpng','-r300');
            end
            close
        end
        if fnum>1; close; end
    end
    
    % update progress bar
    cntr = cntr + 1;
    waitbar(cntr/tot,wb,['file: ' num2str(m) ' of ' num2str(fnum)]) 
    
    dat(m).hdr = hdr;
    dat(m).lbl = lbl;
    dat(m).D = D;
    dat(m).vol = vol;
    dat(m).name = fnm(m).name;
    dat(m).folder = fnm(m).folder;
    dat(m).mu = mu; % only for original clustering
%     dat(m).sigma = sigma; % only for original clustering
%     dat(m).cnt = cnt;
end

%% meta cluster means
if meta && fnum>1
    MU = [];
%     lti = tril(ones(length(ord)))>0;
%     ltn = sum(lti,'all');
%     SIGMA = [];
    for m = 1:fnum
        MU = [MU; dat(m).mu]; % concatenate cluster means
%         k = size(dat(m).sigma,3);
%         LTI = repmat(lti,1,1,k);
%         sigma = dat(m).sigma(LTI);
%         SIGMA = [SIGMA; reshape(sigma,[k ltn])]; % concatenate unwrapped covariance matrices as row vectors
    end
%     MU = zscore([MU SIGMA]);
    mlb = wtrshd(MU,mr,mminsz,0); % <------------------- watershed
    clustdplot4(MU, mlb, 'lim', 7.23, 'hdr', hdr);
    if ischar(plt)
        if ~isfolder(plt)
            mkdir(plt);
        end
        print([plt filesep 'metacluster.png'],'-dpng','-r300'); 
    elseif islogical(plt) && plt
        print('metacluster.png','-dpng','-r300'); 
    end
    close; 

    % relabel
    pos = 1;
    for m = 1:fnum
        if getappdata(wb,'canceling') % check if cancelled
            delete(wb)
            error('canceled')
        end
        
        n = size(dat(m).mu,1);
        p = mlb(pos:pos+n-1);
        lbl = dat(m).lbl;
        rlb = zeros(size(lbl)); 
        for n = unique(lbl(lbl>0)')
            rlb(lbl==n) = p(n);
        end
        pos = pos + n;

        if plt
            X = FC_read([fnm(m).folder filesep fnm(m).name], ord);
            clustdplot4(X, rlb,'lim', 7.23, 'hdr', hdr);
            if ischar(plt)
                if ~isfolder(plt)
                    mkdir(plt);
                end
                j = max(strfind(fnm(m).folder,filesep))+1;
                print([plt '\meta_' fnm(m).folder(j:end) '_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300'); 
            else
                if ~isfolder([fnm(m).folder '\plots'])
                    mkdir([fnm(m).folder '\plots']);
                end
    %             clustplot4(X, rlb{m}, 'lim', [2 2.9 0 0; 7.23 7.23 7.23 7.23], 'hdr', hdr);
    %             print(['plots\meta_' fnm(m).folder(j:end) '_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300'); 
                print([fnm(m).folder '\plots\meta_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300'); 
            end
            close 
        end
%         if plt; print(['train\meta_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300'); close; end
%         print(['plots\meta_' fnm(m).name(5:end-4) '_' fnm(m).name(1:3) '.png'],'-dpng','-r300');
%         close

        % update progress bar
        cntr = cntr + 1;
        waitbar(cntr/tot,wb,['file: ' num2str(m) ' of ' num2str(fnum)]) 
        
        dat(m).lbl = rlb;
    end
end

%% compute cells per mL and area data, recompute mu and sigma for relabelled groups
ngrp = max(arrayfun(@(a) max(a.lbl), dat));
nd = length(ord);
cpm = zeros(fnum, ngrp);
apm = zeros(fnum, ngrp, nd);
for m = 1:fnum
    mu = zeros(ngrp,nd);
    sigma = zeros(nd,nd,ngrp);
    cnt = zeros(ngrp,1);
    [X,vol,hdr] = FC_read([fnm(m).folder filesep fnm(m).name], 1:12, 'A');
    X(:,ord) = log10(1+X(:,ord));
    for n = unique(dat(m).lbl(dat(m).lbl>0)')
        idx = dat(m).lbl==n;
        mu(n,:) = mean(X(idx,ord));
        cnt(n) = sum(idx);
        if cnt(n)>=minsz
            sigma(:,:,n) = cov(X(idx,ord));
        else
            dat(m).lbl(idx) = 0;
            mu(n,:) = zeros(1,nd);
            cnt(n) = 0;
            sigma(:,:,n) = zeros(nd);
            idx = [];
%             sigma(:,:,n) = eye(nd)*var(X(idx,ord));
        end
        cpm(m,n) = cnt(n)/vol*1000;
        apm1 = sum(X(idx,ord-6),1)/vol*1000; % sum area per mL for each channel
        apm(m,n,:) = reshape(apm1,[1 1 nd]);
    end
%     dat(m).gm = gmdistribution(mu,sigma,cnt/sum(cnt));
    dat(m).mu = mu;
    dat(m).sigma = sigma;
    dat(m).cnt = cnt;
    dat(m).hdr = hdr(ord);
end

delete(wb)

function [X,vol,hdr] = FC_read(name,ord,typ)

[dat,vol,hdr] = readAccuri(name);

if istable(dat); dat = table2array(dat); end

if exist('typ','var') && typ=='A'
    hdr = hdr(ord);
    X = dat(:,ord);
else
    hdr = hdr(ord);
    X = log10(1+dat(:,ord)); % log convert
end
    

