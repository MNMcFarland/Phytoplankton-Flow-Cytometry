function [dat,cpm,apm] = FCclsfy(fin,mdl,varargin)

% Classify flow cytometer data using trained ML model
% Use with FCnet.m, Random Forest, Discriminant Analysis, etc.
% Model (mdl) must have predict function
%
% USAGE: [dat,cpm,apm] = FCclsfy(fin,mdl,'property',value,... )
%
% INPUT:
%   fin - .fcs file names, string, cell array of strings, or directory structure
%   mdl - trained model object (class)
%   Optional 'property', value pairs:
%       'ord' - parameters/columns and order of input from files, default = [7 11 10 12 8]
%       'minprb' - minimum classification probablility, default = 0
%       'printplots' - optional logical flag to print plots to disk, true (default) = write to disk,
%           or directory as character string to print plots to
%
% OUTPUT:
%   dat - structure with fields for each file:
%           lbl - cluster label for each event
%           prb - probability for each event
%           vol - volume
%           name - file names
%           folder - folder for each file
%           hdr - parameter names from .fcs data files
%           mu - cluster means
%           sigma - covariance matrix for each cluster
%           cnt - event count for each cluster
%   cpm - cells per mL for every file x cluster
%   apm - area per mL for every file x cluster x parameter

% defaults
ord = [7 11 10 12 8];
minprb = 0;
plt = true;

% parse optional input
if rem(length(varargin),2)~=0
    error('unmatched property, value inputs')
end
for m = 1:2:length(varargin)
    if ischar(varargin{m})
        switch lower(varargin{m})
            case 'minprb'
                minprb = varargin{m+1};
            case 'printplots'
                plt = varargin{m+1};
            case 'ord'
                ord = varargin{m+1};
            otherwise 
                error('FCclsfy.m input error: Unrecognized input')
        end
    else
        error('FCclsfy.m input error');
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

fnum = length(fnm);

%% set up a progress bar with cancel button
if plt
    wb = waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
        'name','classifying data');%,'position',[985 374 270 62.75]);
    wb.Position = [785 374 270 62.75];
    hnd = findobj(wb,'type','axes');
    set(get(hnd,'title'),'interpreter','none')
    setappdata(wb,'canceling',0)
    cntr = 0; % counter for progress bar
    tot = length(fnm); % total # of files to process
end
%%
for m = 1:length(fnm)
%     if getappdata(wb,'canceling') % check if cancelled
%         delete(wb)
%         error('cancelled')
%     end
        
    [X1,vol,hdr] = FC_read([fnm(m).folder filesep fnm(m).name],ord); % read data
    X = X1/7.2247; % normalize
    z = any(X1<=0,2); % identify events with zeros
    [lbl1,prb1] = predict(mdl,X(~z,:)); % classify
    if iscell(lbl1); lbl1 = str2double(lbl1); end % RF model output is cell array of char vectors
    prb1 = max(prb1,[],2);
    lbl = zeros(size(X,1),1); % initialize label vector for this file
    prb = zeros(size(X,1),1); % initialize probability vector for this file
    lbl(~z) = lbl1; % avoid events with zeros but don't remove
    prb(~z) = prb1;
    
    lbl(prb<minprb) = 0;

%     mu = zeros(max(lbl),size(X1,2));
%     for n = unique(lbl(lbl>0)')
%         mu(n,:) = mean(X1(lbl==n,:),1); % cluster means
%     end
    
%     clustdplot4(X1, lbl, 'lim', 7.23, 'hdr', hdr);
%         clustplot4(X, lbl{m}, 'lim', [2 2.9 0 0; 7.23 7.23 7.23 7.23], 'hdr', hdr);
    if plt
        if getappdata(wb,'canceling') % check if cancelled
            delete(wb)
            error('cancelled')
        end

        clustdplot4(X1, lbl, 'lim', 7.23, 'hdr', hdr);
        if ischar(plt)
            if ~isfolder(plt)
                mkdir(plt);
            end
            j = max(strfind(fnm(m).folder,filesep))+1;
            print([plt '\pred_' fnm(m).folder(j:end) '_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300'); 
        else
            if ~isfolder([fnm(m).folder '\plots'])
                mkdir([fnm(m).folder '\plots']);
            end
%             print(['plots\' fnm(m).folder(j:end) '_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300');
            print([fnm(m).folder '\plots\pred_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300');
        end
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(m) ' of ' num2str(length(fnm))]) 
    end
    if fnum>1; close; end
    
%     print(['plots\mdl_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300');
%     if ~isfolder([fnm(m).folder '\plots'])
%         mkdir([fnm(m).folder '\plots']);
%     end
%     print([fnm(m).folder '\plots\mdl_' fnm(m).name(1:end-4) '.png'],'-dpng','-r300');
%     close
    
    dat(m).lbl = lbl;
    dat(m).prb = prb;
    dat(m).vol = vol;
    dat(m).name = fnm(m).name;
    dat(m).folder = fnm(m).folder;
    dat(m).hdr = hdr;
%     dat(m).mu = mu;
end

%% compute group means, cells per mL, and area per mL
if nargout>1
    ngrp = max(arrayfun(@(a) max(a.lbl), dat));
    nd = length(ord);
    cpm = zeros(fnum, ngrp);
    apm = zeros(fnum, ngrp, nd);
    for m = 1:fnum
        mu = zeros(ngrp,nd);
        sigma = zeros(nd,nd,ngrp);
        cnt = zeros(ngrp,1);
        X = FC_read([fnm(m).folder filesep fnm(m).name], 1:12, 'A'); 
        X(:,ord) = log10(1+X(:,ord));
        for n = unique(dat(m).lbl(dat(m).lbl>0)') % 1:length(mdl.lblnames)
            idx = dat(m).lbl==n;
            mu(n,:) = mean(X(idx,ord),1);
            cnt(n) = sum(idx);
            if cnt(n)>1
                sigma(:,:,n) = cov(X(idx,ord));
            else
                sigma(:,:,n) = zeros(nd);
    %             sigma(:,:,n) = eye(nd);
            end
            cpm(m,n) = cnt(n)/dat(m).vol*1000;
            apm1 = sum(X(idx,ord-6),1)/dat(m).vol*1000; % sum area per mL for each channel
            apm(m,n,:) = reshape(apm1,[1 1 nd]);
%             for p = 1:size(X,2)
%                 apm(m,n,p) = sum(X(dat(m).lbl==n,p))/dat(m).vol*1000; % sum area per mL for each channel
%             end
        end
        dat(m).mu = mu;
        dat(m).sigma = sigma;
        dat(m).cnt = cnt;
%         dat(m).hdr = hdr(ord);
    end
end
if exist('wb','var'); delete(wb); end

function [X,vol,hdr] = FC_read(name,ord,typ)

[dat,vol,hdr] = readAccuri(name);

if exist('typ','var') && typ=='A'
    hdr = hdr(ord);
    X = dat(:,ord);
else
    hdr = hdr(ord);
    X = log10(1+dat(:,ord)); % log convert
end


