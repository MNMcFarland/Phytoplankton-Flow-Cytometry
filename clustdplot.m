function ax = clustdplot(dat,idx,varargin)

% Show results of cluster analysis for 2 parameters.
% Requires mu2cmap.m
%
% USAGE: ax = clustdplot(dat,idx,'property',value,...)
%
% INPUT:
%   dat - clustered data, n by >=2
%   idx - cluster labels, n by 1
%   Optional 'property', value pairs:
%       'lim' - max value for each axis (scalar), axes limits (2 element vector),
%           or matrix of axes limits (2 x 2) where min = row 1 and max = row 2
%       'hdr' - axes labels (cell)
%       'scale' - axes scales, 'lin' (default) or 'log'
%       'msize' - marker size, default = 5
%       'cmap' - colors for clusters
%       'pltlbl' - plot cluster labels (logical), default = true
%
% OUTPUT: plot
%   ax - axes handles

% Malcolm McFarland, version 2021
% Harbor Branch Oceanographic Institute, Florida Atlantic University
% mmcfarland@fau.edu

%% set defaults
if ~exist('idx','var') || isempty(idx); idx = zeros(size(dat,1),1); end
hdr = {'par 1','par 2','par 3','par 4'};
scl = 'lin';
ms = 5;
plb = true;

%% parse input
if rem(length(varargin),2)~=0
    error('clustdplot: unmatched property, value inputs')
end
for m = 1:2:length(varargin)
    if ischar(varargin{m})
        switch lower(varargin{m})
            case 'lim'
                lim = varargin{m+1};
            case 'hdr'
                hdr = varargin{m+1};
            case 'scale'
                scl = varargin{m+1};
                if ~any(strcmpi(scl,{'lin','log'}))
                    error('clustdplot input error: scale must be lin or log'); 
                end
            case 'msize'
                ms = varargin{m+1};
            case 'cmap'
                cmap = varargin{m+1};
            case 'pltlbl'
                plb = varargin{m+1};
            otherwise 
                error('clustdplot input error: Unrecognized property')
        end
    else
        error('clustdplot input error');
    end
end

% if ~exist('lim','var') || isempty(lim) 
    mn = min(dat)*.95;
    mx = max(dat)*1.05; 
% elseif numel(lim)==1
if exist('lim','var') && ~isempty(lim) 
    if numel(lim)==1
    %     mn=zeros(1,2);
    %     mn = min(dat)*.9;
        mx = repmat(lim,1,size(dat,2));
    elseif numel(lim)==2
        mn = repmat(lim(1),1,size(dat,2)); 
        mx = repmat(lim(2),1,size(dat,2));
    elseif size(lim,1)==2 && size(lim,2)>=2
        mn(1:size(lim,2)) = lim(1,:);
        mx(1:size(lim,2)) = lim(2,:);
    end
end
% if ~exist('idx','var') || isempty(idx); idx = zeros(size(dat,1),1); end
% if ~exist('hdr','var') || isempty(hdr); hdr = {'par 1','par 2','par 3','par 4'}; end
% if ~exist('scl','var') || isempty(scl); scl = 'lin'; end
% if ~exist('ms','var') || isempty(ms); ms = 5; end
% if ~exist('plb','var') || isempty(plb); plb = true; end

%%
grps = unique(idx(idx>0)); % groups vector

mu = zeros(max(grps),size(dat,2));
for m = grps'%1:length(grps)
%     mu(m,:) = mean(dat(idx==grps(m),:),1); % cluster means
    mu(m,:) = mean(dat(idx==m,:),1); % cluster means
end

if ~exist('cmap','var') || isempty(cmap)
    if size(dat,2)>2
        if strcmpi(scl,'lin')
            cmap = mu2cmap(mu(:,1:3),mx(1:3),mn(1:3)); % assign cluster colors according to means
        elseif strcmpi(scl,'log')
            cmap = mu2cmap(log10(1+mu(:,1:3)),log10(1+mx(1:3)),log10(1+mn(1:3))); % assign cluster colors according to means for log scale
        else
            error('clustdplot input error: scale must be lin or log')
        end
    else
        cmap = lines(max(grps));
    end
end

% figure('position',[20 40 900 700])%,'renderer','opengl');
ax = gca;

set(ax,'xscale',scl,'yscale',scl) % for log scale

[h,~,D0] = dplot(dat(:,1), dat(:,2),'xscale', scl, 'yscale', scl,'msize',ms);
if ~isempty(h)
    D0 = log10(D0+1);
    h.CData = [1 1 1] .* D0/max(D0);
end
hold on
% delete(h) % temporary tweek, delete this!
for n = grps'%1:length(grps)
%     [h,~,D] = dplot(dat(idx==grps(n),1), dat(idx==grps(n),2), 'xscale', scl, 'yscale', scl,'msize',ms);
    [h,~,D] = dplot(dat(idx==n,1), dat(idx==n,2), 'xscale', scl, 'yscale', scl,'msize',ms);
    if ~isempty(h)
        D = log10(D+1);
        c = cmap(n,:);
        h.CData = c + (1-c) .* D/max(D0);
    end
end
set(ax,'color','none','box','on','xlim',[mn(1) mx(1)],'ylim',[mn(2) mx(2)],'layer','top')
if plb
    text(mu(grps,1),mu(grps,2),num2str(grps),'color','k','horizontalalignment','center',...
        'backgroundcolor','none','margin',0.01,'fontweight','bold','fontsize',10);
end
xlabel(hdr{1},'interpreter','none')
ylabel(hdr{2},'interpreter','none')
hold off

