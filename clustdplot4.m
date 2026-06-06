function ax = clustdplot4(dat,idx,varargin)

% Show results of cluster analysis for 4 parameters.
% Plots all 6 parameter combinations in a 2x3 layout.
% Requires dplot.m, mu2cmap.m
%
% USAGE: ax = clustdplot4(dat,idx,'property',value,...)
%
% INPUT:
%   dat - clustered data, n by >=4
%   idx - cluster labels, n by 1
%   Optional 'property', value pairs:
%       'lim' - max value for all axes (scalar), axes limits (2 element vector),
%           or matrix of axes limits (2 x 4) where min = row 1 and max = row 2
%       'hdr' - axes labels (cell)
%       'scale' - axes scales, 'lin' (default) or 'log'
%
% OUTPUT: plot
%   ax - axes handles

% Malcolm McFarland, version 2021
% Harbor Branch Oceanographic Institute, Florida Atlantic University
% mmcfarland@fau.edu

% parse input
if rem(length(varargin),2)~=0
    error('clustdplot4: unmatched property, value inputs')
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
            otherwise 
                error('clustdplot4.m input error: Unrecognized input')
        end
    else
        error('clustdplot4.m input error');
    end
end

ms=3;
if ~exist('lim','var') || isempty(lim) 
    mn=min(dat);
    mx=max(dat); 
elseif numel(lim)==1
%     mn=zeros(1,4);
    mn=min(dat)*.9;
    mx=repmat(lim,1,4);
elseif numel(lim)==2
    mn=repmat(lim(1),1,4); 
    mx=repmat(lim(2),1,4);
elseif size(lim,1)==2 && size(lim,2)>=4
    mn=lim(1,:);
    mx=lim(2,:);
end
if ~exist('hdr','var') || isempty(hdr); hdr = {'par 1','par 2','par 3','par 4'}; end
if ~exist('scl','var') || isempty(scl); scl = 'lin'; end
if ~exist('idx','var') || isempty(idx); idx = zeros(size(dat,1),1); end
    
grps=unique(idx(idx>0)); % groups vector for text labels

mu=zeros(length(grps),4);
for m=1:length(grps)
    mu(m,:)=mean(dat(idx==grps(m),1:4),1); % cluster means
end

if strcmpi(scl,'lin')
    cmap=mu2cmap(mu(:,1:3),mx(1:3),mn(1:3)); % assign cluster colors according to means
elseif strcmpi(scl,'log')
    cmap=mu2cmap(log10(1+mu(:,1:3)),log10(1+mx(1:3)),log10(1+mn(1:3))); % assign cluster colors according to means for log scale
else
    error('clusPlot4 input error: scale must be lin or log')
end

aw=.27;
ah=.4;

figure('position',[20 40 1150 700])%,'renderer','opengl');
ax(1)=axes('position',[.045 .58 aw ah]);
ax(2)=axes('position',[.3825 .58 aw ah]);
ax(3)=axes('position',[.72 .58 aw ah]);
ax(4)=axes('position',[.045 .075 aw ah]);
ax(5)=axes('position',[.3825 .075 aw ah]);
ax(6)=axes('position',[.72 .075 aw ah]);

set(ax,'xscale',scl,'yscale',scl) % for log scale

par=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % parameter combinations for each plot

for m=1:6
    axes(ax(m))
    [h,~,D0] = dplot(dat(:,par(m,1)), dat(:,par(m,2)),'xscale', scl, 'yscale', scl,'msize', ms);
    if ~isempty(h)
        D0 = log10(D0+1);
        h.CData = [1 1 1] .* D0/max(D0);
    end
    hold on
    for n=1:length(grps)
        [h,~,D] = dplot(dat(idx==grps(n),par(m,1)), dat(idx==grps(n),par(m,2)), 'xscale', scl, 'yscale', scl,...
            'msize', ms);
        if ~isempty(h)
            D = log10(D+1);
            c = cmap(n,:);
            h.CData = c + (1-c) .* D/max(D0);
        end
    end
    set(ax(m),'color','none','box','on','xlim',[mn(par(m,1)) mx(par(m,1))],'ylim',[mn(par(m,2)) mx(par(m,2))],'layer','top')
    text(mu(:,par(m,1)),mu(:,par(m,2)),num2str(grps),'color','k','horizontalalignment','center',...
        'backgroundcolor','none','margin',0.01,'fontweight','bold','fontsize',10);
    xlabel(hdr{par(m,1)},'interpreter','none')
    ylabel(hdr{par(m,2)},'interpreter','none')
    hold off
end

