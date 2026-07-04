function [h,H,Z] = dplot(X,Y,varargin)

% Bivariate density scatter plot
% Uses histcounts2.m
%
% USAGE: [h,H,Z]=dplot(X,Y,'Xedg',xedg,'Yedg',yedg,'Xscale',xscale,'Yscale',yscale)
%
% INPUT:
%   X - X axis data
%   Y - Y axis data
% Optional 'name',value pairs:
%   'Xedg' - X bin edges vector or scalar number of bins, default = 128
%   'Yedg' - Y bin edges vector or scalar number of bins, default = 128
%   'Xscale' - 'lin' or 'log' to indicate axis scaling, default = 'lin'
%   'Yscale' - 'lin' or 'log' to indicate axis scaling, default = 'lin'
%   'msize' - marker size as scalar or vector for each X,Y pair, default = 5
%   'parent' - parent axes handle
%
% OUTPUT: 
%   h - scattergroup handle
%   H - bivariate histogram
%   Z - histogram bin density for each X Y pair
%
% see also: kde2, bivHist, histcounts2

% parse input
for m=1:2:length(varargin)
    switch lower(varargin{m})
        case 'xedg'
            Xedg = varargin{m+1};
        case 'yedg'
            Yedg = varargin{m+1};
        case 'xscale'
            xscale = varargin{m+1};
        case 'yscale'
            yscale = varargin{m+1};
        case 'msize'
            msz = varargin{m+1};
        case 'parent'
            ax = varargin{m+1};
        otherwise 
            error('dplot.m input error: Unexpected input')
    end
end


% for m=1:length(varargin)
%     if ischar(varargin{m})
%         scale=varargin{m};
%     elseif isvector(varargin{m})
%         edg{m}=varargin{m};
%     end
% end
% 
% if exist('edg','var')
%     Xedg=edg{1}; 
%     if length(edg)>1
%         Yedg=edg{2};
%     end
% end
if ~exist('xscale','var'); xscale = 'lin'; end
if ~exist('yscale','var'); yscale = 'lin'; end
if ~exist('msz','var'); msz = 5; end
if ~exist('ax','var'); ax = gca; end
if ~exist('X','var') || isempty(X) || ~exist('Y','var') || isempty(Y)
    h = [];
    H = [];
    Z = [];
    return; 
end

if exist('Xedg','var') && ~exist('Yedg','var')
    Yedg = Xedg;
end

if exist('Xedg','var') && isscalar(Xedg)
    nXbins = Xedg;
else
    nXbins = 128;
end

if exist('Yedg','var') && isscalar(Yedg)
    nYbins = Yedg;
else 
    nYbins = 128;
end

% nz=all([X Y]>0,2); % find all non-zero
% X=X(nz); 
% Y=Y(nz);

% calculate bin vectors from scalar or missing input
if ~exist('Xedg','var') || isscalar(Xedg) || isempty(Xedg)
    if strcmp(xscale,'log')
%         Xedg=logspace(floor(log10(min(X+1))),log10(max(X+1)),nXbins);
        Xedg = min(X)-1+logspace(0,log10(max(X)-min(X)+1),nXbins);
        Xedg(end) = max(X);
    elseif strcmp(xscale,'lin')
%         Xedg=linspace(min(X+1),max(X+1),nXbins);
        Xedg = linspace(min(X),max(X),nXbins);
    end
end
if ~exist('Yedg','var') || isscalar(Yedg) || isempty(Yedg)
    if strcmp(yscale,'log')
%         Yedg=logspace(floor(log10(min(Y+1))),log10(max(Y+1)),nYbins);
        Yedg = min(Y)-1+logspace(0,log10(max(Y)-min(Y)+1),nYbins);
        Yedg(end) = max(Y);
    elseif strcmp(yscale,'lin')
%         Yedg=linspace(min(Y+1),max(Y+1),nYbins);
        Yedg = linspace(min(Y),max(Y),nYbins);
    end
end
% [H,Z]=bivHist(X+1,Y+1,Xedg,Yedg);
[H,~,~,bx,by] = histcounts2(X,Y,Xedg,Yedg);
% Z = H(sub2ind(size(H),bx,by));
% fnt = isfinite(X) & isfinite(Y);
fnt = isfinite(X) & isfinite(Y) & bx>0 & by>0; % modified from above 2026-07-04, MMcFarland
Z1 = H(sub2ind(size(H),bx(fnt),by(fnt)));
Z = NaN(size(X));
Z(fnt) = Z1;

%figure
h = scatter(ax,X,Y,msz,log10(Z+1),'filled');

set(gca,'xscale',xscale);
set(gca,'yscale',yscale);

grid on
% colormap(flipud(hot)) % good for black backgrounds

