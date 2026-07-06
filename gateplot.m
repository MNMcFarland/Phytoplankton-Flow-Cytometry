function [idx, lbv] = gateplot(dat,cols,G,varargin)

% Plot gated data set.
% Use with log transformed flow cytometry data
% Requires dplot.m
%
% USAGE: [idx,lbv] = gateplot(dat,cols,G,'property',value,...)
%
% INPUT:
%   dat - flow cyt data set: events x parameters
%   cols - columns of dat to plot (2 element vector)
%   G - gates structure (from polygate.m) with fields: 
%       cols - columns of dat  
%       x, y - gate vertices (vectors). Vertices must be ordered.
%       name - gate name ('string') 
%       parent - parent gate number, events within gate must also be inside
%           parent gate. Parent gates must come before children in gates
%           structure.
%       color - 3 element vector of RGB values (0-1)
%
%   Property, Value - paired arguments (optional)
%      'xscale' - 'log' or 'lin' x axis scale, default = 'lin'
%      'yscale' - 'log' or 'lin' y axis scale, default = 'lin'
%      'lbl'    - x and y axes labels (2 element cell array of strings)
%
% OUTPUT: 
%   plot
%   idx - logical index array for each gate in G 
%   lbv - label vector showing last gate membership

% Malcolm McFarland, version 2020
% Harbor Branch Oceanographic Institute, Florida Atlantic University
% mmcfarland@fau.edu

% parse optional input
if length(varargin)/2 ~= round(length(varargin)/2)
    error('gateplot: unpaired property - value inputs')
end
for m=1:2:length(varargin)
    if ischar(varargin{m})
        switch lower(varargin{m})
            case 'xscale'
                xscale = varargin{m+1};
            case 'yscale'
                yscale = varargin{m+1};
            case 'lbl'
                lbl = varargin{m+1};
                ovrd = false; % don't override specified labels
            otherwise 
                error('gateplot.m input error: Unexpected input property!')
        end
    end
end

% set defaults
if ~exist('xscale','var') || isempty(xscale); xscale='lin'; end
if ~exist('yscale','var') || isempty(yscale); yscale='lin'; end
if ~exist('lbl','var') || isempty(lbl)
    lbl = {['parameter ' num2str(cols(1))], ['parameter ' num2str(cols(2))]}; 
    ovrd = true; % override these auto generated labels if possible
end

%%
% [~,D0,~,h] = kde2(dat(:,cols(1)),dat(:,cols(2)),'sd',1,'sz',256,'p',2);
[h,~,D0] = dplot(dat(:,cols(1)),dat(:,cols(2)),'xscale',xscale,'yscale',yscale); % density plot all the data
mn1 = min(dat(:,cols(1)));
mx1 = max(dat(:,cols(1)));
pd1 = 0.1*(mx1 - mn1); % pad from data range for x axis
mn2 = min(dat(:,cols(2)));
mx2 = max(dat(:,cols(2)));
pd2 = 0.1*(mx2 - mn2); % pad from data range for y axis
set(gca,'xlim',[mn1-pd1 mx1+pd1],'ylim',[mn2-pd2 mx2+pd2])
% set(gca,'xlim',[min(dat(:,cols(1)))*.9 max(dat(:,cols(1)))*1.1],...
%     'ylim',[min(dat(:,cols(2)))*.9 max(dat(:,cols(2)))*1.1]) % this scaling method fails for negative data
% set(gca,'xlim',[floor(min(dat(:,cols(1)))) ceil(max(dat(:,cols(1))))],...
%     'ylim',[floor(min(dat(:,cols(2)))) ceil(max(dat(:,cols(2))))])
% colormap(bone) % use neutral colormap
D0 = log10(D0+1);
h.CData = [1 1 1] .* D0/max(D0);
hold on
% if exist('lbl','var') && ~isempty(lbl) && iscell(lbl)
    xlabel(lbl{1})
    ylabel(lbl{2})
% end

% cmap = lines(length(G));
idx = false(size(dat,1),length(G));
% idx = zeros(size(dat,1),1);
h = zeros(length(G),1);
for n = 1:length(G)
    if ~isempty(G(n).x)
        ind = inpolygon(dat(:,G(n).cols(1)), dat(:,G(n).cols(2)), G(n).x, G(n).y); % determine events inside gate
        if sum(ind)>0
            if G(n).parent
                if G(n).parent >= n
                    error('gateplot: parents must precede children in gates structure')
                else
                    ind = ind & idx(:,G(n).parent); % restrict to inside parent gate
                end
            end
%             [~,D,~,h(n)] = kde2(dat(ind,cols(1)), dat(ind,cols(2)),'sd',1,'sz',256,'p',2); 
            [h1,~,D] = dplot(dat(ind,cols(1)), dat(ind,cols(2)),'xscale',xscale,'yscale',yscale);
            if ~isempty(D)
                h(n) = h1;
                D = log10(D+1);
                c = G(n).color;
                set(h(n), 'CData', c + (1-c) .* D/max(D0) );
            end
%             set(h(n), 'CData', c + (1-c) .* log10(D+1)/max(D0) );
%             set(h(n), 'CData', cmap(n,:) + (1-cmap(n,:)) .* D/max(D) );
%             h(n) = scatter(dat(ind,cols(1)), dat(ind,cols(2)), 5, cmap(n,:), 'filled'); % plot events inside this gate
        else
            ind = 0;
        end
        if G(n).cols(1)==cols(1) && G(n).cols(2)==cols(2) % if gate columns match plot axes
            line([G(n).x; G(n).x(1)],[G(n).y; G(n).y(1)],'color',G(n).color,'linewidth',1.5); % draw gate as line
            if ovrd && isfield(G,'xlbl') && ~isempty(G(n).xlbl); xlabel(G(n).xlbl); end
            if ovrd && isfield(G,'ylbl') && ~isempty(G(n).ylbl); ylabel(G(n).ylbl); end
        end
        idx(:,n) = ind; % add gated data indeces to output
%         idx(ind) = n; % add gated data indeces to output
    end
end
hold off
if any(h>0)
    legend(h(h>0),{G(h>0).name},'location','northwest','interpreter','none');
end

% convert logical index matrix to label vector
[r,c] = ind2sub(size(idx),find(idx));
lbv = zeros(size(idx,1),1);
for p=unique(c)'; lbv(r(c==p)) = p; end


