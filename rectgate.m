function [G,idx,h] = rectgate(dat,cols,varargin)

% Interactively draw rectangular gate on plot and return vertices.
% Returns gate vertex coordinates in gates structure G.
% Appends to existing gates structure if provided.
% Requires image processing toolbox, gateplot.m, and dplot.m
%
% USAGE: [G,ind,h] = rectgate(dat,cols,'property',value,...)
%
%   [G,ind,h] = rectgate(dat,cols,'property1',value1,'property2',value2,...) 
%   Generates a new plot from tabular data in 'dat'. 
%   Returns indices of gated data in 'ind'. 
%   Also returns handle of gate object in 'h'. 
%
% INPUT:
%   dat - data matrix (events x parameters)
%   cols - columns of data to gate, 2 element vector
%
%   Property - Value pair arguments (optional)
%   ----------------------------------------------------------------
%       'gates'  - gate structure with fields: cols, x, y, parent, name. 
%               Parent gates must precede children in gates structure. 
%       'lbl'    - x and y axes labels (2 element cell array of char vectors)
%       'number' - gate number (scalar)
%       'parent' - parent gate number (scalar), events in the gated region must
%               also be inside the parent gate
%       'name'   - gate name (string)
%       'color' - gate color (3 element row vector)
%
% OUTPUT:
%   G - gates structure with fields: cols, x, y, parent, name, xlbl, and ylbl for each gate
%   ind - logical index array showing membership for each gate in G (events x gates)
%   h - handle of gate polygon object in plot

% Malcolm McFarland, version 2020
% Harbor Branch Oceanographic Institute, Florida Atlantic University
% mmcfarland@fau.edu

lbl = {['parameter ' num2str(cols(1))],['parameter '  num2str(cols(2))]}; 
par = [];

% parse input
if nargin/2~=round(nargin/2)
    error('polygate: unpaired inputs')
end
for m=1:2:length(varargin)
    if ischar(varargin{m})
        switch lower(varargin{m})
            case 'color'
                c = varargin{m+1};
            case 'lbl'
                lbl = varargin{m+1};
            case 'gates'
                G = varargin{m+1};
            case 'number'
                n = varargin{m+1};
            case 'parent'
                par = varargin{m+1};
            case 'name'
                nm = varargin{m+1};
            otherwise 
                error('polygate.m input error: Unexpected input property!')
        end
    end
end

% set defaults
if ~exist('G','var') || isempty(G)
    if ~exist('nm','var'); nm = 'gate 1'; end
    G = struct('cols',cols,'x',[],'y',[],'name',nm,'parent',[],'color',[0 0 0]);
    if ~exist('n','var'); n = 1; end
end
if ~exist('n','var'); n = length(G)+1; end % increment gate number if not provided
if ~exist('nm','var'); nm = ['gate ' num2str(n)]; end
if ~exist('c','var'); c = lines(n); end
% if ~exist('c','var'); cmap = lines; c = cmap(n,:); end

%% plot all data and existing gates
figure
gateplot(dat,cols,G,'xscale','lin','yscale','lin','lbl',lbl);
% set(gca,'xlim',[min(dat(:,cols(1)))*.9 max(dat(:,cols(1)))*1.1],...
%     'ylim',[min(dat(:,cols(2)))*.9 max(dat(:,cols(2)))*1.1])

% set up or append new gate
G(n).cols = cols;
G(n).parent = par;
G(n).name = nm;
G(n).xlbl = lbl{1};
G(n).ylbl = lbl{2};
G(n).color = c;

if G(n).parent >= n
    error('polygate: parents must precede children in gates structure')
end

% get new gate
% h = drawpolygon('FaceAlpha',0,'FaceSelectable',false,'color',c,'interactionsallowed','none'); 
h = drawrectangle('FaceAlpha',0,'FaceSelectable',false,'color',c,'interactionsallowed','none'); 
% h.InteractionsAllowed = 'none';
% G(n).x = h.Position(:,1);
% G(n).y = h.Position(:,2);
G(n).x = [h.Position(1); h.Position(1); h.Position(1)+h.Position(3); h.Position(1)+h.Position(3)];
G(n).y = [h.Position(2); h.Position(2)+h.Position(4); h.Position(2)+h.Position(4); h.Position(2)];

% plot gated data and get label vector
idx = gateplot(dat,cols,G,'lbl',lbl);

% lgd = findobj(gcf,'type','legend');
% if isempty(lgd)
%     legend(lh,'location','northwest','interpreter','none')
% end


