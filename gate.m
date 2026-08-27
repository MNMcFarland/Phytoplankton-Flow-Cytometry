function [idx,lbv] = gate(dat,G)

% Gate data set hierarchically using gates in G
% Returns gate membership for all observations in dat
%
% USAGE: [idx,lbv] = gate(dat,G)
%
% INPUT:
%   dat - flow cyt data set: events x parameters
%   G - gates structure (from polygate.m) with fields: 
%       cols - columns of dat  
%       x, y - gate vertices (vectors). Vertices must be ordered.
%       name - gate name ('string') 
%       parent - parent gate number, events within gate must also be inside
%           parent gate. Parent gates must come before children in gates
%           structure.
%       color - 3 element vector of RGB values (0-1)
%
% OUTPUT: 
%   idx - logical index array for each gate in G 
%   lbv - label vector showing last child gate membership

% Malcolm McFarland, version 2021
% Harbor Branch Oceanographic Institute, Florida Atlantic University
% mmcfarland@fau.edu

idx = false(size(dat,1),length(G));
for n = 1:length(G)
    if ~isempty(G(n).x)
        % ind = inpolygon(dat(:,G(n).cols(1)), dat(:,G(n).cols(2)), G(n).x, G(n).y); % determine events inside gate
        ind = inpolygon(dat.(G(n).cols{1}), dat.(G(n).cols{2}), G(n).x, G(n).y); % determine events inside gate
        if sum(ind)>0
            if G(n).parent
                if G(n).parent >= n
                    error('gate: parents must precede children in gates structure')
                else
                    ind = ind & idx(:,G(n).parent); % restrict to inside parent gate
                end
            end
        else
            ind = 0;
        end
        idx(:,n) = ind; % add gated data indeces to output
    end
end
% convert logical index matrix to label vector
[r,c] = ind2sub(size(idx),find(idx));
lbv = zeros(size(idx,1),1);
for p=unique(c)'; lbv(r(c==p)) = p; end
