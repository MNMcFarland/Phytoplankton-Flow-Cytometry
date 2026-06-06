function cmap=lbl2cmap(lbl,mu,mx)

% compute cluster color map from 3 parameter means
%
% USAGE: cmap=lbl2cmap(lbl,mu,mx)
%
% INPUT:
%   lbl - label vector from cluster operation (integers)
%       values should correspond to rows of mu
%   mu - means for each label ID, cluster x 3 parameters
%   mx - (optional) maximum for each parameter (vector) or all (scalar)
%
% OUTPUT:
%   cmap - rgb color map for each lbl

if size(mu,2)~=3
    error('mu must have 3 columns')
end

if ~exist('mx','var')
    mx=max(mu);
end

u=unique(lbl);
if size(u,1)>size(u,2)
    u=u';
end
u(u<1)=[]; % ignore labels<1

colors=bsxfun(@rdivide,mu,mx);
cmap=zeros(length(lbl),3);
for m=u
    n=sum(lbl==m);
    cmap(lbl==m,:)=repmat(colors(m,:),n,1);
end

