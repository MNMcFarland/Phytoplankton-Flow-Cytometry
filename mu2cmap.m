function cmap=mu2cmap(mu,mx,mn)

% compute cluster color map from 3 parameter means
%
% USAGE: cmap=mu2cmap(mu,mx,mn)
%
% INPUT:
%   mu - means for each label ID, cluster x 3 parameters
%   mx - (optional) maximum for each parameter (vector) or all (scalar)
%   mn - (optional) minimum for each parameter (vector) or all (scalar)
%
% OUTPUT:
%   cmap - rgb color map for each row of mu

if size(mu,2)~=3
    error('mu must have 3 columns')
end

if ~exist('mx','var') || isempty(mx); mx=max(mu); end
if ~exist('mn','var') || isempty(mn); mn=min(mu); end

mu=bsxfun(@minus,mu,mn);
cmap=bsxfun(@rdivide,mu,mx-mn);
% cmap=bsxfun(@rdivide,mu,mx);
% cmap=[0 0 0; cmap]; % add zeros to account for unlabeled events
cmap(cmap>1)=1;
cmap(cmap<0)=0;
