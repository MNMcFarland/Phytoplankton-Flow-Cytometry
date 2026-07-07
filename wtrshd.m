function [lbl,D,mu,sigma,p] = wtrshd(X,r,mn,t)

% Watershed algorithm to cluster unorganized observations.
% Uses a range search to determine density as number of neighbors.
%
% USAGE: [lbl,D,mu,sigma,p] = wtrshd(X,r,mn,t)
%
% EXAMPLE:
%   [lbl,D,mu,sigma,p] = wtrshd(X,r);
%   gm = gmdistribution(mu,sigma,p/sum(p)); % create gaussian mixture model from results
%
% INPUT:
%   X - data matrix, observations (rows) by parameters (columns)
%   r - distance to neighbors for density determination (i.e. neighborhood radius),
%       determines the amount of smoothing of the distribution, default = 0.05*max(X(:))
%   mn - (optional) minimum number of neighbors to start a new cluster, default = 3
%   t - (optional) merge unlabeled observations with D>=t with nearest cluster
%
% OUTPUT:
%   lbl - watershed labels for each observation in X
%   D - density for every observation (number of neighbors within distance r)
%   mu - cluster means
%   sigma - covariance matrix for each cluster
%   p - number of observations in each cluster (mixing proportions = p/sum(p))

%   gm - Gaussian mixture model for all clusters

% Malcolm McFarland 2026

% this can quickly run out of memory
if ~exist('r','var') || isempty(r); r = max(X(:))*0.05; end
if ~exist('mn','var') || isempty(mn); mn = 3; end
% if ~exist('t','var'); t = mn; end

%% use KD tree range search to determine density and neighbors
N = rangesearch(X, X, r); % ,'BucketSize',50); % this might leak memory!
% N = rangesearch(X, X, r, 'Distance','mahalanobis'); % might work better but much slower, use with larger radius
D = cellfun(@length, N);
N = cellfun(@(x) x(1:min([length(x) 1000])), N, 'uniformoutput', false); % only store closest 1000 neighbors

%% watershed algorithm
[~,i] = sort(-D); % sort in order of decreasing density
[~,j] = sort(i);
Ns = N(i,:);

% sz = round(logspace(log10(20),3,length(Ns))); % variable max neighborhood size small -> large
lbl = zeros(size(Ns));
for m = 1:length(Ns)
    n = Ns{m}; % original indices of neighbors within range r
%     n = Ns{m}(1:min([length(Ns{m}) sz(m)])); % original indices of neighbors within range r, limited by sz
    n = j(n); % indeces of sorted neighbors within range r
    nlb = unique(lbl(n)); % unique labels within this neighborhood
    if length(nlb)==1 && nlb==0 && length(n)>=mn % if no labeled neighbors
        lbl(m) = max(lbl)+1; % start new cluster
    elseif length(nlb)==2 % if only one labeled neighbor
        lbl(m) = nlb(2); % propagate cluster label
    elseif length(nlb)>2 % if multiple different labeled neighbors
        d = find(lbl(n),1,'first'); % find index of closest labeled neighbor  
        lbl(m) = lbl(n(d)); % choose label of closest neighbor
    end
end
lbl = lbl(j); % undo density sort

%% assign unlabelled observations greater than t
nd = size(X,2);
k = max(lbl);
if exist('t','var') && ~isempty(t)
    wi = lbl==0 & D>=t; % indices of unlabeled observations
    dst = zeros(sum(wi),k);
    for m = 1:k % distance to each centroid for unlabeled observations
        if sum(lbl==m) > nd
            dst(:,m) = mahal(X(wi,:),X(lbl==m,:)); % squared mahalanobis distance
        else
            muk = repmat(mean(X(lbl==m,:),1),sum(wi),1);
            dst(:,m) = sum((X(wi,:)-muk).^2,2); % squared euclidean distance
        end
    end
    [~,d] = min(dst,[],2);
    lbl(wi) = d;
end

%% compute cluster means (mu), covariance (sigma), and counts (p)
if nargout>2
    mu = zeros(k,nd);
    p = zeros(k,1);
    sigma = zeros(nd,nd,k);
    for m = 1:k
        mu(m,:) = mean(X(lbl==m,:),1);
        p(m) = sum(lbl==m);
        if p(m)>1
            sigma(:,:,m) = cov(X(lbl==m,:));
        else
            sigma(:,:,m) = eye(nd); % zeros(nd);
        end
    end
%     gm = gmdistribution(mu,sigma,p/sum(p)); % create gaussian mixture model
    % NlogL=-sum(log(pdf(gm,X(lbl>0,:)))); % negative log Likelihood
    % aic=2*NlogL+2*k; % Akaike's Information Criterion
    % aic=2*NlogL+2*k*(n/(n-K-1)); % Corrected Akaike's Information Criterion, n = sample size
end

% % zero clusters with < mn members
% idx = p<mn;
% lbl(ismember(lbl,find(idx))) = 0;
% mu(idx,:) = 0;
% sigma(:,:,idx) = 0;
% p(idx) = 0;

