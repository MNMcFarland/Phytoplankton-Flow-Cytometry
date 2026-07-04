function [X,Y,Xv,Yv] = maketrn(fnames,ws,ID,ord)

% Make ML training data set from .fcs files
% labels for data are in ws.lbl
% Randomly assigns 20% of data to the validation set
%
% USAGE: [X,Y,Xv,Yv] = maketrn(fnames,ws,ord)
%
% INPUT:
%   fnames - cell array of file names
%   ws - output from FCwtrshd
%   ID - cell array of cluster IDs
%   ord - (optional) columns to include in output
%
% OUTPUT:
%   X - training data set values
%   Y - training data set labels
%   Xv - validation data set values
%   Yv - validation data set labels

if ~exist('ord','var') || isempty(ord); ord = [7 11 10 12 8 9]; end

X = []; Y = [];
for m = 1:length(fnames)
    n = strcmp({ws.name},fnames{m});
%     lbl = ws(n).lbl;
    nzl = ws(n).lbl>0; % non-zero labelled events
    Y = [Y; ws(n).lbl(nzl)]; % concatenate and remove zero labels
    
    fcdat = readAccuri(fnames{m});
    if istable(fcdat); fcdat = table2array(fcdat); end
    X1 = log10(fcdat(:,ord)+1);
    X = [X; X1(nzl,:)]; % concatenate and remove zero labels
end
clear X1

X = X/7.2247; % normalize training data

% % remove less abundant classes
% unq = unique(Y);
% cnt = zeros(size(unq));
% for m = 1:length(unq)
%     cnt(m) = sum(Y==unq(m));
% end
% [cnts,srti] = sort(cnt,'descend');
% plot(cumsum(cnts)/sum(cnts)); % up to 30 >= 99% 
% unqs = unq(srti);
% incl = unqs(1:30); % class labels to include
% linc = ismember(Y, incl);
% Y = Y(linc);
% X = X(linc,:);

% re-label according to cluster IDs
Y1 = zeros(size(Y));
for m = 1:length(ID)
    Y1(ismember(Y,ID{m})) = m;
end
% Y1(ismember(Y,bkg)) = 1;
% Y1(ismember(Y,ded)) = 2;
% Y1(ismember(Y,syn)) = 3;
% Y1(ismember(Y,pek)) = 4;
% Y1(ismember(Y,nek1)) = 5;
% Y1(ismember(Y,nek2)) = 6;
% Y1(ismember(Y,cyb)) = 7;
% Y1(ismember(Y,crp)) = 8;
% Y1(ismember(Y,unk)) = 9;
% Y1(~ismember(Y,[syn pek nek cyb crp unk])) = 7; % other
Y = Y1;
clear Y1

% split into training and validation data sets
n = size(X,1);
% p = randperm(n,5e4);
p = randperm(n,round(0.2*n));
Xv = X(p,:);
X(p,:) = [];
Yv = Y(p);
Y(p) = [];



