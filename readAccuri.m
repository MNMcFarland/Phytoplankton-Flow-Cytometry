function [dat,vol,hdr]=readAccuri(fin,p,lc)

% Read in data from Accuri C6 flow cytometer
% Concatenates multiple files
% Requires readfcs.m
%
% USAGE: [dat,vol,hdr] = readAccuri(fin,p,lc)
%
% INPUT:
%   fin - input file name (string) or directory structure
%   p - optional parameter type indicator: 'H' for height or 'A' for area
%   lc - optional boolean flag to log convert data, default = true
%
% OUTPUT:
%   dat - listmode data
%   vol - volume analyzed
%   hdr - dat column headers

if ischar(fin)
    d=dir(fin);
elseif isstruct(fin)
    d=fin;
end

if ~exist('p','var') || isempty(p)
    cid=1:14;
elseif p=='H'
    cid=[7 8 9 10 11 12]; % Height parameters
elseif p=='A'
    cid=[1 2 3 4 5 6]; % Area parameters
end

if ~exist('lc','var') || isempty(lc)
    lc = true;
end

dat=[];
vol=0;
for n=1:length(d)
    [dat1,hdr,fcshdr] = readfcs([d(n).folder filesep d(n).name]);
    % dat1 = dat1(:,cid);
    dat = [dat; dat1];
    vol = vol + fcshdr.vol;
end

if lc
    dat(:,1:12) = log10(dat(:,1:12)+1); % log convert
end

dat = dat(:,cid); % select parameter type
hdr = hdr(cid);

hdr = matlab.lang.makeValidName(hdr);
dat = array2table(dat,'VariableNames',hdr); % make table
