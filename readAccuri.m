function [dat,vol,hdr]=readAccuri(fin,p)

% Read in data from Accuri C6 flow cytometer
% Concatenates multiple files
% Requires readfcs.m
%
% USAGE: [dat,vol,hdr] = readAccuri(fin,p)
%
% INPUT:
%   fin - input file name (string) or directory structure
%   p - optional parameter type indicator: 'H' for height or 'A' for area
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

dat=[];
vol=0;
for n=1:length(d)
    [dat1,hdr,fcshdr] = readfcs([d(n).folder filesep d(n).name]);
    % dat1 = dat1(:,cid);
    dat = [dat; dat1];
    vol = vol + fcshdr.vol;
end
dat(:,1:12) = log10(dat(:,1:12)+1); % log convert
dat = dat(:,cid); % select parameter type
hdr = hdr(cid);

dat = array2table(dat,'VariableNames',matlab.lang.makeValidName(hdr)); % make table
