function key = dir2key(d,prj,f)

% Generate sample ID from directory structure
% to use as key variable for merging tables.
% File names in d must contain site and bottle IDs
% separated by '_' characters. Folder names must
% contain sample dates in 'yyyy-mm-dd' format.
%
% USAGE:
%   key = dir2key(d,prj,f)
%
% INPUT:
%   d - directory structure containing the files of interest
%   prj - (char or string) project name prefix
%   f - (optional vector) field indeces for site, bottle, and other
%       ID info, default = [2 3]
%
% OUTPUT:
%   key - string vector of sampleIDs for each file in d

prj = string(prj);
if ~exist('f','var') || isempty(f); f = [2 3]; end

% d = dir('C:\Users\mmcfarland\Documents\PROJECTS\NEP\FCdata\**\*.fcs');
% idx = contains({d.name},{'_son',' Pre',' Post',' Dec','clog'},IgnoreCase=true); % ignore these files
% d = d(~idx);

sb = regexp({d.name},'\s0?([A-Za-z0-9]+)_?0?([A-Za-z0-9]*)_?0?([A-Za-z0-9]*)_?0?([A-Za-z0-9]*)','tokens'); % extract site and bottle IDs as tokens
sb = vertcat(sb{:});
sb = vertcat(sb{:});
sb(contains(sb,'son')) = {''}; % get rid of son/unson
sb = string(sb);

sd = regexp({d.folder},'_?(\d{4})-?(\d{2})-?(\d{2})','tokens'); % get sample date from folder names
sd = vertcat(sd{:});
sd = vertcat(sd{:});
sd = strcat(sd(:,1),sd(:,2),sd(:,3)); % concatenate year month day
sd = string(sd);

key = prj + "_" + sd + "_" + sb(:,f(1)); % concatenate project, date, and siteID

for m = 2:length(f)
    i = sb(:,f(m))~=""; % find bottleIDs
    key(i) = key(i) + "_" + sb(i,f(m)); % concatenate bottleIDs
end
