function [in,x,y] = inrect(xq,yq,r)

% determine points located inside or on edge of rectangle region
%
% USAGE: [in,on] = inrect(xq,yq,r)
%
% INPUT:
%   xq - x-coordinates of query points (vector)
%   yq - y-coordinates of query points (vector)
%   r - rectangle object from drawrectangle function
%
% OUTPUT:
%   in - logical vector indicating points in or on edge of rectangle
%   x - rectangle vertices x-coordinates
%   y - rectangle vertices y-coordinates

x = [r.Position(1) r.Position(1)+r.Position(3) r.Position(1)+r.Position(3) r.Position(1)];
y = [r.Position(2) r.Position(2) r.Position(2)+r.Position(4) r.Position(2)+r.Position(4)];
in = inpolygon(xq,yq,x,y);
