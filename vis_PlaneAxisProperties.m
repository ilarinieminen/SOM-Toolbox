function vis_PlaneAxisProperties(ax,lattice,msize,pos)

% VIS_PLANEAXISPROPERTIES Set axis properties for SOM_CPLANE, 
%                         SOM_PIEPLANE, SOM_BARPLANE and SOM_PLOTPLANE.
%
% vis_PlaneAxisProperties(ax,lattice,msize,pos)
%
%  Input arguments: 
%   ax        (scalar) axis handle     
%   lattice   (string) 'hexa', 'rect', 'hexaU' or 'rectU'
%             (matrix) defines the patch, see e.g. help vis_patch
%   msize     (vector) a 1x2 vector defining the grid size
%   pos       (vector) a 1x2 vector that determines position of
%                      origin or NaN which means default operation: 
%                      origin to [1 1] and tighten axis limits 
%                      according to the grid size.
% 
% This is a subfunction for SOM_CPLANE, SOM_PIEPLANE, SOM_BARPLANE and
% SOM_PLOTPLANE. This subfunction sets the proper values for axis.

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 060799

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xdim=msize(1);ydim=msize(2);
set(ax,'Visible','off');
set(get(ax,'Title'),'Visible','on');
set(ax,'XaxisLocation','Top');            % axis orientation
set(ax,'xdir','normal');                  % = axis ij = matrix mode
set(ax,'ydir','reverse'); 

switch lattice
case {'rect', 'rectU'}
  lelim=-.51; rilim=.51; uplim=-.51; lolim=.51;  % axis limits
  set(ax,'DataAspectRatio', [1 1 1]);            % =axis equal
case {'hexa','hexaU'}
  lelim=-.51; rilim=1.01; uplim=-.67; lolim=.67; % axis limits
  set(ax,'DataAspectRatio',[0.9015 1 1]);        % this corrects hexagons
end

% Nan: default origin [1 1] & tighten the axis
if isnan(pos)
  set(ax,'XLim',[1+lelim ydim+rilim],'YLim',[1+uplim xdim+lolim], ...
      'XLimMode','manual','YLimMode','manual'); % tighten the axis
end

