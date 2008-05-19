function S=vis_grid_struct

% VIS_GRID_STRUCT Initiates and returns a default grid struct
%
%  S = vis_grid_struct
%
%  Output arguments:
%   S       (struct) som_grid struct 
%
% Fields and their initial values: 
%  S.type='som_grid';
%  S.neigh=neigh;
%  S.shape='sheet';
%  S.msize=msize;
%  S.coord=[];
%  S.line='-';
%  S.linecolor=[.8 .8 .8];
%  S.linewidth=0.5;
%  S.marker='o';
%  S.markersize=6;
%  S.markercolor='k';
%  S.surf=[];
%  S.label=[];
%  S.labelcolor='g';
%  S.labelsize=12;
%
% This is a subfunction for SOM_GRID. 

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 090499

S.type='som_grid';
S.neigh='hexa';
S.shape='sheet';
S.msize=[1 1];
S.coord=[];
S.line='-';
S.linecolor=[.8 .8 .8];
S.linewidth=0.5;
S.marker='o';
S.markersize=6;
S.markercolor='k';
S.surf=[];
S.label=[];
S.labelcolor='g';
S.labelsize=12;
